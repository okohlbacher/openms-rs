// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded native idXML 1.5 interchange. See `docs/IDXML_SUPPORT.md` for the
//! supported source encodings and explicit errors for unrepresentable state.

use crate::chemistry::{AASequence, ModificationsDB};
use crate::comparison::Tolerance;
use crate::identification::{
    AnalysisResult, EnzymeTermSpecificity, FlankingResidue, PeakAnnotation, PeakMassType,
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification, SearchParameters,
};
use crate::metadata::{CompletionTime, MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use quick_xml::{NsReader, events::Event, name::ResolveResult};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};

const IDENTIFIER: &str = "openms-rust:run_identifier";
const RANK: &str = "openms-rust:rank";

/// Flat native records linked by run identifier, plus otherwise unused search
/// parameter blocks. XML IDs are transport references and are regenerated.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IdXmlDocument {
    pub document_id: String,
    pub protein_identifications: Vec<ProteinIdentification>,
    pub peptide_identifications: Vec<PeptideIdentification>,
    pub unreferenced_search_parameters: Vec<SearchParameters>,
}

#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    pub max_xml_bytes: u64,
    /// Total XML elements, including UserParam and modification entries.
    pub max_records: usize,
    /// Maximum entries in one evidence or metadata list.
    pub max_list_items: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
            max_list_items: 1_000_000,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct WriteOptions {
    pub max_xml_bytes: usize,
    pub max_records: usize,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_records: 1_000_000,
        }
    }
}
fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
fn number<T: std::str::FromStr>(value: &str) -> Result<T> {
    value
        .trim()
        .parse()
        .map_err(|_| bad(format!("invalid number {value:?}")))
}
fn finite(value: &str) -> Result<f64> {
    let value = number::<f64>(value)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("idXML numeric values must be finite"))
    }
}
fn boolean(value: &str) -> Result<bool> {
    match value {
        "true" | "1" => Ok(true),
        "false" | "0" => Ok(false),
        _ => Err(bad("invalid XML boolean")),
    }
}
fn xml_text(value: &str) -> Result<()> {
    if value.chars().all(|c| {
        matches!(c, '\t' | '\n' | '\r')
            || ('\u{20}'..='\u{d7ff}').contains(&c)
            || ('\u{e000}'..='\u{fffd}').contains(&c)
            || c >= '\u{10000}'
    }) {
        Ok(())
    } else {
        Err(bad("invalid XML 1.0 character"))
    }
}
fn escape(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}
fn date(value: &str) -> Result<()> {
    // The native wall-clock validator checks calendar fields; preserve optional
    // xs:dateTime fractional seconds and explicit timezone without conversion.
    if value.len() < 19 || !value.is_ascii() || value.as_bytes()[10] != b'T' {
        return Err(bad("invalid idXML dateTime"));
    }
    value[..19].parse::<CompletionTime>()?;
    let mut suffix = &value[19..];
    if let Some(fraction) = suffix.strip_prefix('.') {
        let count = fraction.bytes().take_while(u8::is_ascii_digit).count();
        if count == 0 {
            return Err(bad("empty fractional seconds"));
        }
        suffix = &fraction[count..];
    }
    if suffix.is_empty() || suffix == "Z" {
        return Ok(());
    }
    let b = suffix.as_bytes();
    if b.len() != 6
        || !matches!(b[0], b'+' | b'-')
        || b[3] != b':'
        || !b[1..3].iter().chain(&b[4..]).all(u8::is_ascii_digit)
    {
        return Err(bad("invalid dateTime timezone"));
    }
    let h: u8 = number(&suffix[1..3])?;
    let m: u8 = number(&suffix[4..])?;
    if h > 14 || m > 59 || (h == 14 && m != 0) {
        return Err(bad("invalid dateTime timezone offset"));
    }
    Ok(())
}

#[derive(Default)]
struct Node {
    name: String,
    attrs: BTreeMap<String, String>,
    children: Vec<Node>,
}
impl Node {
    fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
    fn attr(&mut self, name: &str, value: impl ToString) {
        self.attrs.insert(name.into(), value.to_string());
    }
    fn get(&self, key: &str) -> Result<&str> {
        self.attrs
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| bad(format!("{} requires {key}", self.name)))
    }
    fn optional(&self, key: &str) -> Option<&str> {
        self.attrs.get(key).map(String::as_str)
    }
    fn check(&self, attrs: &[&str], children: &[&str]) -> Result<()> {
        for attr in self.attrs.keys() {
            if !attrs.contains(&attr.as_str()) {
                return Err(unsupported(format!("{} attribute {attr}", self.name)));
            }
        }
        let mut previous = 0;
        for child in &self.children {
            let index = children
                .iter()
                .position(|&name| name == child.name)
                .ok_or_else(|| unsupported(format!("{} inside {}", child.name, self.name)))?;
            if index < previous {
                return Err(bad(format!("invalid child order in {}", self.name)));
            }
            previous = index;
        }
        Ok(())
    }
}

fn declaration(value: &str) -> Result<()> {
    xml_text(value)?;
    let trim = |s: &str| s.trim_start_matches([' ', '\t', '\n', '\r']).len();
    let mut rest = value
        .strip_prefix("xml")
        .ok_or_else(|| bad("invalid XML declaration"))?;
    let mut phase = 0;
    loop {
        let remaining = trim(rest);
        if remaining == rest.len() {
            if rest.is_empty() && phase > 0 {
                break;
            }
            return Err(bad("XML declaration fields require whitespace separators"));
        }
        rest = &rest[rest.len() - remaining..];
        if rest.is_empty() {
            break;
        }
        let end = rest
            .find(['=', ' ', '\t', '\n', '\r'])
            .ok_or_else(|| bad("invalid XML declaration field"))?;
        let key = &rest[..end];
        rest = &rest[end..];
        rest = &rest[rest.len() - trim(rest)..];
        rest = rest
            .strip_prefix('=')
            .ok_or_else(|| bad("missing declaration equals sign"))?;
        rest = &rest[rest.len() - trim(rest)..];
        let quote = rest
            .chars()
            .next()
            .filter(|c| matches!(c, '\'' | '"'))
            .ok_or_else(|| bad("declaration value requires quotes"))?;
        rest = &rest[1..];
        let end = rest
            .find(quote)
            .ok_or_else(|| bad("unterminated declaration value"))?;
        let field = &rest[..end];
        rest = &rest[end + 1..];
        match (phase, key) {
            (0, "version") => {
                if field != "1.0" {
                    return Err(unsupported("XML version other than 1.0"));
                }
                phase = 1;
            }
            (1, "encoding") => {
                if !field.eq_ignore_ascii_case("UTF-8") {
                    return Err(unsupported("idXML reader requires UTF-8"));
                }
                phase = 2;
            }
            (1 | 2, "standalone") => {
                if !matches!(field, "yes" | "no") {
                    return Err(bad("invalid standalone declaration value"));
                }
                phase = 3;
            }
            _ => return Err(bad("unknown, duplicate or misplaced XML declaration field")),
        }
    }
    if phase == 0 {
        return Err(bad("missing XML declaration version"));
    }
    Ok(())
}

fn parse_xml(input: impl BufRead, options: &ReadOptions) -> Result<Node> {
    if options.max_records == 0 || options.max_list_items == 0 {
        return Err(bad("idXML limits must be positive"));
    }
    let limit = options
        .max_xml_bytes
        .checked_add(1)
        .ok_or_else(|| bad("invalid XML byte limit"))?;
    let mut reader = NsReader::from_reader(input.take(limit));
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    let mut stack: Vec<Node> = Vec::new();
    let mut root = None;
    let mut count = 0usize;
    let mut declared = false;
    let mut at_start = true;
    loop {
        let decoder = reader.decoder();
        let (ns, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| bad(e.to_string()))?;
        let unbound = matches!(ns, ResolveResult::Unbound);
        if reader.buffer_position() > options.max_xml_bytes {
            return Err(bad("idXML exceeds XML byte limit"));
        }
        match event {
            Event::Start(start) => {
                if !unbound {
                    return Err(unsupported("idXML elements must have no namespace"));
                }
                if stack.is_empty() && root.is_some() {
                    return Err(bad("multiple XML roots"));
                }
                count += 1;
                if count > options.max_records || stack.len() >= 8 {
                    return Err(bad("idXML element/depth limit exceeded"));
                }
                let name = std::str::from_utf8(start.name().as_ref())
                    .map_err(|e| bad(e.to_string()))?
                    .to_owned();
                let mut node = Node::new(&name);
                for attr in start.attributes() {
                    let attr = attr.map_err(|e| bad(e.to_string()))?;
                    let key = std::str::from_utf8(attr.key.as_ref())
                        .map_err(|e| bad(e.to_string()))?
                        .to_owned();
                    if attr.value.contains(&b'<') {
                        return Err(bad("unescaped < in XML attribute"));
                    }
                    let raw = decoder
                        .decode(&attr.value)
                        .map_err(|e| bad(e.to_string()))?;
                    // XML 1.0 normalizes literal whitespace before resolving
                    // character references; &#10; must remain a newline.
                    let normalized = raw.replace("\r\n", "\n").replace(['\r', '\n', '\t'], " ");
                    let value = quick_xml::escape::unescape(&normalized)
                        .map_err(|e| bad(e.to_string()))?
                        .into_owned();
                    xml_text(&value)?;
                    if node.attrs.insert(key, value).is_some() {
                        return Err(bad("duplicate XML attribute"));
                    }
                }
                stack.push(node);
            }
            Event::End(_) => {
                let node = stack.pop().ok_or_else(|| bad("unexpected XML end"))?;
                if let Some(parent) = stack.last_mut() {
                    parent.children.push(node);
                } else {
                    root = Some(node);
                }
            }
            Event::Text(text) => {
                let text = text.xml_content().map_err(|e| bad(e.to_string()))?;
                if !text.chars().all(|c| matches!(c, ' ' | '\t' | '\n' | '\r')) {
                    return Err(bad("idXML elements cannot contain text"));
                }
            }
            Event::Comment(text) => {
                xml_text(&text.decode().map_err(|e| bad(e.to_string()))?)?;
            }
            Event::PI(text) => {
                xml_text(std::str::from_utf8(text.as_ref()).map_err(|e| bad(e.to_string()))?)?;
                let target = std::str::from_utf8(text.target()).map_err(|e| bad(e.to_string()))?;
                xml_name(target, true)?;
                if target.eq_ignore_ascii_case("xml") {
                    return Err(bad("reserved XML processing-instruction target"));
                }
            }
            Event::Decl(decl) => {
                if !at_start || declared || !stack.is_empty() || root.is_some() {
                    return Err(bad("misplaced XML declaration"));
                }
                declared = true;
                declaration(std::str::from_utf8(decl.as_ref()).map_err(|e| bad(e.to_string()))?)?;
            }
            Event::Eof => break,
            _ => return Err(unsupported("DTD, CDATA and text entity nodes in idXML")),
        }
        at_start = false;
        buffer.clear();
    }
    if !stack.is_empty() {
        return Err(bad("truncated idXML"));
    }
    let root = root.ok_or_else(|| bad("missing IdXML root"))?;
    if root.name != "IdXML" {
        return Err(bad("expected IdXML root"));
    }
    Ok(root)
}

fn list<'a>(value: &'a str, options: &ReadOptions) -> Result<Vec<&'a str>> {
    let inner = value
        .strip_prefix('[')
        .and_then(|v| v.strip_suffix(']'))
        .ok_or_else(|| bad("metadata list requires [...]"))?;
    if inner.is_empty() {
        return Ok(Vec::new());
    }
    let mut values = Vec::new();
    for value in inner.split(',') {
        if values.len() == options.max_list_items {
            return Err(bad("metadata list exceeds limit"));
        }
        values.push(value);
    }
    Ok(values)
}
fn read_meta(node: &Node, options: &ReadOptions) -> Result<MetaInfo> {
    let mut meta = MetaInfo::new();
    for child in node.children.iter().filter(|c| c.name == "UserParam") {
        child.check(&["name", "type", "value"], &[])?;
        let name = child.get("name")?;
        let value = child.get("value")?;
        let value = match child.get("type")? {
            "string" => MetaValue::from(value),
            "int" => MetaValue::from(number::<i64>(value)?),
            "float" => MetaValue::try_from(finite(value)?)?,
            "stringList" => MetaValue::from(
                list(value, options)?
                    .into_iter()
                    .map(|v| v.replace("\\|", ","))
                    .collect::<Vec<_>>(),
            ),
            "intList" => MetaValue::from(
                list(value, options)?
                    .into_iter()
                    .map(|v| number::<i64>(v.trim()))
                    .collect::<Result<Vec<_>>>()?,
            ),
            "floatList" => MetaValue::try_from(
                list(value, options)?
                    .into_iter()
                    .map(|v| finite(v.trim()))
                    .collect::<Result<Vec<_>>>()?,
            )?,
            other => return Err(unsupported(format!("UserParam type {other}"))),
        };
        if meta.insert(name.into(), value).is_some() {
            return Err(bad(format!("duplicate UserParam {name}")));
        }
    }
    Ok(meta)
}
fn take_text(meta: &mut MetaInfo, key: &str) -> Result<Option<String>> {
    meta.remove(key)
        .map(|v| v.as_str().map(str::to_owned))
        .transpose()
}
fn rank(meta: &mut MetaInfo) -> Result<u32> {
    meta.remove(RANK)
        .map(|v| u32::try_from(v.as_i64()?).map_err(|_| bad("invalid rank")))
        .transpose()
        .map(|v| v.unwrap_or(0))
}
fn specificity(value: &str) -> Result<EnzymeTermSpecificity> {
    match value {
        "unknown" => Ok(EnzymeTermSpecificity::Unknown),
        "full" => Ok(EnzymeTermSpecificity::Full),
        "semi" => Ok(EnzymeTermSpecificity::Semi),
        "none" => Ok(EnzymeTermSpecificity::None),
        _ => Err(bad("invalid EnzymeTermSpecificity")),
    }
}
fn read_search(node: &Node, options: &ReadOptions) -> Result<SearchParameters> {
    node.check(
        &[
            "id",
            "db",
            "db_version",
            "taxonomy",
            "mass_type",
            "charges",
            "enzyme",
            "missed_cleavages",
            "precursor_peak_tolerance",
            "precursor_peak_tolerance_ppm",
            "peak_mass_tolerance",
            "peak_mass_tolerance_ppm",
        ],
        &["FixedModification", "VariableModification", "UserParam"],
    )?;
    let tolerance = |key, ppm| -> Result<Tolerance> {
        let value = finite(node.get(key)?)?;
        Ok(if boolean(node.optional(ppm).unwrap_or("false"))? {
            Tolerance::Ppm(value)
        } else {
            Tolerance::Absolute(value)
        })
    };
    let mut value = SearchParameters {
        database: node.get("db")?.into(),
        database_version: node.get("db_version")?.into(),
        taxonomy: node.optional("taxonomy").unwrap_or("").into(),
        charges: node.get("charges")?.into(),
        mass_type: match node.get("mass_type")? {
            "average" => PeakMassType::Average,
            "monoisotopic" => PeakMassType::Monoisotopic,
            _ => return Err(bad("invalid mass_type")),
        },
        digestion_enzyme: node.optional("enzyme").unwrap_or("unknown_enzyme").into(),
        missed_cleavages: number(node.optional("missed_cleavages").unwrap_or("0"))?,
        fragment_tolerance: tolerance("peak_mass_tolerance", "peak_mass_tolerance_ppm")?,
        precursor_tolerance: tolerance("precursor_peak_tolerance", "precursor_peak_tolerance_ppm")?,
        metadata: read_meta(node, options)?,
        ..Default::default()
    };
    if value.metadata.contains_key("modification_definitions") {
        return Err(unsupported(
            "tool-defined modification registration from idXML",
        ));
    }
    if let Some(spec) = take_text(&mut value.metadata, "EnzymeTermSpecificity")? {
        value.enzyme_specificity = specificity(&spec)?;
    }
    for child in node.children.iter().filter(|c| c.name != "UserParam") {
        child.check(&["name"], &[])?;
        let name = child.get("name")?;
        if name.is_empty() {
            return Err(bad("empty modification name"));
        }
        if child.name == "FixedModification" {
            value.fixed_modifications.push(name.into());
        } else {
            value.variable_modifications.push(name.into());
        }
    }
    value.validate()?;
    Ok(value)
}

fn groups(
    meta: &mut MetaInfo,
    prefix: &str,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
) -> Result<Vec<ProteinGroup>> {
    let mut result = Vec::new();
    for index in 0.. {
        let key = format!("{prefix}_{index}");
        let Some(value) = take_text(meta, &key)? else {
            break;
        };
        let mut fields = value.split(',');
        let probability = finite(fields.next().unwrap_or(""))?;
        let mut accessions = Vec::new();
        for id in fields {
            if accessions.len() == options.max_list_items {
                return Err(bad("protein group accession list exceeds limit"));
            }
            accessions.push(
                refs.get(id)
                    .cloned()
                    .ok_or_else(|| bad(format!("unknown protein group reference {id}")))?,
            );
        }
        if accessions.is_empty() {
            return Err(bad("protein group has no protein references"));
        }
        result.push(ProteinGroup {
            probability,
            accessions,
            ..Default::default()
        });
    }
    if meta
        .keys()
        .any(|key| key.starts_with(&format!("{prefix}_")))
    {
        return Err(bad("protein group indices must be contiguous from zero"));
    }
    Ok(result)
}
fn paths(meta: &mut MetaInfo, key: &str) -> Result<Vec<String>> {
    meta.remove(key)
        .map(|v| v.as_string_list().map(<[String]>::to_vec))
        .transpose()
        .map(|v| v.unwrap_or_default())
}
fn read_protein(
    node: &Node,
    mut run: ProteinIdentification,
    refs: &mut BTreeMap<String, String>,
    all_ids: &mut BTreeSet<String>,
    options: &ReadOptions,
) -> Result<ProteinIdentification> {
    node.check(
        &[
            "score_type",
            "higher_score_better",
            "significance_threshold",
        ],
        &["ProteinHit", "UserParam"],
    )?;
    run.score_type = node.get("score_type")?.into();
    run.higher_score_better = boolean(node.get("higher_score_better")?)?;
    run.significance_threshold = finite(node.optional("significance_threshold").unwrap_or("0"))?;
    run.metadata = read_meta(node, options)?;
    for child in node.children.iter().filter(|c| c.name == "ProteinHit") {
        child.check(
            &["id", "accession", "score", "coverage", "sequence"],
            &["UserParam"],
        )?;
        let id = child.get("id")?;
        xml_id(id)?;
        if !all_ids.insert(id.into()) {
            return Err(bad(format!("duplicate XML ID {id}")));
        }
        let accession = child.get("accession")?.to_owned();
        refs.insert(id.into(), accession.clone());
        let mut hit = ProteinHit {
            accession,
            score: finite(child.get("score")?)?,
            sequence: child.optional("sequence").unwrap_or("").into(),
            coverage: child.optional("coverage").map(finite).transpose()?,
            metadata: read_meta(child, options)?,
            ..Default::default()
        };
        // The source has a -1 unknown coverage sentinel; native records use None.
        if hit.coverage == Some(-1.0) {
            hit.coverage = None;
        }
        hit.rank = rank(&mut hit.metadata)?;
        hit.validate()?;
        run.hits.push(hit);
    }
    run.protein_groups = groups(&mut run.metadata, "protein_group", refs, options)?;
    run.indistinguishable_groups = groups(
        &mut run.metadata,
        "indistinguishable_proteins",
        refs,
        options,
    )?;
    run.primary_ms_run_paths = paths(&mut run.metadata, "spectra_data")?;
    run.raw_ms_run_paths = paths(&mut run.metadata, "spectra_data_raw")?;
    if let Some(identifier) = take_text(&mut run.metadata, IDENTIFIER)? {
        run.identifier = identifier;
    }
    run.validate()?;
    Ok(run)
}
fn words<'a>(text: Option<&'a str>, options: &ReadOptions) -> Result<Vec<&'a str>> {
    let mut result = Vec::new();
    for word in text.unwrap_or("").split_ascii_whitespace() {
        if result.len() == options.max_list_items {
            return Err(bad("evidence list exceeds limit"));
        }
        result.push(word);
    }
    Ok(result)
}
fn position(value: &str) -> Result<Option<usize>> {
    let value: i32 = number(value)?;
    if value == -1 {
        Ok(None)
    } else {
        usize::try_from(value)
            .map(Some)
            .map_err(|_| bad("invalid negative evidence position"))
    }
}
fn flank(value: &str) -> Result<FlankingResidue> {
    let mut chars = value.chars();
    let first = chars.next().ok_or_else(|| bad("empty flank"))?;
    if chars.next().is_some() {
        return Err(bad("flanking residue requires one character"));
    }
    FlankingResidue::from_code(first)
}
fn read_evidence(
    node: &Node,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
) -> Result<Vec<PeptideEvidence>> {
    let columns = ["protein_refs", "aa_before", "aa_after", "start", "end"]
        .map(|key| words(node.optional(key), options));
    let [refs_column, before, after, starts, ends] = columns;
    let (refs_column, before, after, starts, ends) =
        (refs_column?, before?, after?, starts?, ends?);
    let len = [
        refs_column.len(),
        before.len(),
        after.len(),
        starts.len(),
        ends.len(),
    ]
    .into_iter()
    .max()
    .unwrap();
    let mut result = vec![PeptideEvidence::default(); len];
    // Short optional columns populate their prefix, matching the source loader.
    for (i, value) in refs_column.into_iter().enumerate() {
        result[i].protein_accession = refs
            .get(value)
            .cloned()
            .ok_or_else(|| bad(format!("unknown protein reference {value}")))?;
    }
    for (i, value) in before.into_iter().enumerate() {
        result[i].aa_before = flank(value)?;
    }
    for (i, value) in after.into_iter().enumerate() {
        result[i].aa_after = flank(value)?;
    }
    for (i, value) in starts.into_iter().enumerate() {
        result[i].start = position(value)?;
    }
    for (i, value) in ends.into_iter().enumerate() {
        result[i].end = position(value)?;
    }
    Ok(result)
}
fn quoted_split(value: &str, separator: char, max_parts: usize) -> Result<Vec<&str>> {
    let (mut quoted, mut escaped, mut start) = (false, false, 0);
    let mut parts = Vec::new();
    for (index, c) in value.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if c == '\\' && quoted {
            escaped = true;
        } else if c == '"' {
            quoted = !quoted;
        } else if c == separator && !quoted {
            if parts.len() >= max_parts.saturating_sub(1) {
                return Err(bad("fragment annotation list exceeds limit"));
            }
            parts.push(&value[start..index]);
            start = index + c.len_utf8();
        }
    }
    if quoted || escaped {
        return Err(bad("unterminated fragment annotation quote"));
    }
    parts.push(&value[start..]);
    Ok(parts)
}
fn annotations(value: &str, options: &ReadOptions) -> Result<Vec<PeakAnnotation>> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    let records = quoted_split(value, '|', options.max_list_items)?;
    if records.len() > options.max_list_items {
        return Err(bad("fragment annotation limit exceeded"));
    }
    records
        .into_iter()
        .map(|record| {
            let fields = quoted_split(record, ',', 4)?;
            if fields.len() != 4 {
                return Err(bad("fragment annotation needs four fields"));
            }
            let quoted = fields[3]
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .ok_or_else(|| bad("fragment annotation must be quoted"))?;
            let mut annotation = String::new();
            let mut chars = quoted.chars();
            while let Some(c) = chars.next() {
                if c == '\\' {
                    let c = chars
                        .next()
                        .ok_or_else(|| bad("incomplete annotation escape"))?;
                    if !matches!(c, '"' | '\\') {
                        return Err(bad("invalid annotation escape"));
                    }
                    annotation.push(c);
                } else {
                    annotation.push(c);
                }
            }
            Ok(PeakAnnotation {
                mz: finite(fields[0])?,
                intensity: finite(fields[1])?,
                charge: number(fields[2])?,
                annotation,
            })
        })
        .collect()
}
fn analysis_results(meta: &mut MetaInfo) -> Result<Vec<AnalysisResult>> {
    let keys: Vec<_> = meta
        .keys()
        .filter(|k| k.starts_with("_ar_"))
        .cloned()
        .collect();
    let mut results: BTreeMap<usize, BTreeMap<String, MetaValue>> = BTreeMap::new();
    for key in keys {
        let (index, field) = key[4..]
            .split_once('_')
            .ok_or_else(|| bad("invalid analysis result key"))?;
        let numeric_index: usize = number(index)?;
        if index != numeric_index.to_string() {
            return Err(bad(
                "analysis result indices must use canonical unsigned decimal notation",
            ));
        }
        if results
            .entry(numeric_index)
            .or_default()
            .insert(field.into(), meta.remove(&key).unwrap())
            .is_some()
        {
            return Err(bad("duplicate normalized analysis result field"));
        }
    }
    let mut out = Vec::new();
    for (index, mut fields) in results {
        if index != out.len() {
            return Err(bad("analysis result indices must be contiguous from zero"));
        }
        let mut result = AnalysisResult {
            score_type: take_text(&mut fields, "score_type")?
                .ok_or_else(|| bad("analysis result needs score_type"))?,
            main_score: fields
                .remove("score")
                .ok_or_else(|| bad("analysis result needs score"))?
                .as_f64()?,
            higher_is_better: fields
                .remove("higher_is_better")
                .map(|v| v.to_bool())
                .transpose()?
                .unwrap_or(true),
            ..Default::default()
        };
        for (field, value) in fields {
            let name = field
                .strip_prefix("subscore_")
                .ok_or_else(|| unsupported(format!("analysis result field {field}")))?;
            result.sub_scores.insert(name.into(), value.as_f64()?);
        }
        out.push(result);
    }
    Ok(out)
}
fn read_peptide(
    node: &Node,
    identifier: &str,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<PeptideIdentification> {
    node.check(
        &[
            "score_type",
            "higher_score_better",
            "significance_threshold",
            "RT",
            "MZ",
            "spectrum_reference",
        ],
        &["PeptideHit", "UserParam"],
    )?;
    let mut value = PeptideIdentification {
        identifier: identifier.into(),
        score_type: node.get("score_type")?.into(),
        higher_score_better: boolean(node.get("higher_score_better")?)?,
        significance_threshold: finite(node.optional("significance_threshold").unwrap_or("0"))?,
        rt: node.optional("RT").map(finite).transpose()?,
        mz: node.optional("MZ").map(finite).transpose()?,
        metadata: read_meta(node, options)?,
        ..Default::default()
    };
    if let Some(reference) = node.optional("spectrum_reference") {
        insert_meta(&mut value.metadata, "spectrum_reference", reference.into())?;
    }
    for child in node.children.iter().filter(|c| c.name == "PeptideHit") {
        child.check(
            &[
                "score",
                "sequence",
                "charge",
                "protein_refs",
                "aa_before",
                "aa_after",
                "start",
                "end",
            ],
            &["UserParam"],
        )?;
        let mut hit = PeptideHit {
            sequence: AASequence::parse_with_registry(child.get("sequence")?, registry)?,
            score: finite(child.get("score")?)?,
            charge: number(child.get("charge")?)?,
            evidences: read_evidence(child, refs, options)?,
            metadata: read_meta(child, options)?,
            ..Default::default()
        };
        hit.rank = rank(&mut hit.metadata)?;
        if let Some(value) = take_text(&mut hit.metadata, "fragment_annotation")? {
            hit.peak_annotations = annotations(&value, options)?;
        }
        hit.analysis_results = analysis_results(&mut hit.metadata)?;
        value.hits.push(hit);
    }
    value.validate()?;
    Ok(value)
}
fn xml_id(value: &str) -> Result<()> {
    xml_name(value, false)
}
fn xml_name(value: &str, allow_colon: bool) -> Result<()> {
    // IDs use NCName; processing-instruction targets use Name, also allowing ':'.
    let start = |c: char| {
        (allow_colon && c == ':')
            || c == '_'
            || c.is_ascii_alphabetic()
            || ('\u{c0}'..='\u{d6}').contains(&c)
            || ('\u{d8}'..='\u{f6}').contains(&c)
            || ('\u{f8}'..='\u{2ff}').contains(&c)
            || ('\u{370}'..='\u{37d}').contains(&c)
            || ('\u{37f}'..='\u{1fff}').contains(&c)
            || ('\u{200c}'..='\u{200d}').contains(&c)
            || ('\u{2070}'..='\u{218f}').contains(&c)
            || ('\u{2c00}'..='\u{2fef}').contains(&c)
            || ('\u{3001}'..='\u{d7ff}').contains(&c)
            || ('\u{f900}'..='\u{fdcf}').contains(&c)
            || ('\u{fdf0}'..='\u{fffd}').contains(&c)
            || ('\u{10000}'..='\u{effff}').contains(&c)
    };
    let mut chars = value.chars();
    if !chars.next().is_some_and(start)
        || !chars.all(|c| {
            start(c)
                || c.is_ascii_digit()
                || matches!(c, '-' | '.' | '\u{b7}')
                || ('\u{300}'..='\u{36f}').contains(&c)
                || ('\u{203f}'..='\u{2040}').contains(&c)
        })
    {
        return Err(bad("invalid XML ID"));
    }
    Ok(())
}

pub fn read(reader: impl BufRead) -> Result<IdXmlDocument> {
    read_with_options(reader, &ReadOptions::default())
}
/// Reads a bounded XML tree and returns native records only after full validation.
/// Unknown tags/attributes and unsupported encodings are errors, never dropped.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<IdXmlDocument> {
    read_with_registry(reader, options, ModificationsDB::global())
}
/// Read with caller-supplied modification names. Returned sequences retain owned
/// handles, so the registry may be dropped immediately after successful reading.
pub fn read_with_registry(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<IdXmlDocument> {
    let root = parse_xml(reader, options)?;
    root.check(
        &[
            "version",
            "id",
            "xmlns:xsi",
            "xsi:noNamespaceSchemaLocation",
        ],
        &["SearchParameters", "IdentificationRun"],
    )?;
    let version = finite(root.optional("version").unwrap_or("1.0"))?;
    if !(1.0..=1.5).contains(&version) {
        return Err(unsupported("idXML versions outside 1.0 through 1.5"));
    }
    if root.optional("xsi:noNamespaceSchemaLocation").is_some()
        && root.optional("xmlns:xsi").is_none()
    {
        return Err(bad("unbound xsi attribute prefix"));
    }
    if root
        .optional("xmlns:xsi")
        .is_some_and(|v| v != "http://www.w3.org/2001/XMLSchema-instance")
    {
        return Err(bad("invalid xsi namespace"));
    }
    let mut document = IdXmlDocument {
        document_id: root.optional("id").unwrap_or("").into(),
        ..Default::default()
    };
    let mut search = BTreeMap::new();
    let mut search_order = Vec::new();
    let mut ids = BTreeSet::new();
    let mut used = BTreeSet::new();
    let mut run_ids = BTreeSet::new();
    for node in &root.children {
        if node.name == "SearchParameters" {
            let id = node.get("id")?;
            xml_id(id)?;
            if !ids.insert(id.into()) {
                return Err(bad(format!("duplicate XML ID {id}")));
            }
            search.insert(id.to_owned(), read_search(node, options)?);
            search_order.push(id.to_owned());
        } else {
            node.check(
                &[
                    "date",
                    "search_engine",
                    "search_engine_version",
                    "search_parameters_ref",
                ],
                &["ProteinIdentification", "PeptideIdentification"],
            )?;
            let reference = node.get("search_parameters_ref")?;
            let parameters = search
                .get(reference)
                .ok_or_else(|| bad(format!("unknown search parameters {reference}")))?
                .clone();
            used.insert(reference.to_owned());
            let date_time = node.get("date")?;
            date(date_time)?;
            let engine = node.get("search_engine")?;
            let mut run = ProteinIdentification {
                identifier: format!(
                    "{engine}_{date_time}_{}",
                    document.protein_identifications.len()
                ),
                search_engine: engine.into(),
                search_engine_version: node.get("search_engine_version")?.into(),
                search_parameters: parameters,
                date_time: Some(date_time.into()),
                ..Default::default()
            };
            let proteins: Vec<_> = node
                .children
                .iter()
                .filter(|c| c.name == "ProteinIdentification")
                .collect();
            if proteins.len() > 1 {
                return Err(unsupported(
                    "multiple ProteinIdentification blocks in one run",
                ));
            }
            let mut refs = BTreeMap::new();
            if let Some(protein) = proteins.first() {
                run = read_protein(protein, run, &mut refs, &mut ids, options)?;
            }
            if !run_ids.insert(run.identifier.clone()) {
                return Err(bad("duplicate identification run identifier"));
            }
            for peptide in node
                .children
                .iter()
                .filter(|c| c.name == "PeptideIdentification")
            {
                document.peptide_identifications.push(read_peptide(
                    peptide,
                    &run.identifier,
                    &refs,
                    options,
                    registry,
                )?);
            }
            document.protein_identifications.push(run);
        }
    }
    if document.protein_identifications.is_empty() {
        return Err(bad("idXML needs at least one IdentificationRun"));
    }
    for id in search_order {
        if !used.contains(&id) {
            document
                .unreferenced_search_parameters
                .push(search.remove(&id).unwrap());
        }
    }
    Ok(document)
}

fn insert_meta(meta: &mut MetaInfo, key: &str, value: MetaValue) -> Result<()> {
    if let Some(existing) = meta.get(key) {
        if existing != &value {
            return Err(bad(format!("conflicting reserved metadata {key}")));
        }
    } else {
        meta.insert(key.into(), value);
    }
    Ok(())
}
fn write_meta(node: &mut Node, meta: &MetaInfo) -> Result<()> {
    for (key, value) in meta {
        value.validate()?;
        if value.unit().is_some() {
            return Err(unsupported("idXML UserParam cannot represent units"));
        }
        let (kind, text) = match value.data() {
            MetaValueData::Empty => {
                return Err(unsupported(
                    "idXML cannot distinguish Empty metadata from empty string",
                ));
            }
            MetaValueData::String(value) => ("string", value.clone()),
            MetaValueData::Integer(value) => ("int", value.to_string()),
            MetaValueData::Float(value) => ("float", value.to_string()),
            MetaValueData::IntegerList(value) => (
                "intList",
                format!(
                    "[{}]",
                    value
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            ),
            MetaValueData::FloatList(value) => (
                "floatList",
                format!(
                    "[{}]",
                    value
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            ),
            MetaValueData::StringList(value) => {
                if value.iter().any(|v| v.contains("\\|"))
                    || (value.len() == 1 && value[0].is_empty())
                {
                    return Err(unsupported(
                        "ambiguous idXML stringList: literal backslash-pipe or single empty entry",
                    ));
                }
                (
                    "stringList",
                    format!(
                        "[{}]",
                        value
                            .iter()
                            .map(|v| v.replace(',', "\\|"))
                            .collect::<Vec<_>>()
                            .join(",")
                    ),
                )
            }
        };
        let mut child = Node::new("UserParam");
        child.attr("name", key);
        child.attr("type", kind);
        child.attr("value", text);
        node.children.push(child);
    }
    Ok(())
}
fn write_search(value: &SearchParameters, id: &str) -> Result<Node> {
    value.validate()?;
    if value.metadata.contains_key("modification_definitions") {
        return Err(unsupported(
            "tool-defined modification definitions in idXML",
        ));
    }
    if !value.digestion_regex.is_empty() {
        return Err(unsupported(
            "idXML does not encode a custom digestion_regex",
        ));
    }
    let mut node = Node::new("SearchParameters");
    node.attr("id", id);
    node.attr("db", &value.database);
    node.attr("db_version", &value.database_version);
    node.attr("taxonomy", &value.taxonomy);
    node.attr("charges", &value.charges);
    node.attr(
        "mass_type",
        match value.mass_type {
            PeakMassType::Monoisotopic => "monoisotopic",
            PeakMassType::Average => "average",
        },
    );
    node.attr("enzyme", &value.digestion_enzyme);
    node.attr("missed_cleavages", value.missed_cleavages);
    for (tolerance, name, ppm) in [
        (
            value.precursor_tolerance,
            "precursor_peak_tolerance",
            "precursor_peak_tolerance_ppm",
        ),
        (
            value.fragment_tolerance,
            "peak_mass_tolerance",
            "peak_mass_tolerance_ppm",
        ),
    ] {
        let (Tolerance::Absolute(amount) | Tolerance::Ppm(amount)) = tolerance;
        node.attr(name, amount);
        node.attr(ppm, matches!(tolerance, Tolerance::Ppm(_)));
    }
    for (names, tag) in [
        (&value.fixed_modifications, "FixedModification"),
        (&value.variable_modifications, "VariableModification"),
    ] {
        for name in names {
            if name.is_empty() {
                return Err(bad("empty modification name"));
            }
            let mut child = Node::new(tag);
            child.attr("name", name);
            node.children.push(child);
        }
    }
    let mut meta = value.metadata.clone();
    if value.enzyme_specificity != EnzymeTermSpecificity::Unknown {
        insert_meta(
            &mut meta,
            "EnzymeTermSpecificity",
            match value.enzyme_specificity {
                EnzymeTermSpecificity::Full => "full",
                EnzymeTermSpecificity::Semi => "semi",
                EnzymeTermSpecificity::None => "none",
                EnzymeTermSpecificity::Unknown => unreachable!(),
            }
            .into(),
        )?;
    } else if meta.contains_key("EnzymeTermSpecificity") {
        return Err(bad(
            "EnzymeTermSpecificity metadata must be represented by its typed field",
        ));
    }
    write_meta(&mut node, &meta)?;
    Ok(node)
}
fn write_groups(
    meta: &mut MetaInfo,
    groups: &[ProteinGroup],
    prefix: &str,
    refs: &BTreeMap<String, String>,
) -> Result<()> {
    if meta
        .keys()
        .any(|key| key.starts_with(&format!("{prefix}_")))
    {
        return Err(bad("reserved protein group metadata collision"));
    }
    for (index, group) in groups.iter().enumerate() {
        if !group.float_data_arrays.is_empty()
            || !group.integer_data_arrays.is_empty()
            || !group.string_data_arrays.is_empty()
        {
            return Err(unsupported(
                "idXML protein groups do not encode sample arrays",
            ));
        }
        if group.accessions.is_empty() {
            return Err(bad("idXML protein group needs accessions"));
        }
        let mut values = vec![group.probability.to_string()];
        for accession in &group.accessions {
            values.push(refs.get(accession).cloned().ok_or_else(|| {
                bad(format!(
                    "protein group accession absent from run: {accession}"
                ))
            })?);
        }
        meta.insert(format!("{prefix}_{index}"), values.join(",").into());
    }
    Ok(())
}
fn write_protein(
    value: &ProteinIdentification,
    count: &mut usize,
) -> Result<(Node, BTreeMap<String, String>)> {
    let mut node = Node::new("ProteinIdentification");
    node.attr("score_type", &value.score_type);
    node.attr("higher_score_better", value.higher_score_better);
    node.attr("significance_threshold", value.significance_threshold);
    let mut refs = BTreeMap::new();
    for hit in &value.hits {
        if !hit.modifications.is_empty() {
            return Err(unsupported(
                "idXML does not encode ProteinHit modification positions",
            ));
        }
        let id = format!("PH_{count}");
        *count = count
            .checked_add(1)
            .ok_or_else(|| bad("protein reference count overflow"))?;
        if refs.insert(hit.accession.clone(), id.clone()).is_some() {
            return Err(bad("duplicate protein accession within run is ambiguous"));
        }
        let mut child = Node::new("ProteinHit");
        child.attr("id", id);
        child.attr("accession", &hit.accession);
        child.attr("sequence", &hit.sequence);
        child.attr("score", hit.score);
        if let Some(coverage) = hit.coverage {
            child.attr("coverage", coverage);
        }
        let mut meta = hit.metadata.clone();
        if meta.contains_key(RANK) {
            return Err(bad("reserved rank metadata collision"));
        }
        if hit.rank != 0 {
            meta.insert(RANK.into(), hit.rank.into());
        }
        write_meta(&mut child, &meta)?;
        node.children.push(child);
    }
    let mut meta = value.metadata.clone();
    if meta.contains_key(IDENTIFIER) {
        return Err(bad("reserved run identifier metadata collision"));
    }
    meta.insert(IDENTIFIER.into(), value.identifier.clone().into());
    for (key, paths) in [
        ("spectra_data", &value.primary_ms_run_paths),
        ("spectra_data_raw", &value.raw_ms_run_paths),
    ] {
        if meta.contains_key(key) {
            return Err(bad(format!(
                "{key} metadata must use the typed run-path field"
            )));
        }
        if !paths.is_empty() {
            meta.insert(key.into(), paths.clone().into());
        }
    }
    write_groups(&mut meta, &value.protein_groups, "protein_group", &refs)?;
    write_groups(
        &mut meta,
        &value.indistinguishable_groups,
        "indistinguishable_proteins",
        &refs,
    )?;
    write_meta(&mut node, &meta)?;
    Ok((node, refs))
}
fn write_evidence(
    node: &mut Node,
    values: &[PeptideEvidence],
    refs: &BTreeMap<String, String>,
) -> Result<()> {
    if values.is_empty() {
        return Ok(());
    }
    // XML IDREFS has no empty-entry placeholder. Only a suffix of empty
    // accessions is representable without shifting positional correspondence.
    let prefix = values
        .iter()
        .take_while(|e| !e.protein_accession.is_empty())
        .count();
    if values[prefix..]
        .iter()
        .any(|e| !e.protein_accession.is_empty())
    {
        return Err(unsupported(
            "evidence has an empty accession before a nonempty one",
        ));
    }
    if prefix > 0 {
        let ids = values[..prefix]
            .iter()
            .map(|e| {
                refs.get(&e.protein_accession).cloned().ok_or_else(|| {
                    bad(format!(
                        "evidence accession absent from run: {}",
                        e.protein_accession
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        node.attr("protein_refs", ids.join(" "));
    }
    // Emit both marker lists even if all unknown, retaining empty-accession
    // evidence entries that C++ would otherwise silently drop.
    node.attr(
        "aa_before",
        values
            .iter()
            .map(|e| e.aa_before.code().to_string())
            .collect::<Vec<_>>()
            .join(" "),
    );
    node.attr(
        "aa_after",
        values
            .iter()
            .map(|e| e.aa_after.code().to_string())
            .collect::<Vec<_>>()
            .join(" "),
    );
    for (name, positions) in [
        ("start", values.iter().map(|e| e.start).collect::<Vec<_>>()),
        ("end", values.iter().map(|e| e.end).collect::<Vec<_>>()),
    ] {
        if positions.iter().any(Option::is_some) {
            let values = positions
                .into_iter()
                .map(|p| {
                    p.map(|v| i32::try_from(v).map_err(|_| bad("evidence position exceeds xs:int")))
                        .transpose()
                        .map(|v| v.unwrap_or(-1).to_string())
                })
                .collect::<Result<Vec<_>>>()?;
            node.attr(name, values.join(" "));
        }
    }
    Ok(())
}
fn write_peptide(
    value: &PeptideIdentification,
    refs: &BTreeMap<String, String>,
    registry: &ModificationsDB,
) -> Result<Node> {
    value.validate()?;
    let mut node = Node::new("PeptideIdentification");
    node.attr("score_type", &value.score_type);
    node.attr("higher_score_better", value.higher_score_better);
    node.attr("significance_threshold", value.significance_threshold);
    if let Some(rt) = value.rt {
        node.attr("RT", rt);
    }
    if let Some(mz) = value.mz {
        node.attr("MZ", mz);
    }
    let mut meta = value.metadata.clone();
    if let Some(reference) = meta.remove("spectrum_reference") {
        if reference.unit().is_some() {
            return Err(unsupported("spectrum_reference unit"));
        }
        node.attr("spectrum_reference", reference.as_str()?);
    }
    for hit in &value.hits {
        let mut child = Node::new("PeptideHit");
        let sequence_text = hit.sequence.to_string();
        // Some source generator paths attach terminal modifications to residue
        // slots. idXML's sequence syntax cannot preserve every such typed state.
        if !matches!(AASequence::parse_with_registry(&sequence_text, registry), Ok(ref parsed) if parsed == &hit.sequence)
        {
            return Err(unsupported(
                "peptide modification placement in idXML sequence syntax",
            ));
        }
        child.attr("sequence", sequence_text);
        child.attr("score", hit.score);
        child.attr("charge", hit.charge);
        write_evidence(&mut child, &hit.evidences, refs)?;
        let mut meta = hit.metadata.clone();
        if meta.contains_key(RANK)
            || meta.contains_key("fragment_annotation")
            || meta.keys().any(|k| k.starts_with("_ar_"))
        {
            return Err(bad("reserved peptide hit metadata collision"));
        }
        if hit.rank != 0 {
            meta.insert(RANK.into(), hit.rank.into());
        }
        if !hit.peak_annotations.is_empty() {
            let values = hit
                .peak_annotations
                .iter()
                .map(|a| {
                    format!(
                        "{},{},{},\"{}\"",
                        a.mz,
                        a.intensity,
                        a.charge,
                        a.annotation.replace('\\', "\\\\").replace('"', "\\\"")
                    )
                })
                .collect::<Vec<_>>();
            meta.insert("fragment_annotation".into(), values.join("|").into());
        }
        for (i, result) in hit.analysis_results.iter().enumerate() {
            meta.insert(
                format!("_ar_{i}_score_type"),
                result.score_type.clone().into(),
            );
            meta.insert(
                format!("_ar_{i}_score"),
                MetaValue::try_from(result.main_score)?,
            );
            meta.insert(
                format!("_ar_{i}_higher_is_better"),
                result.higher_is_better.to_string().into(),
            );
            for (key, score) in &result.sub_scores {
                meta.insert(
                    format!("_ar_{i}_subscore_{key}"),
                    MetaValue::try_from(*score)?,
                );
            }
        }
        write_meta(&mut child, &meta)?;
        node.children.push(child);
    }
    write_meta(&mut node, &meta)?;
    Ok(node)
}
struct Output {
    bytes: Vec<u8>,
    records: usize,
    options: WriteOptions,
}
impl Output {
    fn append(&mut self, text: &str) -> Result<()> {
        if self
            .bytes
            .len()
            .checked_add(text.len())
            .is_none_or(|n| n > self.options.max_xml_bytes)
        {
            return Err(bad("idXML output byte limit exceeded"));
        }
        self.bytes.extend_from_slice(text.as_bytes());
        Ok(())
    }
    fn node(&mut self, node: &Node, depth: usize) -> Result<()> {
        self.records += 1;
        if self.records > self.options.max_records {
            return Err(bad("idXML output element limit exceeded"));
        }
        self.append(&"  ".repeat(depth))?;
        self.append(&format!("<{}", node.name))?;
        for (key, value) in &node.attrs {
            xml_text(value)?;
            if value.len() > self.options.max_xml_bytes {
                return Err(bad("XML attribute exceeds output limit"));
            }
            self.append(&format!(" {key}=\"{}\"", escape(value)))?;
        }
        if node.children.is_empty() {
            self.append("/>\n")?;
        } else {
            self.append(">\n")?;
            for child in &node.children {
                self.node(child, depth + 1)?;
            }
            self.append(&format!("{}</{}>\n", "  ".repeat(depth), node.name))?;
        }
        Ok(())
    }
}
pub fn write(writer: impl Write, document: &IdXmlDocument) -> Result<()> {
    write_with_options(writer, document, &WriteOptions::default())
}
/// Builds and validates the complete output before touching the caller's writer.
/// An I/O error during the final write can still leave a partial external file.
pub fn write_with_options(
    writer: impl Write,
    document: &IdXmlDocument,
    options: &WriteOptions,
) -> Result<()> {
    write_with_registry(writer, document, options, ModificationsDB::global())
}
/// Validate sequence reconstruction against an explicitly supplied registry before
/// writing. A reader needs equivalent registry chemistry to recover these records.
pub fn write_with_registry(
    mut writer: impl Write,
    document: &IdXmlDocument,
    options: &WriteOptions,
    registry: &ModificationsDB,
) -> Result<()> {
    if document.protein_identifications.is_empty() {
        return Err(bad("idXML writer needs at least one protein run"));
    }
    if options.max_records == 0 || options.max_xml_bytes == 0 {
        return Err(bad("idXML write limits must be positive"));
    }
    let mut run_ids = BTreeSet::new();
    for run in &document.protein_identifications {
        run.validate()?;
        if !run_ids.insert(run.identifier.as_str()) {
            return Err(bad("duplicate identification run identifier"));
        }
        date(
            run.date_time
                .as_deref()
                .ok_or_else(|| bad("IdentificationRun requires date_time"))?,
        )?;
    }
    let mut peptides: BTreeMap<&str, Vec<&PeptideIdentification>> = BTreeMap::new();
    for peptide in &document.peptide_identifications {
        if !run_ids.contains(peptide.identifier.as_str()) {
            return Err(bad("peptide identifier does not match any protein run"));
        }
        peptides
            .entry(&peptide.identifier)
            .or_default()
            .push(peptide);
    }
    let mut root = Node::new("IdXML");
    root.attr("version", "1.5");
    root.attr("id", &document.document_id);
    root.attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.attr(
        "xsi:noNamespaceSchemaLocation",
        "https://www.openms.de/xml-schema/IdXML_1_5.xsd",
    );
    // Keep one parameter block per run: unlike C++ equality-based deduplication,
    // this cannot collapse search parameter metadata from distinct runs.
    for (i, search) in document
        .protein_identifications
        .iter()
        .map(|p| &p.search_parameters)
        .chain(&document.unreferenced_search_parameters)
        .enumerate()
    {
        root.children
            .push(write_search(search, &format!("SP_{i}"))?);
    }
    let mut protein_count = 0;
    for (i, run) in document.protein_identifications.iter().enumerate() {
        let mut node = Node::new("IdentificationRun");
        node.attr("date", run.date_time.as_deref().unwrap());
        node.attr("search_engine", &run.search_engine);
        node.attr("search_engine_version", &run.search_engine_version);
        node.attr("search_parameters_ref", format!("SP_{i}"));
        let (protein, refs) = write_protein(run, &mut protein_count)?;
        node.children.push(protein);
        for peptide in peptides.get(run.identifier.as_str()).into_iter().flatten() {
            node.children.push(write_peptide(peptide, &refs, registry)?);
        }
        root.children.push(node);
    }
    let mut output = Output {
        bytes: Vec::new(),
        records: 0,
        options: *options,
    };
    output.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
    output.node(&root, 0)?;
    writer.write_all(&output.bytes)?;
    writer.flush()?;
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared bounded identification XML leaf codecs; dialects own run/group layout.

use crate::chemistry::{AASequence, ModificationsDB};
use crate::comparison::Tolerance;
use crate::identification::{
    AnalysisResult, EnzymeTermSpecificity, FlankingResidue, PeakAnnotation, PeakMassType,
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinIdentification, SearchParameters,
};
use crate::metadata::{CompletionTime, MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use quick_xml::{
    NsReader, Reader,
    events::{BytesStart, Event},
    name::ResolveResult,
};
use std::collections::BTreeMap;
use std::io::{BufRead, Read};

pub(crate) const IDENTIFIER: &str = "openms-rust:run_identifier";
pub(crate) const RANK: &str = "openms-rust:rank";
/// Ceilings the shared identification-XML parser applies to one document.
///
/// A dialect whose ceilings grow with the document it is reading computes these
/// from the decoded size and passes the result; see
/// `src/format/featurexml_scaling.rs`.
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
/// Ceilings the shared identification-XML renderer applies to one document.
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
pub(crate) fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
pub(crate) fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
pub(crate) fn number<T: std::str::FromStr>(value: &str) -> Result<T> {
    value
        .trim()
        .parse()
        .map_err(|_| bad(format!("invalid number {value:?}")))
}
pub(crate) fn finite(value: &str) -> Result<f64> {
    let value = number::<f64>(value)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("idXML numeric values must be finite"))
    }
}
pub(crate) fn boolean(value: &str) -> Result<bool> {
    match value {
        "true" | "1" => Ok(true),
        "false" | "0" => Ok(false),
        _ => Err(bad("invalid XML boolean")),
    }
}
pub(crate) fn xml_text(value: &str) -> Result<()> {
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
pub(crate) fn date(value: &str) -> Result<()> {
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

#[derive(Clone, Debug, Default)]
pub(crate) struct Node {
    pub(crate) name: String,
    pub(crate) text: String,
    pub(crate) attrs: BTreeMap<String, String>,
    pub(crate) children: Vec<Node>,
}
impl Node {
    pub(crate) fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
    pub(crate) fn attr(&mut self, name: &str, value: impl ToString) {
        self.attrs.insert(name.into(), value.to_string());
    }
    pub(crate) fn get(&self, key: &str) -> Result<&str> {
        self.attrs
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| bad(format!("{} requires {key}", self.name)))
    }
    pub(crate) fn optional(&self, key: &str) -> Option<&str> {
        self.attrs.get(key).map(String::as_str)
    }
    pub(crate) fn check(&self, attrs: &[&str], children: &[&str]) -> Result<()> {
        if !self.text.trim_matches([' ', '\t', '\r', '\n']).is_empty() {
            return Err(bad(format!("{} cannot contain text", self.name)));
        }
        self.check_with_text(attrs, children)
    }
    pub(crate) fn check_with_text(&self, attrs: &[&str], children: &[&str]) -> Result<()> {
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

/// Read at most `limit` bytes, growing the buffer in bounded steps whose
/// failure is a checked error rather than an allocation abort.
///
/// `Read::read_to_end` doubles its buffer with the infallible allocator, which
/// aborts the process when a multi-gigabyte document does not fit. Reserving
/// each step with `try_reserve` keeps a document larger than the host can hold
/// a refusal, as every other ceiling in this reader is.
fn read_all_bounded(mut input: impl Read, limit: usize) -> Result<Vec<u8>> {
    /// Bytes read per step; also the initial reservation.
    const STEP: usize = 1 << 20;
    let mut bytes: Vec<u8> = Vec::new();
    loop {
        if bytes.len() == bytes.capacity() {
            let step = STEP.min(limit.saturating_add(1).saturating_sub(bytes.len()));
            if step == 0 {
                break;
            }
            bytes
                .try_reserve(step)
                .map_err(|_| bad("identification XML allocation failed"))?;
        }
        let read = {
            let spare = bytes.capacity() - bytes.len();
            let start = bytes.len();
            bytes.resize(start + spare, 0);
            let read = input.read(&mut bytes[start..])?;
            bytes.truncate(start + read);
            read
        };
        if read == 0 {
            break;
        }
        if bytes.len() > limit {
            return Err(bad("identification XML byte limit exceeded"));
        }
    }
    if bytes.len() > limit {
        return Err(bad("identification XML byte limit exceeded"));
    }
    Ok(bytes)
}

fn document(input: impl Read, limit: usize) -> Result<String> {
    let mut bytes = read_all_bounded(input, limit)?;
    let utf16 = if bytes.starts_with(&[0xff, 0xfe]) {
        Some((true, 2))
    } else if bytes.starts_with(&[0xfe, 0xff]) {
        Some((false, 2))
    } else if bytes.starts_with(&[b'<', 0, b'?', 0]) {
        Some((true, 0))
    } else if bytes.starts_with(&[0, b'<', 0, b'?']) {
        Some((false, 0))
    } else {
        None
    };
    let text = if let Some((little, offset)) = utf16 {
        if (bytes.len() - offset) % 2 != 0 {
            return Err(bad("odd UTF-16 byte count"));
        }
        let units = bytes[offset..].chunks_exact(2).map(|b| {
            if little {
                u16::from_le_bytes([b[0], b[1]])
            } else {
                u16::from_be_bytes([b[0], b[1]])
            }
        });
        let mut decoded = String::new();
        for value in char::decode_utf16(units) {
            let c = value.map_err(|_| bad("invalid UTF-16 identification XML"))?;
            if decoded.len().saturating_add(c.len_utf8()) > limit {
                return Err(bad("decoded identification XML byte limit exceeded"));
            }
            decoded.push(c);
        }
        decoded
    } else {
        let bom = bytes.starts_with(&[0xef, 0xbb, 0xbf]);
        if bom {
            bytes.drain(..3);
        }
        let encoding = encoding(&bytes)?;
        if bom
            && encoding
                .as_deref()
                .is_some_and(|s| !matches!(s, "utf-8" | "utf8"))
        {
            return Err(bad("XML declaration conflicts with UTF-8 byte order mark"));
        }
        // Every branch that can take the buffer as it stands does, because a
        // copy of a multi-gigabyte document doubles the reader's peak memory.
        match encoding.as_deref().unwrap_or("utf-8") {
            "utf-8" | "utf8" => {
                String::from_utf8(bytes).map_err(|_| bad("invalid UTF-8 identification XML"))?
            }
            "us-ascii" | "ascii" | "iso-8859-1" | "iso8859-1" | "latin1" if bytes.is_ascii() => {
                String::from_utf8(bytes).map_err(|_| bad("invalid ASCII identification XML"))?
            }
            "iso-8859-1" | "iso8859-1" | "latin1" => {
                // Two UTF-8 bytes per Latin-1 byte is the exact worst case, so
                // one fallible reservation covers the whole decode.
                let mut decoded = String::new();
                decoded
                    .try_reserve_exact(bytes.len().saturating_mul(2).min(limit))
                    .map_err(|_| bad("identification XML allocation failed"))?;
                for &byte in &bytes {
                    let c = char::from(byte);
                    if decoded.len().saturating_add(c.len_utf8()) > limit {
                        return Err(bad("decoded identification XML byte limit exceeded"));
                    }
                    decoded.push(c);
                }
                decoded
            }
            value => return Err(unsupported(format!("identification XML encoding {value}"))),
        }
    };
    let declared = encoding(text.as_bytes())?;
    if let Some((little, _)) = utf16 {
        if let Some(name) = declared.as_deref() {
            if !matches!(name, "utf-16" | "utf16")
                && name != if little { "utf-16le" } else { "utf-16be" }
            {
                return Err(bad("XML declaration conflicts with UTF-16 byte encoding"));
            }
        }
    }

    // XML 1.0 line ending normalization happens before attribute normalization.
    // `str::replace` always allocates, so a document that already has no
    // carriage return keeps its single buffer instead of being copied twice.
    if text.as_bytes().contains(&b'\r') {
        return Ok(text.replace("\r\n", "\n").replace('\r', "\n"));
    }
    Ok(text)
}

fn encoding(bytes: &[u8]) -> Result<Option<String>> {
    let mut reader = Reader::from_reader(bytes);
    match reader.read_event().map_err(|e| bad(e.to_string()))? {
        Event::Decl(declaration) => {
            let raw = std::str::from_utf8(declaration.as_ref())
                .map_err(|_| bad("invalid XML declaration"))?;
            let element = BytesStart::from_content(raw, 3);
            attribute_spacing(&element)?;
            let mut phase = 0;
            for attr in element.attributes() {
                let attr = attr.map_err(|e| bad(e.to_string()))?;
                match attr.key.as_ref() {
                    b"version" if phase == 0 && attr.value.as_ref() == b"1.0" => phase = 1,
                    b"encoding"
                        if phase == 1
                            && attr.value.first().is_some_and(u8::is_ascii_alphabetic)
                            && attr.value.iter().all(|b| {
                                b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-')
                            }) =>
                    {
                        phase = 2
                    }
                    b"standalone"
                        if matches!(phase, 1 | 2)
                            && matches!(attr.value.as_ref(), b"yes" | b"no") =>
                    {
                        phase = 3
                    }
                    _ => return Err(bad("invalid XML declaration attributes or order")),
                }
            }
            if phase == 0 {
                return Err(bad("XML declaration needs a version"));
            }
            if declaration
                .version()
                .map_err(|e| bad(e.to_string()))?
                .as_ref()
                != b"1.0"
            {
                return Err(unsupported("only XML version 1.0 is supported"));
            }
            declaration
                .encoding()
                .transpose()
                .map_err(|e| bad(e.to_string()))?
                .map(|value| {
                    std::str::from_utf8(&value)
                        .map(str::to_ascii_lowercase)
                        .map_err(|_| bad("invalid XML encoding name"))
                })
                .transpose()
        }
        _ => Ok(None),
    }
}

fn attribute_spacing(element: &BytesStart<'_>) -> Result<()> {
    let tail = &element.as_ref()[element.name().as_ref().len()..];
    let space = |b| matches!(b, b' ' | b'\t' | b'\r' | b'\n');
    if tail.first().is_some_and(|b| !space(*b)) {
        return Err(bad("XML attributes require whitespace separators"));
    }
    let mut quote = None;
    let mut closed = false;
    for &byte in tail {
        if let Some(delimiter) = quote {
            if byte == delimiter {
                quote = None;
                closed = true;
            }
        } else if closed {
            if !space(byte) {
                return Err(bad("XML attributes require whitespace separators"));
            }
            closed = false;
        } else if matches!(byte, b'\'' | b'"') {
            quote = Some(byte);
        }
    }
    Ok(())
}
/// Cumulative decoded tree and parser-work limits, separate from input bytes.
#[cfg(feature = "idxml")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct XmlLimits {
    pub(crate) max_depth: usize,
    pub(crate) max_payload_bytes: usize,
    pub(crate) max_work: usize,
}
#[cfg(feature = "idxml")]
impl Default for XmlLimits {
    fn default() -> Self {
        Self {
            max_depth: 260,
            max_payload_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
        }
    }
}
struct XmlMeter<'a> {
    bytes: &'a mut usize,
    work: &'a mut usize,
}
impl XmlMeter<'_> {
    fn add(&mut self, bytes: usize, work: usize) -> Result<()> {
        *self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| bad("XML payload limit exceeded"))?;
        *self.work = self
            .work
            .checked_sub(work)
            .ok_or_else(|| bad("XML work limit exceeded"))?;
        Ok(())
    }
}
/// One root child whose own children are handed over as they close, instead of
/// being retained in the returned tree.
///
/// A dialect whose payload is a long flat list — featureXML's `featureList`,
/// whose `feature` children are 99% of every real document — converts each
/// child and drops it, so the tree in memory stays the size of one child
/// rather than of the whole document. The payload charged for the subtree is
/// refunded once `take` returns, because that storage is no longer held; the
/// work charged for it is not, because the effort was really spent.
pub(crate) struct Detach<'a> {
    /// Name of the root child whose children are detached, e.g. `featureList`.
    pub(crate) container: &'a str,
    /// Name of the child element handed to `take`, e.g. `feature`.
    pub(crate) element: &'a str,
    /// Receives the root as parsed so far, the completed element, and the
    /// remaining work and payload budgets, so that the conversion is charged
    /// against the same ceilings. The root carries every child that closed
    /// before the container opened, which is all of a schema-valid featureMap's
    /// metadata and identification data.
    pub(crate) take: &'a mut dyn FnMut(&Node, Node, &mut usize, &mut usize) -> Result<()>,
}

/// Decode a whole document, or only its prefix through the opening `stop_tag`,
/// to text bounded by `limit` bytes.
///
/// Splitting this from [`parse_text_with_budget`] lets a caller learn the
/// decoded size before it decides the cumulative ceilings that size earns; see
/// `src/format/featurexml_scaling.rs`.
pub(crate) fn decode_document(
    input: impl BufRead,
    limit: usize,
    stop_tag: Option<&str>,
) -> Result<String> {
    match stop_tag {
        Some(stop) => document(read_xml_prefix(input, limit, stop)?.as_slice(), limit),
        None => document(input, limit),
    }
}

pub(crate) fn parse_xml_with_budget(
    input: impl BufRead,
    options: &ReadOptions,
    max_depth: usize,
    stop_tag: Option<&str>,
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
) -> Result<Node> {
    // Checked before the input is touched, as it was before decoding and
    // parsing became separable: invalid ceilings refuse without reading.
    if options.max_records == 0 || options.max_list_items == 0 || max_depth == 0 || max_depth > 512
    {
        return Err(bad("invalid identification XML limits"));
    }
    let max_bytes =
        usize::try_from(options.max_xml_bytes).map_err(|_| bad("XML byte limit overflows"))?;
    let decode_limit = max_bytes.min(*remaining_bytes / 8).min(*remaining_work / 4);
    let text = decode_document(input, decode_limit, stop_tag)?;
    parse_text_with_budget(
        &text,
        options,
        max_depth,
        stop_tag,
        remaining_work,
        remaining_bytes,
        None,
    )
}

/// Parse already-decoded document text into a [`Node`] tree under the shared
/// cumulative budgets, optionally detaching one container's children.
pub(crate) fn parse_text_with_budget(
    text: &str,
    options: &ReadOptions,
    max_depth: usize,
    stop_tag: Option<&str>,
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
    mut detach: Option<Detach<'_>>,
) -> Result<Node> {
    if options.max_records == 0 || options.max_list_items == 0 || max_depth == 0 || max_depth > 512
    {
        return Err(bad("invalid identification XML limits"));
    }
    let mut meter = XmlMeter {
        bytes: remaining_bytes,
        work: remaining_work,
    };
    // Account input/decoded storage and normalization scratch before tree copies.
    meter.add(text.len().saturating_mul(8), text.len().saturating_mul(4))?;
    let mut reader = NsReader::from_reader(text.as_bytes());
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    let mut stack: Vec<Node> = Vec::new();
    let mut root = None;
    let mut count = 0usize;
    let mut declared = false;
    let mut at_start = true;
    // Payload budget as it stood when the open detached element began, so the
    // subtree's storage charge can be returned once the element is handed over.
    let mut detached_from: Option<usize> = None;
    loop {
        let (ns, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| bad(e.to_string()))?;
        let unbound = matches!(ns, ResolveResult::Unbound);
        meter.add(0, 1)?;
        match event {
            Event::Start(start) => {
                let payload_before_element = *meter.bytes;
                if !unbound {
                    return Err(unsupported(
                        "identification XML elements must have no namespace",
                    ));
                }
                if stack.is_empty() && root.is_some() {
                    return Err(bad("multiple XML roots"));
                }
                count += 1;
                if count > options.max_records || stack.len() >= max_depth {
                    return Err(bad("XML element/depth limit exceeded"));
                }
                attribute_spacing(&start)?;
                let mut name = std::str::from_utf8(start.name().as_ref())
                    .map_err(|e| bad(e.to_string()))?
                    .to_owned();
                if name == "userParam" {
                    name = "UserParam".into();
                }
                xml_name(&name, false)?;
                meter.add(
                    size_of::<Node>()
                        .saturating_mul(5)
                        .saturating_add(name.len()),
                    start.len(),
                )?;
                let mut node = Node::new(&name);
                for attr in start.attributes() {
                    let attr = attr.map_err(|e| bad(e.to_string()))?;
                    let key =
                        std::str::from_utf8(attr.key.as_ref()).map_err(|e| bad(e.to_string()))?;
                    xml_name(key, true)?;
                    if attr.value.contains(&b'<') {
                        return Err(bad("unescaped < in XML attribute"));
                    }
                    meter.add(
                        key.len()
                            .saturating_add(attr.value.len().saturating_mul(6))
                            .saturating_add(256),
                        attr.value.len().saturating_mul(4).saturating_add(key.len()),
                    )?;
                    let raw = std::str::from_utf8(&attr.value).map_err(|e| bad(e.to_string()))?;
                    let normalized = raw.replace(['\r', '\n', '\t'], " ");
                    let value = quick_xml::escape::unescape(&normalized)
                        .map_err(|e| bad(e.to_string()))?
                        .into_owned();
                    xml_text(&value)?;
                    if node.attrs.insert(key.to_owned(), value).is_some() {
                        return Err(bad("duplicate XML attribute"));
                    }
                }
                if stack.len() == 1 && stop_tag == Some(node.name.as_str()) {
                    let mut root = stack.pop().ok_or_else(|| bad("unexpected XML end"))?;
                    root.children.push(node);
                    return Ok(root);
                }
                if let Some(detach) = detach.as_ref() {
                    if stack.len() == 2
                        && stack[1].name == detach.container
                        && node.name == detach.element
                    {
                        detached_from = Some(payload_before_element);
                    }
                }
                stack.push(node);
            }
            Event::End(_) => {
                let node = stack.pop().ok_or_else(|| bad("unexpected XML end"))?;
                let handed_over = match (detach.as_mut(), detached_from) {
                    (Some(detach), Some(payload))
                        if stack.len() == 2 && node.name == detach.element =>
                    {
                        detached_from = None;
                        let before = *meter.bytes;
                        (detach.take)(&stack[0], node, meter.work, meter.bytes)?;
                        // Return the subtree's storage charge and keep only what
                        // `take` charged for what it retained.
                        let retained = before.saturating_sub(*meter.bytes);
                        *meter.bytes = payload.saturating_sub(retained);
                        None
                    }
                    _ => Some(node),
                };
                if let Some(node) = handed_over {
                    if let Some(parent) = stack.last_mut() {
                        parent.children.push(node);
                    } else {
                        root = Some(node);
                    }
                }
            }
            Event::Text(text) => {
                let decoded = text.xml_content().map_err(|e| bad(e.to_string()))?;
                meter.add(
                    decoded.len().saturating_mul(4),
                    decoded.len().saturating_mul(2),
                )?;
                let decoded =
                    quick_xml::escape::unescape(&decoded).map_err(|e| bad(e.to_string()))?;
                xml_text(&decoded)?;
                if let Some(node) = stack.last_mut() {
                    node.text.push_str(&decoded);
                } else if !decoded.trim_matches([' ', '\t', '\n', '\r']).is_empty() {
                    return Err(bad("text outside XML root"));
                }
            }
            Event::GeneralRef(reference) => {
                // quick-xml resolves the reference, as in the mzML readers.
                // `resolve_char_ref` applies the XML 1.0 CharRef grammar, so
                // `&#X2E;`, a signed number and `&#0;` are refused. Without a
                // DTD, which this reader refuses, only the five predefined
                // entities exist. `resolve_predefined_entity` is avoided because
                // quick-xml's `escape-html` feature would switch it to HTML5.
                let name = reference.decode().map_err(|e| bad(e.to_string()))?;
                let mut utf8 = [0u8; 4];
                let text: &str = match reference
                    .resolve_char_ref()
                    .map_err(|e| bad(e.to_string()))?
                {
                    Some(c) => c.encode_utf8(&mut utf8),
                    None => quick_xml::escape::resolve_xml_entity(&name)
                        .ok_or_else(|| unsupported("external XML entity"))?,
                };
                xml_text(text)?;
                meter.add(8, name.len())?;
                stack
                    .last_mut()
                    .ok_or_else(|| bad("entity outside XML root"))?
                    .text
                    .push_str(text);
            }
            Event::CData(text) => {
                let text = text.decode().map_err(|e| bad(e.to_string()))?;
                meter.add(text.len().saturating_mul(2), text.len())?;
                xml_text(&text)?;
                stack
                    .last_mut()
                    .ok_or_else(|| bad("CDATA outside XML root"))?
                    .text
                    .push_str(&text);
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
            Event::Decl(_) => {
                if !at_start || declared || !stack.is_empty() || root.is_some() {
                    return Err(bad("misplaced XML declaration"));
                }
                declared = true;
            }
            Event::Eof => break,
            _ => {
                return Err(unsupported(
                    "DTD and external entities in identification XML",
                ));
            }
        }
        at_start = false;
        buffer.clear();
    }
    if !stack.is_empty() {
        return Err(bad("truncated identification XML"));
    }
    root.ok_or_else(|| bad("missing XML root"))
}

pub(crate) fn list<'a>(value: &'a str, options: &ReadOptions) -> Result<Vec<&'a str>> {
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
pub(crate) fn read_meta(node: &Node, options: &ReadOptions) -> Result<MetaInfo> {
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
pub(crate) fn take_text(meta: &mut MetaInfo, key: &str) -> Result<Option<String>> {
    meta.remove(key)
        .map(|v| v.as_str().map(str::to_owned))
        .transpose()
}
pub(crate) fn rank(meta: &mut MetaInfo) -> Result<u32> {
    meta.remove(RANK)
        .map(|v| u32::try_from(v.as_i64()?).map_err(|_| bad("invalid rank")))
        .transpose()
        .map(|v| v.unwrap_or(0))
}
pub(crate) fn specificity(value: &str) -> Result<EnzymeTermSpecificity> {
    match value {
        "unknown" => Ok(EnzymeTermSpecificity::Unknown),
        "full" => Ok(EnzymeTermSpecificity::Full),
        "semi" => Ok(EnzymeTermSpecificity::Semi),
        "none" => Ok(EnzymeTermSpecificity::None),
        _ => Err(bad("invalid EnzymeTermSpecificity")),
    }
}
pub(crate) fn read_search(node: &Node, options: &ReadOptions) -> Result<SearchParameters> {
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

pub(crate) fn words<'a>(text: Option<&'a str>, options: &ReadOptions) -> Result<Vec<&'a str>> {
    let mut result = Vec::new();
    for word in text.unwrap_or("").split_ascii_whitespace() {
        if result.len() == options.max_list_items {
            return Err(bad("evidence list exceeds limit"));
        }
        result.push(word);
    }
    Ok(result)
}
pub(crate) fn position(value: &str) -> Result<Option<usize>> {
    let value: i32 = number(value)?;
    if value == -1 {
        Ok(None)
    } else {
        usize::try_from(value)
            .map(Some)
            .map_err(|_| bad("invalid negative evidence position"))
    }
}
pub(crate) fn flank(value: &str) -> Result<FlankingResidue> {
    let mut chars = value.chars();
    let first = chars.next().ok_or_else(|| bad("empty flank"))?;
    if chars.next().is_some() {
        return Err(bad("flanking residue requires one character"));
    }
    FlankingResidue::from_code(first)
}
pub(crate) fn read_evidence(
    node: &Node,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
    work: &mut usize,
    bytes: &mut usize,
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
    let mut meter = IdMeter { work, bytes };
    meter.slots(len)?;
    let mut result = vec![PeptideEvidence::default(); len];
    // Short optional columns populate their prefix, matching the source loader.
    for (i, value) in refs_column.into_iter().enumerate() {
        let accession = refs
            .get(value)
            .ok_or_else(|| bad(format!("unknown protein reference {value}")))?;
        meter.text(accession)?;
        result[i].protein_accession = accession.clone();
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
pub(crate) fn quoted_split(value: &str, separator: char, max_parts: usize) -> Result<Vec<&str>> {
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
pub(crate) fn annotations(value: &str, options: &ReadOptions) -> Result<Vec<PeakAnnotation>> {
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
pub(crate) fn analysis_results(meta: &mut MetaInfo) -> Result<Vec<AnalysisResult>> {
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
#[allow(clippy::too_many_arguments)] // One shared counter pair spans the whole file.
pub(crate) fn read_peptide_with_budget(
    node: &Node,
    identifier: &str,
    refs: &BTreeMap<String, String>,
    options: &ReadOptions,
    registry: &ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<PeptideIdentification> {
    measure_node(node, work, bytes)?;
    IdMeter { work, bytes }.text(identifier)?;
    node.check(
        &[
            "score_type",
            "higher_score_better",
            "significance_threshold",
            "RT",
            "MZ",
            "spectrum_reference",
            "identification_run_ref",
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
            sequence: AASequence::parse_with_budget(child.get("sequence")?, registry, work, bytes)?,
            score: finite(child.get("score")?)?,
            charge: number(child.get("charge")?)?,
            evidences: read_evidence(child, refs, options, work, bytes)?,
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
pub(crate) fn xml_id(value: &str) -> Result<()> {
    xml_name(value, false)
}
pub(crate) fn xml_name(value: &str, allow_colon: bool) -> Result<()> {
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

pub(crate) fn insert_meta(meta: &mut MetaInfo, key: &str, value: MetaValue) -> Result<()> {
    if let Some(existing) = meta.get(key) {
        if existing != &value {
            return Err(bad(format!("conflicting reserved metadata {key}")));
        }
    } else {
        meta.insert(key.into(), value);
    }
    Ok(())
}
pub(crate) fn write_meta(node: &mut Node, meta: &MetaInfo) -> Result<()> {
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
pub(crate) fn write_search(value: &SearchParameters, id: &str) -> Result<Node> {
    value.validate()?;
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

pub(crate) fn write_evidence(
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
pub(crate) fn write_peptide_with_budget(
    value: &PeptideIdentification,
    refs: &BTreeMap<String, String>,
    registry: &ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<Node> {
    measure_identifications(&[], std::slice::from_ref(value), work, bytes)?;
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
        if !matches!(AASequence::parse_with_budget(&sequence_text, registry, work, bytes), Ok(ref parsed) if parsed == &hit.sequence)
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
pub(crate) struct Output {
    pub(crate) bytes: Vec<u8>,
    pub(crate) records: usize,
    pub(crate) options: WriteOptions,
}
impl Output {
    pub(crate) fn append(&mut self, text: &str) -> Result<()> {
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
    fn escaped(&mut self, value: &str) -> Result<()> {
        xml_text(value)?;
        let mut start = 0;
        for (index, c) in value.char_indices() {
            let escaped = match c {
                '&' => "&amp;",
                '<' => "&lt;",
                '>' => "&gt;",
                '"' => "&quot;",
                '\n' => "&#10;",
                '\r' => "&#13;",
                '\t' => "&#9;",
                _ => continue,
            };
            self.append(&value[start..index])?;
            self.append(escaped)?;
            start = index + c.len_utf8();
        }
        self.append(&value[start..])
    }
    pub(crate) fn node(&mut self, node: &Node, depth: usize) -> Result<()> {
        if depth >= 512 {
            return Err(bad("XML output depth limit exceeded"));
        }
        self.records += 1;
        if self.records > self.options.max_records {
            return Err(bad("XML output element limit exceeded"));
        }
        xml_name(&node.name, false)?;
        for _ in 0..depth {
            self.append("  ")?;
        }
        self.append("<")?;
        self.append(&node.name)?;
        for (key, value) in &node.attrs {
            xml_name(key, true)?;
            self.append(" ")?;
            self.append(key)?;
            self.append("=\"")?;
            self.escaped(value)?;
            self.append("\"")?;
        }
        if node.children.is_empty() && node.text.is_empty() {
            self.append("/>\n")?;
        } else {
            if !node.children.is_empty()
                && !node.text.trim_matches([' ', '\t', '\r', '\n']).is_empty()
            {
                return Err(unsupported("mixed identification XML content"));
            }
            self.append(">")?;
            if node.children.is_empty() {
                self.escaped(&node.text)?;
            } else {
                self.append("\n")?;
                for child in &node.children {
                    self.node(child, depth + 1)?;
                }
                for _ in 0..depth {
                    self.append("  ")?;
                }
            }
            self.append("</")?;
            self.append(&node.name)?;
            self.append(">\n")?;
        }
        Ok(())
    }
}
#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub(crate) fn render(root: &Node, options: &WriteOptions) -> Result<Vec<u8>> {
    let mut output = Output {
        bytes: Vec::new(),
        records: 0,
        options: *options,
    };
    output.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
    output.node(root, 0)?;
    Ok(output.bytes)
}

/// Precharge copying and serializing legacy IDs before a dialect builds XML.
/// Counters are shared with parsing/registry and map geometry by the caller.
pub(crate) fn measure_identifications(
    proteins: &[ProteinIdentification],
    peptides: &[PeptideIdentification],
    remaining_work: &mut usize,
    remaining_bytes: &mut usize,
) -> Result<()> {
    let mut meter = IdMeter {
        work: remaining_work,
        bytes: remaining_bytes,
    };
    meter.slots(proteins.len().saturating_add(peptides.len()))?;
    for protein in proteins {
        meter.text(&protein.identifier)?;
        meter.text(&protein.search_engine)?;
        meter.text(&protein.search_engine_version)?;
        meter.text(&protein.score_type)?;
        if let Some(date) = &protein.date_time {
            meter.text(date)?;
        }
        meter.meta(&protein.metadata)?;
        measure_search_with_meter(&protein.search_parameters, &mut meter)?;
        meter.slots(protein.hits.len())?;
        for hit in &protein.hits {
            meter.text(&hit.accession)?;
            meter.text(&hit.sequence)?;
            meter.meta(&hit.metadata)?;
            meter.slots(hit.modifications.len())?;
        }
        meter.slots(
            protein
                .protein_groups
                .len()
                .saturating_add(protein.indistinguishable_groups.len()),
        )?;
        for group in protein
            .protein_groups
            .iter()
            .chain(&protein.indistinguishable_groups)
        {
            meter.strings(&group.accessions)?;
            meter.slots(
                group
                    .float_data_arrays
                    .len()
                    .saturating_add(group.integer_data_arrays.len())
                    .saturating_add(group.string_data_arrays.len()),
            )?;
            for array in &group.float_data_arrays {
                if array.has_description_metadata() {
                    return Err(unsupported(
                        "protein-group array description metadata or processing is not represented",
                    ));
                }
                meter.text(&array.name)?;
                meter.slots(array.data.len())?;
            }
            for array in &group.integer_data_arrays {
                if array.has_description_metadata() {
                    return Err(unsupported(
                        "protein-group array description metadata or processing is not represented",
                    ));
                }
                meter.text(&array.name)?;
                meter.slots(array.data.len())?;
            }
            for array in &group.string_data_arrays {
                if array.has_description_metadata() {
                    return Err(unsupported(
                        "protein-group array description metadata or processing is not represented",
                    ));
                }
                meter.text(&array.name)?;
                meter.strings(&array.data)?;
            }
        }
        meter.strings(&protein.primary_ms_run_paths)?;
        meter.strings(&protein.raw_ms_run_paths)?;
    }
    for peptide in peptides {
        meter.text(&peptide.identifier)?;
        meter.text(&peptide.score_type)?;
        meter.meta(&peptide.metadata)?;
        meter.slots(peptide.hits.len())?;
        for hit in &peptide.hits {
            meter.slots(hit.sequence.len())?;
            meter.add(
                hit.sequence.generation_payload_bytes()?.saturating_mul(4),
                hit.sequence.len(),
            )?;
            for annotation in hit
                .sequence
                .n_terminal_modification()
                .into_iter()
                .chain(hit.sequence.c_terminal_modification())
                .chain(
                    (0..hit.sequence.len())
                        .filter_map(|i| hit.sequence.residue_modification(i).ok().flatten()),
                )
            {
                meter.text(annotation.full_id())?;
                if let Some(record) = annotation.known() {
                    let size = record.payload_bytes()?;
                    meter.add(size.saturating_mul(4), size.saturating_mul(4))?;
                    meter.text(record.name())?;
                }
            }
            meter.meta(&hit.metadata)?;
            meter.slots(
                hit.evidences
                    .len()
                    .saturating_add(hit.peak_annotations.len())
                    .saturating_add(hit.analysis_results.len()),
            )?;
            for evidence in &hit.evidences {
                meter.text(&evidence.protein_accession)?;
            }
            for annotation in &hit.peak_annotations {
                meter.text(&annotation.annotation)?;
            }
            for analysis in &hit.analysis_results {
                meter.text(&analysis.score_type)?;
                meter.slots(analysis.sub_scores.len())?;
                for name in analysis.sub_scores.keys() {
                    meter.text(name)?;
                }
            }
        }
    }
    Ok(())
}
pub(crate) fn measure_search(
    parameters: &SearchParameters,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<()> {
    measure_search_with_meter(parameters, &mut IdMeter { work, bytes })
}
fn measure_search_with_meter(parameters: &SearchParameters, meter: &mut IdMeter<'_>) -> Result<()> {
    for text in [
        &parameters.database,
        &parameters.database_version,
        &parameters.taxonomy,
        &parameters.charges,
        &parameters.digestion_enzyme,
        &parameters.digestion_regex,
    ] {
        meter.text(text)?;
    }
    meter.strings(&parameters.fixed_modifications)?;
    meter.strings(&parameters.variable_modifications)?;
    meter.meta(&parameters.metadata)
}
pub(crate) fn clone_registry(
    registry: &ModificationsDB,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<ModificationsDB> {
    registry.clone_with_budget(work, bytes)
}
struct IdMeter<'a> {
    work: &'a mut usize,
    bytes: &'a mut usize,
}
impl IdMeter<'_> {
    fn add(&mut self, bytes: usize, work: usize) -> Result<()> {
        *self.work = self
            .work
            .checked_sub(work)
            .ok_or_else(|| bad("identification XML work limit exceeded"))?;
        *self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| bad("identification XML payload limit exceeded"))?;
        Ok(())
    }
    fn slots(&mut self, n: usize) -> Result<()> {
        self.add(n.saturating_mul(1024), n)
    }
    fn text(&mut self, text: &str) -> Result<()> {
        self.add(
            text.len().saturating_mul(16).saturating_add(128),
            text.len().saturating_mul(4).saturating_add(1),
        )
    }
    fn strings(&mut self, values: &[String]) -> Result<()> {
        self.slots(values.len())?;
        for value in values {
            self.text(value)?;
        }
        Ok(())
    }
    fn meta(&mut self, metadata: &MetaInfo) -> Result<()> {
        self.slots(metadata.len())?;
        for (name, value) in metadata {
            if value.unit().is_some() {
                return Err(unsupported("identification XML metadata units"));
            }
            self.text(name)?;
            match value.data() {
                MetaValueData::String(value) => self.text(value)?,
                MetaValueData::StringList(values) => self.strings(values)?,
                MetaValueData::IntegerList(values) => self.slots(values.len())?,
                MetaValueData::FloatList(values) => self.slots(values.len())?,
                _ => self.slots(1)?,
            }
        }
        Ok(())
    }
}

/// Read only through the requested root-child opening tag. Lexical scanning is
/// ASCII-compatible for UTF-8/Latin-1 and code-unit-aware for UTF-16; the shared
/// parser subsequently validates the complete retained prefix and its encoding.
fn read_xml_prefix(mut input: impl BufRead, limit: usize, stop: &str) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    let mut first = [0u8; 4];
    let mut initial = 0;
    while initial < 4 {
        match input.read(&mut first[initial..initial + 1])? {
            0 => break,
            _ => {
                initial += 1;
                if initial > limit {
                    return Err(bad("XML prefix byte limit exceeded"));
                }
            }
        }
    }
    bytes.extend_from_slice(&first[..initial]);
    let (width, little, mut position) = if bytes.starts_with(&[0xff, 0xfe]) {
        (2, true, 2)
    } else if bytes.starts_with(&[0xfe, 0xff]) {
        (2, false, 2)
    } else if bytes.starts_with(&[b'<', 0, b'?', 0]) {
        (2, true, 0)
    } else if bytes.starts_with(&[0, b'<', 0, b'?']) {
        (2, false, 0)
    } else {
        (1, false, 0)
    };
    let mut token = Vec::new();
    let mut quote = None;
    let mut depth = 0usize;
    loop {
        while bytes.len() < position + width {
            if bytes.len() == limit {
                return Err(bad("XML prefix byte limit exceeded"));
            }
            let mut next = [0u8; 1];
            if input.read(&mut next)? == 0 {
                return Ok(bytes);
            }
            bytes.push(next[0]);
        }
        let value = if width == 1 {
            u16::from(bytes[position])
        } else if little {
            u16::from_le_bytes([bytes[position], bytes[position + 1]])
        } else {
            u16::from_be_bytes([bytes[position], bytes[position + 1]])
        };
        position += width;
        let c = u8::try_from(value).unwrap_or(0xff);
        if token.is_empty() {
            if c == b'<' {
                token.push(c);
            }
            continue;
        }
        token.push(c);
        if token.starts_with(b"<!--") {
            if token.ends_with(b"-->") {
                token.clear();
            }
            continue;
        }
        if token.starts_with(b"<![CDATA[") {
            if token.ends_with(b"]]>") {
                token.clear();
            }
            continue;
        }
        if token.starts_with(b"<?") {
            if token.ends_with(b"?>") {
                token.clear();
            }
            continue;
        }
        if let Some(q) = quote {
            if c == q {
                quote = None;
            }
            continue;
        }
        if matches!(c, b'\'' | b'"') {
            quote = Some(c);
            continue;
        }
        if c != b'>' {
            continue;
        }
        if token.starts_with(b"<!") {
            return Err(unsupported("DTD in identification XML"));
        }
        if token.starts_with(b"</") {
            depth = depth.saturating_sub(1);
        } else {
            let end = token[1..]
                .iter()
                .position(|b| matches!(b, b' ' | b'\t' | b'\r' | b'\n' | b'/' | b'>'))
                .map_or(token.len(), |i| i + 1);
            if depth == 1 && &token[1..end] == stop.as_bytes() {
                return Ok(bytes);
            }
            if !token.ends_with(b"/>") {
                depth = depth.saturating_add(1);
            }
        }
        token.clear();
    }
}

/// Cumulative typed-conversion allowance, including list slot amplification.
pub(crate) fn measure_node(node: &Node, work: &mut usize, bytes: &mut usize) -> Result<()> {
    let mut meter = IdMeter { work, bytes };
    fn visit(node: &Node, meter: &mut IdMeter<'_>, depth: usize) -> Result<()> {
        if depth > 512 {
            return Err(bad("XML conversion depth limit exceeded"));
        }
        meter.slots(1)?;
        meter.text(&node.name)?;
        meter.text(&node.text)?;
        for (key, value) in &node.attrs {
            meter.text(key)?;
            meter.text(value)?;
            meter.add(
                value.len().saturating_mul(128),
                value.len().saturating_mul(8),
            )?;
        }
        for child in &node.children {
            visit(child, meter, depth + 1)?;
        }
        Ok(())
    }
    visit(node, &mut meter, 0)
}

#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
pub(crate) fn measure_metadata(
    metadata: &MetaInfo,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<()> {
    IdMeter { work, bytes }.meta(metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn resolved_reference_and_run_identifier_copies_share_the_remaining_budget() {
        let refs = [("PH_0".to_owned(), "A".repeat(50_000))].into();
        let mut node = Node::new("PeptideHit");
        node.attr("protein_refs", "PH_0");
        assert_eq!(
            read_evidence(
                &node,
                &refs,
                &ReadOptions::default(),
                &mut 1_000_000,
                &mut 1_000_000
            )
            .unwrap()
            .len(),
            1
        );
        node.attr("protein_refs", "PH_0 PH_0");
        let before = node.clone();
        assert!(
            read_evidence(
                &node,
                &refs,
                &ReadOptions::default(),
                &mut 1_000_000,
                &mut 1_000_000
            )
            .is_err()
        );
        assert_eq!(node.attrs, before.attrs);
        let mut peptide = Node::new("PeptideIdentification");
        peptide.attr("score_type", "score");
        peptide.attr("higher_score_better", true);
        assert!(
            read_peptide_with_budget(
                &peptide,
                &"R".repeat(50_000),
                &BTreeMap::new(),
                &ReadOptions::default(),
                &ModificationsDB::default(),
                &mut 1_000_000,
                &mut 100_000
            )
            .is_err()
        );
    }

    fn root_text(xml: &str) -> Result<String> {
        parse_xml_with_budget(
            xml.as_bytes(),
            &ReadOptions::default(),
            8,
            None,
            &mut 1_000_000,
            &mut 1_000_000,
        )
        .map(|root| root.text)
    }

    /// References in element text follow the XML 1.0 `CharRef` and predefined
    /// entity rules: a lowercase `x` for hexadecimal, no sign, and no NUL.
    #[test]
    fn text_references_follow_the_xml_1_0_grammar() {
        let mut wrong = Vec::new();
        for (xml, expected) in [
            ("<r>1&#46;5</r>", "1.5"),
            ("<r>1&#x2E;5</r>", "1.5"),
            ("<r>1&#x2e;5</r>", "1.5"),
            ("<r>&#0046;&#x0002E;</r>", ".."),
            ("<r>&amp;&lt;&gt;&apos;&quot;</r>", "&<>'\""),
            ("<r>&#9;&#xA;&#13;&#x10FFFF;</r>", "\t\n\r\u{10FFFF}"),
        ] {
            match root_text(xml) {
                Ok(text) if text == expected => {}
                other => wrong.push(format!("{xml}: {other:?}")),
            }
        }
        for xml in [
            "<r>1&#X2E;5</r>",
            "<r>1&#+46;5</r>",
            "<r>1&#x+2E;5</r>",
            "<r>1&#-46;5</r>",
            "<r>&#0;</r>",
            "<r>&#x0;</r>",
            "<r>&#;</r>",
            "<r>&#x;</r>",
            "<r>&#xD800;</r>",
            "<r>&#x110000;</r>",
            "<r>&#1;</r>",
            "<r>&#xFFFE;</r>",
            "<r/>&#46;",
        ] {
            match root_text(xml) {
                Err(Error::Parse { .. }) => {}
                other => wrong.push(format!("{xml}: {other:?}")),
            }
        }
        for xml in ["<r>&bogus;</r>", "<r>&x2E;</r>", "<r>&AMP;</r>"] {
            match root_text(xml) {
                Err(Error::Unsupported(_)) => {}
                other => wrong.push(format!("{xml}: {other:?}")),
            }
        }
        assert!(wrong.is_empty(), "{wrong:#?}");
    }
}

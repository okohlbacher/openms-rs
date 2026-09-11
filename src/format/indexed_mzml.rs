// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded indexed-mzML footer discovery and spectrum/chromatogram offsets.
//!
//! This reads the plain-file index only; it does not validate the mzML payload,
//! checksum, or the XML element located at each reported offset. Source contracts
//! are from Core SDK `54a232fe2cae9c590d5c997fa49d20e7769860fb`.

use crate::{Error, Result};
use quick_xml::{Reader, events::Event};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

/// Native ID and byte offset, in the order stored in an index section.
pub type OffsetVector = Vec<(String, u64)>;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct IndexOffsets {
    pub spectra: OffsetVector,
    pub chromatograms: OffsetVector,
}

/// Limits apply before suffix allocation and cumulatively to parsed offsets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IndexReadLimits {
    pub max_footer_bytes: usize,
    pub max_index_bytes: usize,
    pub max_offsets: usize,
    pub max_id_bytes: usize,
    pub max_depth: usize,
    pub max_attributes: usize,
}

impl Default for IndexReadLimits {
    fn default() -> Self {
        Self {
            max_footer_bytes: 1 << 20,
            max_index_bytes: 16 << 20,
            max_offsets: 1_000_000,
            max_id_bytes: 65_536,
            max_depth: 64,
            max_attributes: 64,
        }
    }
}

/// The source decoder's two public operations, with owned, atomic results.
#[derive(Clone, Copy, Debug, Default)]
pub struct IndexedMzMLDecoder {
    pub limits: IndexReadLimits,
}

impl IndexedMzMLDecoder {
    /// Search the final 1023 bytes, as in the source default.
    pub fn find_index_list_offset(&self, path: impl AsRef<Path>) -> Result<Option<u64>> {
        self.find_index_list_offset_with_buffer_size(path, 1023)
    }

    /// A small file searches its whole content. Zero bytes gives no match.
    pub fn find_index_list_offset_with_buffer_size(
        &self,
        path: impl AsRef<Path>,
        buffer_size: usize,
    ) -> Result<Option<u64>> {
        if buffer_size > self.limits.max_footer_bytes {
            return Err(bad("footer search exceeds byte limit"));
        }
        let mut file = File::open(path)?;
        let length = file.seek(SeekFrom::End(0))?;
        let amount = length.min(buffer_size as u64);
        file.seek(SeekFrom::Start(length - amount))?;
        let bytes = read_exact_suffix(&mut file, amount, self.limits.max_footer_bytes)?;
        probe(&bytes)
    }

    /// Read the suffix beginning at `index_offset` and decode its index sections.
    /// Missing sections are empty; repeated sections replace earlier sections.
    /// Duplicate IDs and index order are preserved. Invalid input returns no data.
    pub fn parse_offsets(&self, path: impl AsRef<Path>, index_offset: u64) -> Result<IndexOffsets> {
        let mut file = File::open(path)?;
        let length = file.seek(SeekFrom::End(0))?;
        let amount = length
            .checked_sub(index_offset)
            .ok_or_else(|| bad("index offset exceeds file length"))?;
        file.seek(SeekFrom::Start(index_offset))?;
        let bytes = read_exact_suffix(&mut file, amount, self.limits.max_index_bytes)?;
        parse(&bytes, self.limits)
    }
}

/// Source `MzMLFile::hasIndex`: a footer-offset probe, not index validation.
pub fn has_index(path: impl AsRef<Path>) -> Result<bool> {
    Ok(IndexedMzMLDecoder::default()
        .find_index_list_offset(path)?
        .is_some())
}

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn read_exact_suffix(file: &mut File, amount: u64, limit: usize) -> Result<Vec<u8>> {
    let amount =
        usize::try_from(amount).map_err(|_| bad("suffix size exceeds addressable memory"))?;
    if amount > limit {
        return Err(bad("index suffix exceeds byte limit"));
    }
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(amount)
        .map_err(|_| bad("cannot allocate index suffix"))?;
    file.take(amount as u64).read_to_end(&mut bytes)?;
    if bytes.len() != amount {
        return Err(bad("file truncated while reading index suffix"));
    }
    Ok(bytes)
}

fn whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | 9..=13)
}

fn position(value: &str) -> Result<u64> {
    value
        .parse::<u64>()
        .ok()
        .filter(|&v| v <= i64::MAX as u64)
        .ok_or_else(|| bad("invalid nonnegative 63-bit file offset"))
}

// Linear scan equivalent to the source's broad footer regex. In particular,
// prefixes are accepted and the first matching tag with no digits stops probing.
fn probe(bytes: &[u8]) -> Result<Option<u64>> {
    let bytes = &bytes[..bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len())];
    let mut cursor = 0;
    while cursor < bytes.len() {
        let Some(start) = bytes[cursor..].iter().position(|&b| b == b'<') else {
            break;
        };
        let start = cursor + start + 1;
        let Some(end) = bytes[start..].iter().position(|&b| b == b'>' || b == b'/') else {
            break;
        };
        let end = start + end;
        cursor = end + 1;
        if bytes[end] != b'>' {
            continue;
        }
        let header = &bytes[start..end];
        let trimmed = header
            .iter()
            .rposition(|&b| !whitespace(b))
            .map_or(0, |i| i + 1);
        if !header[..trimmed].ends_with(b"indexListOffset") {
            continue;
        }
        while cursor < bytes.len() && whitespace(bytes[cursor]) {
            cursor += 1;
        }
        let first = cursor;
        while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
            cursor += 1;
        }
        return if first == cursor {
            Ok(None)
        } else {
            position(std::str::from_utf8(&bytes[first..cursor]).expect("ASCII digits")).map(Some)
        };
    }
    Ok(None)
}

fn xml_text(text: &str) -> Result<()> {
    if text.chars().any(|c| !matches!(c, '\t' | '\n' | '\r' | '\u{20}'..='\u{d7ff}' | '\u{e000}'..='\u{fffd}' | '\u{10000}'..='\u{10ffff}')) {
        return Err(bad("invalid XML character in index"));
    }
    Ok(())
}

fn xml_name(value: &str) -> Result<()> {
    let start = |c: char| {
        c.is_ascii_alphabetic()
            || matches!(c, ':' | '_'
        | '\u{c0}'..='\u{d6}' | '\u{d8}'..='\u{f6}' | '\u{f8}'..='\u{2ff}'
        | '\u{370}'..='\u{37d}' | '\u{37f}'..='\u{1fff}' | '\u{200c}'..='\u{200d}'
        | '\u{2070}'..='\u{218f}' | '\u{2c00}'..='\u{2fef}' | '\u{3001}'..='\u{d7ff}'
        | '\u{f900}'..='\u{fdcf}' | '\u{fdf0}'..='\u{fffd}' | '\u{10000}'..='\u{effff}')
    };
    let mut chars = value.chars();
    if !chars.next().is_some_and(start) || !chars.all(|c| {
        start(c)
            || c.is_ascii_digit()
            || matches!(c, '-' | '.' | '\u{b7}' | '\u{300}'..='\u{36f}' | '\u{203f}'..='\u{2040}')
    }) {
        return Err(bad("invalid XML name"));
    }
    Ok(())
}

fn xml_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n')
}

fn attribute_spacing(element: &quick_xml::events::BytesStart<'_>) -> Result<()> {
    let tail = &element.as_ref()[element.name().as_ref().len()..];
    if tail.first().is_some_and(|&b| !xml_space(b)) {
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
            if !xml_space(byte) {
                return Err(bad("XML attributes require whitespace separators"));
            }
            closed = false;
        } else if matches!(byte, b'\'' | b'"') {
            quote = Some(byte);
        }
    }
    Ok(())
}

fn parse(bytes: &[u8], limits: IndexReadLimits) -> Result<IndexOffsets> {
    let suffix = std::str::from_utf8(bytes).map_err(|_| bad("index suffix must be UTF-8"))?;
    xml_text(suffix)?;
    let xml = format!("<indexedmzML>{suffix}");
    let mut reader = Reader::from_str(&xml);
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().check_comments = true;
    let mut stack = Vec::<String>::new();
    let mut output = IndexOffsets::default();
    let mut index_count = 0;
    let mut offset_count = 0usize;
    let mut category = None;
    let mut pending = None::<(String, String)>;
    let mut closed = false;
    loop {
        match reader.read_event().map_err(|e| bad(e.to_string()))? {
            Event::Start(element) => {
                if closed || stack.len() >= limits.max_depth {
                    return Err(bad("index XML depth or root limit exceeded"));
                }
                let name = std::str::from_utf8(element.name().as_ref())
                    .map_err(|_| bad("invalid XML name"))?
                    .to_owned();
                xml_name(&name)?;
                attribute_spacing(&element)?;
                let parent = stack.last().map(String::as_str).unwrap_or("");
                if pending.is_some() {
                    return Err(bad("nested element in offset value"));
                }
                let mut section_name = None;
                let mut id = String::new();
                let mut attribute_names = BTreeSet::new();
                for attribute in element.attributes().with_checks(false) {
                    if attribute_names.len() >= limits.max_attributes {
                        return Err(bad("XML attribute count limit exceeded"));
                    }
                    let attribute = attribute.map_err(|e| bad(e.to_string()))?;
                    xml_name(
                        std::str::from_utf8(attribute.key.as_ref())
                            .map_err(|_| bad("invalid attribute name"))?,
                    )?;
                    if !attribute_names.insert(attribute.key.0) {
                        return Err(bad("duplicate XML attribute"));
                    }
                    let raw = std::str::from_utf8(&attribute.value)
                        .map_err(|_| bad("invalid UTF-8 attribute"))?;
                    if raw.contains('<') {
                        return Err(bad("unescaped less-than in XML attribute"));
                    }
                    let value = quick_xml::escape::unescape(
                        &raw.replace("\r\n", "\n").replace(['\t', '\r', '\n'], " "),
                    )
                    .map_err(|e| bad(e.to_string()))?
                    .into_owned();
                    xml_text(&value)?;
                    match attribute.key.as_ref() {
                        b"name" => section_name = Some(value),
                        b"idRef" => {
                            if value.len() > limits.max_id_bytes {
                                return Err(bad("offset ID exceeds byte limit"));
                            }
                            id = value;
                        }
                        _ => {}
                    }
                }
                match (name.as_str(), parent) {
                    ("indexedmzML", "") => {}
                    ("indexList", "indexedmzML") => {
                        index_count += 1;
                        if index_count != 1 {
                            return Err(bad("multiple indexList elements"));
                        }
                    }
                    ("index", "indexList") => {
                        category = Some(match section_name.as_deref() {
                            Some("spectrum") => {
                                output.spectra.clear();
                                true
                            }
                            Some("chromatogram") => {
                                output.chromatograms.clear();
                                false
                            }
                            _ => return Err(bad("unknown index section name")),
                        });
                    }
                    ("offset", "index") => {
                        offset_count = offset_count
                            .checked_add(1)
                            .ok_or_else(|| bad("offset count overflow"))?;
                        if offset_count > limits.max_offsets {
                            return Err(bad("offset count limit exceeded"));
                        }
                        pending = Some((id, String::new()));
                    }
                    (_, "indexList" | "index") | ("indexList" | "index" | "offset", _) => {
                        return Err(bad("misplaced index element"));
                    }
                    _ => {}
                }
                stack.push(name);
            }
            Event::End(element) => {
                let name = std::str::from_utf8(element.name().as_ref())
                    .map_err(|_| bad("invalid XML end name"))?
                    .to_owned();
                if stack.pop().as_deref() != Some(name.as_str()) {
                    return Err(bad("mismatched index XML end"));
                }
                if name == "offset" {
                    let (id, value) = pending.take().ok_or_else(|| bad("missing offset start"))?;
                    let offset = position(&value)?;
                    if category.ok_or_else(|| bad("offset outside index section"))? {
                        output.spectra.push((id, offset));
                    } else {
                        output.chromatograms.push((id, offset));
                    }
                } else if name == "index" {
                    category = None;
                }
                if stack.is_empty() {
                    closed = true;
                }
            }
            Event::Text(text) => {
                let value = text.decode().map_err(|e| bad(e.to_string()))?;
                if value.contains("]]>") {
                    return Err(bad("invalid XML text delimiter"));
                }
                append_text(&mut pending, &value, closed)?;
            }
            Event::CData(text) => {
                if closed {
                    return Err(bad("CDATA after index XML root"));
                }
                append_text(
                    &mut pending,
                    &text.decode().map_err(|e| bad(e.to_string()))?,
                    false,
                )?;
            }
            Event::GeneralRef(reference) => {
                if closed {
                    return Err(bad("entity after index XML root"));
                }
                let value = if let Some(c) = reference
                    .resolve_char_ref()
                    .map_err(|e| bad(e.to_string()))?
                {
                    c.to_string()
                } else {
                    let name = reference.decode().map_err(|e| bad(e.to_string()))?;
                    match name.as_ref() {
                        "amp" => "&",
                        "lt" => "<",
                        "gt" => ">",
                        "apos" => "'",
                        "quot" => "\"",
                        _ => return Err(bad("unknown XML entity")),
                    }
                    .to_owned()
                };
                xml_text(&value)?;
                append_text(&mut pending, &value, closed)?;
            }
            Event::DocType(_) | Event::Decl(_) => {
                return Err(Error::Unsupported(
                    "DTD/XML declaration inside index suffix".into(),
                ));
            }
            Event::Eof => break,
            Event::PI(instruction) => {
                let target = std::str::from_utf8(instruction.target())
                    .map_err(|_| bad("invalid processing instruction"))?;
                xml_name(target)?;
                if target.eq_ignore_ascii_case("xml") {
                    return Err(bad("reserved XML processing instruction target"));
                }
            }
            Event::Comment(_) => {}
            Event::Empty(_) => unreachable!("empty elements are expanded"),
        }
    }
    if !closed || !stack.is_empty() || index_count != 1 {
        return Err(bad("incomplete indexed-mzML suffix"));
    }
    Ok(output)
}

fn append_text(pending: &mut Option<(String, String)>, value: &str, closed: bool) -> Result<()> {
    if closed && !value.bytes().all(xml_space) {
        return Err(bad("text after index XML root"));
    }
    if let Some((_, text)) = pending {
        if text.len().saturating_add(value.len()) > 64 {
            return Err(bad("offset number exceeds text limit"));
        }
        text.push_str(value);
    }
    Ok(())
}

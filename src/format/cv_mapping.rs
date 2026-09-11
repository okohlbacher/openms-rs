// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Bounded source CV mapping XML loading; semantic validation is separate.

use crate::data_structures::{
    CVMappingRule, CVMappingTerm, CVMappings, CVReference, CombinationsLogic, RequirementLevel,
};
use crate::{Error, Result};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::{io::BufRead, mem::size_of, path::Path};

#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    pub strip_namespaces: bool,
    /// Maximum compressed-stream decoded input bytes and normalized UTF-8 bytes.
    pub max_input_bytes: usize,
    pub max_depth: usize,
    pub max_elements: usize,
    pub max_references: usize,
    pub max_rules: usize,
    pub max_terms: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            strip_namespaces: false,
            max_input_bytes: 16 * 1024 * 1024,
            max_depth: 128,
            max_elements: 1_000_000,
            max_references: 100_000,
            max_rules: 100_000,
            max_terms: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 128 * 1024 * 1024,
        }
    }
}

/// Each operation has local parser state. Reusing this loader after an error is safe.
#[derive(Clone, Copy, Debug, Default)]
pub struct CVMappingFile {
    pub options: ReadOptions,
}
impl CVMappingFile {
    pub fn load(&self, path: impl AsRef<Path>) -> Result<CVMappings> {
        self.read(super::path_io::open(path.as_ref())?)
    }
    pub fn read(&self, input: impl BufRead) -> Result<CVMappings> {
        let mut output = CVMappings::default();
        self.read_into(input, &mut output)?;
        Ok(output)
    }
    pub fn load_into(&self, path: impl AsRef<Path>, output: &mut CVMappings) -> Result<()> {
        self.read_into(super::path_io::open(path.as_ref())?, output)
    }
    /// Atomically append parsed references and replace rules, matching successful
    /// source loading into an existing destination. An error leaves output intact.
    pub fn read_into(&self, input: impl BufRead, output: &mut CVMappings) -> Result<()> {
        let mut meter = Meter {
            work: self.options.max_work,
            bytes: self.options.max_bytes,
        };
        let text = document(input, &self.options, &mut meter)?;
        let (references, rules) = parse(&text, &self.options, &mut meter)?;
        let count = output
            .cv_references()
            .len()
            .checked_add(references.len())
            .ok_or_else(limit)?;
        cap(count, self.options.max_references)?;
        // Only old references survive. Existing rules need not be cloned.
        meter.spend(count, 0)?;
        meter.slots::<CVReference>(count.saturating_mul(2))?;
        meter.tree::<String>(count)?;
        let mut joined = Vec::with_capacity(count);
        for r in output.cv_references() {
            joined.push(CVReference {
                name: meter.copy(&r.name)?,
                identifier: meter.copy(&r.identifier)?,
            });
        }
        joined.extend(references);
        // set_cv_references clones index keys even for a repeated identifier.
        // Bound all comparisons by the maximum key length before building it.
        let max_key = joined.iter().map(|r| r.identifier.len()).max().unwrap_or(0);
        let comparisons = (usize::BITS - count.max(1).leading_zeros()) as usize * 12;
        for r in &joined {
            meter.spend(
                comparisons.saturating_mul(max_key.saturating_add(1)),
                r.identifier.len(),
            )?;
        }
        let mut next = CVMappings::default();
        next.mapping_rules = rules;
        next.set_cv_references(joined);
        *output = next;
        Ok(())
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn limit() -> Error {
    bad("CV mapping resource limit exceeded")
}
fn cap(value: usize, maximum: usize) -> Result<()> {
    if value > maximum {
        Err(limit())
    } else {
        Ok(())
    }
}
struct Meter {
    work: usize,
    bytes: usize,
}
impl Meter {
    fn spend(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.work = self.work.checked_sub(work).ok_or_else(limit)?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(limit)?;
        Ok(())
    }
    fn slots<T>(&mut self, n: usize) -> Result<()> {
        self.spend(n, n.checked_mul(size_of::<T>()).ok_or_else(limit)?)
    }
    fn tree<T>(&mut self, n: usize) -> Result<()> {
        if n == 0 {
            return Ok(());
        }
        self.spend(
            n,
            512usize
                .checked_add(
                    n.saturating_mul(3)
                        .max(11)
                        .saturating_mul(size_of::<T>() + 4 * size_of::<usize>()),
                )
                .ok_or_else(limit)?,
        )
    }
    fn copy(&mut self, s: &str) -> Result<String> {
        self.spend(s.len().saturating_add(1), s.len())?;
        Ok(s.into())
    }
    fn push<T>(&mut self, v: &mut Vec<T>, item: T) -> Result<()> {
        if v.len() == v.capacity() {
            let next = v.capacity().saturating_mul(2).max(4);
            self.slots::<T>(next)?;
            self.spend(v.len(), 0)?;
            v.reserve_exact(next - v.len());
        }
        v.push(item);
        Ok(())
    }
}
fn xml_char(c: char) -> bool {
    matches!(c, '\t' | '\n' | '\r')
        || ('\u{20}'..='\u{d7ff}').contains(&c)
        || ('\u{e000}'..='\u{fffd}').contains(&c)
        || c >= '\u{10000}'
}
fn xml_text(s: &str) -> Result<()> {
    if s.chars().all(xml_char) {
        Ok(())
    } else {
        Err(bad("invalid XML 1.0 character"))
    }
}
fn space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\r' | b'\n')
}
fn name_start(c: char) -> bool {
    matches!(c,':'|'_'|'A'..='Z'|'a'..='z'|'\u{c0}'..='\u{d6}'|'\u{d8}'..='\u{f6}'|'\u{f8}'..='\u{2ff}'|'\u{370}'..='\u{37d}'|'\u{37f}'..='\u{1fff}'|'\u{200c}'..='\u{200d}'|'\u{2070}'..='\u{218f}'|'\u{2c00}'..='\u{2fef}'|'\u{3001}'..='\u{d7ff}'|'\u{f900}'..='\u{fdcf}'|'\u{fdf0}'..='\u{fffd}'|'\u{10000}'..='\u{effff}')
}
fn name(s: &str) -> Result<()> {
    let mut cs = s.chars();
    if !cs.next().is_some_and(name_start) || !cs.all(|c| {
        name_start(c)
            || matches!(c,'-'|'.'|'0'..='9'|'\u{b7}'|'\u{300}'..='\u{36f}'|'\u{203f}'..='\u{2040}')
    }) {
        return Err(bad("invalid XML name"));
    }
    Ok(())
}
fn spacing(element: &BytesStart<'_>) -> Result<()> {
    let tail = element.attributes_raw();
    if tail.first().is_some_and(|b| !space(*b)) {
        return Err(bad("missing XML attribute separator"));
    }
    let (mut quote, mut closed) = (None, false);
    for &b in tail {
        if let Some(q) = quote {
            if b == q {
                quote = None;
                closed = true;
            }
        } else {
            if closed && !space(b) {
                return Err(bad("missing XML attribute separator"));
            }
            closed = false;
            if b == b'\'' || b == b'"' {
                quote = Some(b);
            }
        }
    }
    Ok(())
}
type Attributes = Vec<(String, String)>;
fn attrs(element: &BytesStart<'_>, meter: &mut Meter) -> Result<Attributes> {
    spacing(element)?;
    let mut out = Vec::new();
    for attr in element.attributes().with_checks(false) {
        let attr = attr.map_err(|e| bad(e.to_string()))?;
        let key =
            std::str::from_utf8(attr.key.as_ref()).map_err(|_| bad("invalid attribute UTF-8"))?;
        name(key)?;
        for (previous, _) in &out {
            let previous: &String = previous;
            meter.spend(
                previous.len().saturating_add(key.len()).saturating_add(1),
                0,
            )?;
            if previous == key {
                return Err(bad("duplicate XML attribute"));
            }
        }
        if attr.value.contains(&b'<') {
            return Err(bad("raw '<' in XML attribute"));
        }
        let raw = std::str::from_utf8(&attr.value).map_err(|_| bad("invalid attribute UTF-8"))?;
        meter.spend(raw.len().saturating_mul(3), raw.len().saturating_mul(3))?;
        // Normalize literal whitespace before entity expansion; &#10; stays LF.
        let mut normalized = String::with_capacity(raw.len());
        for c in raw.chars() {
            normalized.push(if matches!(c, '\t' | '\n' | '\r') {
                ' '
            } else {
                c
            });
        }
        let value = quick_xml::escape::unescape(&normalized).map_err(|e| bad(e.to_string()))?;
        xml_text(&value)?;
        let pair = (meter.copy(key)?, value.into_owned());
        meter.push(&mut out, pair)?;
    }
    Ok(out)
}
fn get<'a>(a: &'a Attributes, key: &str) -> Option<&'a str> {
    a.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
}
fn required<'a>(a: &'a Attributes, key: &str) -> Result<&'a str> {
    get(a, key).ok_or_else(|| bad(format!("missing CV mapping attribute {key}")))
}
fn boolean(s: &str) -> Result<bool> {
    match s {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(bad("CV mapping boolean requires true or false")),
    }
}

// Read complete bounded input before borrowed XML events. A failed short read
// propagates directly; no uninitialized compression sniff or parser reuse state.
fn document(mut input: impl BufRead, o: &ReadOptions, m: &mut Meter) -> Result<String> {
    let mut bytes = Vec::new();
    loop {
        let buf = input.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        let remaining = o.max_input_bytes.saturating_sub(bytes.len());
        let n = buf.len().min(8192).min(remaining.saturating_add(1));
        m.spend(n, 0)?;
        if n > remaining {
            return Err(limit());
        }
        let wanted = bytes.len().checked_add(n).ok_or_else(limit)?;
        if wanted > bytes.capacity() {
            let capacity = bytes
                .capacity()
                .saturating_mul(2)
                .max(wanted)
                .min(o.max_input_bytes);
            m.spend(bytes.len(), capacity)?;
            bytes.reserve_exact(capacity - bytes.len());
        }
        bytes.extend_from_slice(&buf[..n]);
        input.consume(n);
    }
    let utf16 = if bytes.starts_with(&[255, 254]) {
        Some((true, 2))
    } else if bytes.starts_with(&[254, 255]) {
        Some((false, 2))
    } else if bytes.starts_with(&[b'<', 0, b'?', 0]) {
        Some((true, 0))
    } else if bytes.starts_with(&[0, b'<', 0, b'?']) {
        Some((false, 0))
    } else {
        None
    };
    m.spend(bytes.len().saturating_mul(4), bytes.len().saturating_mul(3))?;
    let mut text = if let Some((little, offset)) = utf16 {
        if (bytes.len() - offset) % 2 != 0 {
            return Err(bad("odd UTF-16 byte count"));
        }
        let mut result = String::with_capacity(bytes.len().saturating_mul(2));
        for c in char::decode_utf16(bytes[offset..].chunks_exact(2).map(|p| {
            if little {
                u16::from_le_bytes([p[0], p[1]])
            } else {
                u16::from_be_bytes([p[0], p[1]])
            }
        })) {
            let c = c.map_err(|_| bad("invalid UTF-16 XML"))?;
            cap(result.len().saturating_add(c.len_utf8()), o.max_input_bytes)?;
            result.push(c);
        }
        result
    } else {
        let content = bytes.strip_prefix(&[239, 187, 191]).unwrap_or(&bytes);
        std::str::from_utf8(content)
            .map_err(|_| {
                Error::Unsupported(
                    "CV mapping input requires UTF-8, UTF-16 or ASCII-compatible bytes".into(),
                )
            })?
            .into()
    };
    // Validate a declaration before accepting its byte encoding. The parser does
    // not reinterpret UTF-16-declared text after explicit decoding to UTF-8.
    let declared = declaration(&text, m)?;
    if let Some((little, _)) = utf16 {
        if declared.as_deref().is_some_and(|v| {
            !matches!(v, "utf-16" | "utf16") && v != if little { "utf-16le" } else { "utf-16be" }
        }) {
            return Err(bad("XML declaration conflicts with UTF-16 bytes"));
        }
    } else if let Some(ref encoding) = declared {
        match encoding.as_str() {
            "utf-8" | "utf8" => {}
            "us-ascii" | "ascii" | "iso-8859-1" | "iso8859-1" | "latin1" if text.is_ascii() => {}
            _ => {
                return Err(Error::Unsupported(format!(
                    "CV mapping XML encoding {encoding}"
                )));
            }
        }
        if bytes.starts_with(&[239, 187, 191]) && !matches!(encoding.as_str(), "utf-8" | "utf8") {
            return Err(bad("XML declaration conflicts with UTF-8 BOM"));
        }
    }
    xml_text(&text)?;
    m.spend(text.len().saturating_mul(2), text.len().saturating_mul(2))?;
    text = text.replace("\r\n", "\n").replace('\r', "\n");
    cap(text.len(), o.max_input_bytes)?;
    Ok(text)
}
fn declaration(text: &str, m: &mut Meter) -> Result<Option<String>> {
    let mut reader = Reader::from_str(text);
    let Event::Decl(d) = reader.read_event().map_err(|e| bad(e.to_string()))? else {
        return Ok(None);
    };
    let raw = std::str::from_utf8(d.as_ref()).map_err(|_| bad("invalid XML declaration"))?;
    let e = BytesStart::from_content(raw, 3);
    spacing(&e)?;
    let mut phase = 0;
    let mut encoding = None;
    for attr in e.attributes() {
        let a = attr.map_err(|e| bad(e.to_string()))?;
        match a.key.as_ref() {
            b"version" if phase == 0 && a.value.as_ref() == b"1.0" => phase = 1,
            b"encoding"
                if phase == 1
                    && a.value.first().is_some_and(u8::is_ascii_alphabetic)
                    && a.value
                        .iter()
                        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-')) =>
            {
                phase = 2;
                m.spend(a.value.len(), a.value.len())?;
                encoding = Some(
                    std::str::from_utf8(&a.value)
                        .map_err(|_| bad("invalid encoding name"))?
                        .to_ascii_lowercase(),
                );
            }
            b"standalone"
                if matches!(phase, 1 | 2) && matches!(a.value.as_ref(), b"yes" | b"no") =>
            {
                phase = 3
            }
            _ => return Err(bad("invalid XML declaration attributes or order")),
        }
    }
    if phase == 0 {
        return Err(bad("XML declaration needs version 1.0"));
    }
    Ok(encoding)
}

/// Correct source namespace stripping (CPP-032). The source intentionally strips
/// only element_path, never scope_path. Slash compaction remains source behavior.
fn strip_path(path: &str, m: &mut Meter) -> Result<String> {
    m.spend(path.len().saturating_mul(3), path.len().saturating_add(1))?;
    let mut out = String::with_capacity(path.len().saturating_add(1));
    for segment in path.split('/').filter(|s| !s.is_empty()) {
        out.push('/');
        if let Some((prefix, local)) = segment.split_once(':') {
            if local.contains(':') {
                return Err(bad("multiple namespace colons in mapping path"));
            }
            if prefix.starts_with('@') {
                out.push('@');
            }
            out.push_str(local);
        } else {
            out.push_str(segment);
        }
    }
    Ok(out)
}
struct Draft {
    references: Vec<CVReference>,
    rules: Vec<CVMappingRule>,
    rule: CVMappingRule,
    terms: usize,
}
impl Draft {
    fn start(&mut self, tag: &str, a: &Attributes, o: &ReadOptions, m: &mut Meter) -> Result<()> {
        // Every attribute lookup is a linear borrowed search, bounded before use.
        m.spend(a.len().saturating_mul(16).saturating_add(1), 0)?;
        match tag {
            "CvReference" => {
                cap(self.references.len().saturating_add(1), o.max_references)?;
                let r = CVReference {
                    name: m.copy(required(a, "cvName")?)?,
                    identifier: m.copy(required(a, "cvIdentifier")?)?,
                };
                m.push(&mut self.references, r)?;
            }
            "CvMappingRule" => {
                self.rule.identifier = m.copy(required(a, "id")?)?;
                let path = required(a, "cvElementPath")?;
                self.rule.element_path = if o.strip_namespaces {
                    strip_path(path, m)?
                } else {
                    m.copy(path)?
                };
                self.rule.requirement_level = match required(a, "requirementLevel")? {
                    "MAY" => RequirementLevel::May,
                    "SHOULD" => RequirementLevel::Should,
                    _ => RequirementLevel::Must,
                };
                self.rule.scope_path = m.copy(required(a, "scopePath")?)?;
                self.rule.combinations_logic = match required(a, "cvTermsCombinationLogic")? {
                    "AND" => CombinationsLogic::And,
                    "XOR" => CombinationsLogic::Xor,
                    _ => CombinationsLogic::Or,
                };
            }
            "CvTerm" => {
                self.terms = self.terms.checked_add(1).ok_or_else(limit)?;
                cap(self.terms, o.max_terms)?;
                let term = CVMappingTerm {
                    accession: m.copy(required(a, "termAccession")?)?,
                    use_term: boolean(required(a, "useTerm")?)?,
                    use_term_name: match get(a, "useTermName") {
                        None | Some("") => false,
                        Some(v) => boolean(v)?,
                    },
                    term_name: m.copy(required(a, "termName")?)?,
                    is_repeatable: match get(a, "isRepeatable") {
                        None | Some("") => true,
                        Some(v) => boolean(v)?,
                    },
                    allow_children: boolean(required(a, "allowChildren")?)?,
                    cv_identifier_ref: m.copy(required(a, "cvIdentifierRef")?)?,
                };
                m.push(&mut self.rule.terms, term)?;
            }
            _ => {}
        }
        Ok(())
    }
    fn end(&mut self, tag: &str, o: &ReadOptions, m: &mut Meter) -> Result<()> {
        if tag == "CvMappingRule" {
            cap(self.rules.len().saturating_add(1), o.max_rules)?;
            let r = std::mem::take(&mut self.rule);
            m.push(&mut self.rules, r)?;
        }
        Ok(())
    }
}
fn parse(
    text: &str,
    o: &ReadOptions,
    m: &mut Meter,
) -> Result<(Vec<CVReference>, Vec<CVMappingRule>)> {
    // The borrowed reader does not build an owned XML tree or copy event bodies.
    m.spend(text.len().saturating_mul(16), 0)?;
    let mut reader = Reader::from_str(text);
    reader.config_mut().check_end_names = false;
    reader.config_mut().check_comments = true;
    let mut stack: Vec<String> = Vec::new();
    let (mut roots, mut elements, mut seen_any) = (0usize, 0usize, false);
    let mut draft = Draft {
        references: Vec::new(),
        rules: Vec::new(),
        rule: CVMappingRule::default(),
        terms: 0,
    };
    loop {
        let event = reader.read_event().map_err(|e| bad(e.to_string()))?;
        let empty = matches!(&event, Event::Empty(_));
        match event {
            Event::Start(e) | Event::Empty(e) => {
                elements = elements.checked_add(1).ok_or_else(limit)?;
                cap(elements, o.max_elements)?;
                cap(stack.len().saturating_add(1), o.max_depth)?;
                if stack.is_empty() {
                    roots += 1;
                    if roots != 1 {
                        return Err(bad("multiple XML roots"));
                    }
                }
                let tag = std::str::from_utf8(e.name().0).map_err(|_| bad("invalid XML tag"))?;
                name(tag)?;
                let a = attrs(&e, m)?;
                draft.start(tag, &a, o, m)?;
                if empty {
                    draft.end(tag, o, m)?;
                } else {
                    let owned = m.copy(tag)?;
                    m.push(&mut stack, owned)?;
                }
            }
            Event::End(e) => {
                let tag = std::str::from_utf8(e.name().0).map_err(|_| bad("invalid XML tag"))?;
                name(tag)?;
                if stack.pop().as_deref() != Some(tag) {
                    return Err(bad("mismatched XML closing tag"));
                }
                draft.end(tag, o, m)?;
            }
            Event::Text(e) => {
                if e.as_ref().windows(3).any(|s| s == b"]]>") {
                    return Err(bad("']]>' outside CDATA"));
                }
                if stack.is_empty() && !e.as_ref().iter().copied().all(space) {
                    return Err(bad("text outside XML root"));
                }
            }
            Event::GeneralRef(e) => {
                if stack.is_empty() {
                    return Err(bad("entity outside XML root"));
                }
                m.spend(
                    e.as_ref().len().saturating_add(2),
                    e.as_ref().len().saturating_add(2).saturating_mul(2),
                )?;
                let raw = std::str::from_utf8(e.as_ref()).map_err(|_| bad("invalid XML entity"))?;
                let encoded = format!("&{raw};");
                let decoded =
                    quick_xml::escape::unescape(&encoded).map_err(|e| bad(e.to_string()))?;
                xml_text(&decoded)?;
            }
            Event::CData(_) => {
                if stack.is_empty() {
                    return Err(bad("CDATA outside XML root"));
                }
            }
            Event::Decl(_) => {
                if seen_any {
                    return Err(bad("misplaced XML declaration"));
                }
            }
            Event::PI(e) => {
                let target =
                    std::str::from_utf8(e.target()).map_err(|_| bad("invalid XML PI target"))?;
                name(target)?;
                if target.eq_ignore_ascii_case("xml") {
                    return Err(bad("reserved XML PI target"));
                }
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported(
                    "DTD and external entities in CV mapping XML".into(),
                ));
            }
            Event::Comment(_) => {}
            Event::Eof => {
                if roots != 1 || !stack.is_empty() {
                    return Err(bad("incomplete XML document"));
                }
                break;
            }
        }
        seen_any = true;
    }
    Ok((draft.references, draft.rules))
}

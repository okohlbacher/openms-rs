// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Shared bounded XML input for CV mapping and general semantic validation.

use crate::{Error, Result};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::{io::BufRead, mem::size_of};

pub(super) struct XmlLimits {
    pub max_input_bytes: usize,
    pub max_depth: usize,
    pub max_elements: usize,
}
fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
pub(super) struct Meter {
    pub work: usize,
    pub bytes: usize,
    context: &'static str,
}
impl Meter {
    pub fn new(work: usize, bytes: usize, context: &'static str) -> Self {
        Self {
            work,
            bytes,
            context,
        }
    }
    pub fn limit(&self) -> Error {
        bad(format!("{} resource limit exceeded", self.context))
    }
    pub fn cap(&self, value: usize, maximum: usize) -> Result<()> {
        if value > maximum {
            Err(self.limit())
        } else {
            Ok(())
        }
    }
    pub fn spend(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.work = self.work.checked_sub(work).ok_or_else(|| self.limit())?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(|| self.limit())?;
        Ok(())
    }
    pub fn slots<T>(&mut self, n: usize) -> Result<()> {
        self.spend(
            n,
            n.checked_mul(size_of::<T>()).ok_or_else(|| self.limit())?,
        )
    }
    pub fn tree<T>(&mut self, n: usize) -> Result<()> {
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
                .ok_or_else(|| self.limit())?,
        )
    }
    pub fn copy(&mut self, s: &str) -> Result<String> {
        self.spend(s.len().saturating_add(1), s.len())?;
        Ok(s.into())
    }
    pub fn push<T>(&mut self, v: &mut Vec<T>, item: T) -> Result<()> {
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
pub(super) type Attributes = Vec<(String, String)>;
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
// Read complete bounded input before borrowed XML events. A failed short read
// propagates directly; no uninitialized compression sniff or parser reuse state.
pub(super) fn document(mut input: impl BufRead, o: &XmlLimits, m: &mut Meter) -> Result<String> {
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
            return Err(m.limit());
        }
        let wanted = bytes.len().checked_add(n).ok_or_else(|| m.limit())?;
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
            m.cap(result.len().saturating_add(c.len_utf8()), o.max_input_bytes)?;
            result.push(c);
        }
        result
    } else {
        let content = bytes.strip_prefix(&[239, 187, 191]).unwrap_or(&bytes);
        std::str::from_utf8(content)
            .map_err(|_| {
                Error::Unsupported(format!(
                    "{} input requires UTF-8, UTF-16 or ASCII-compatible bytes",
                    m.context
                ))
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
                    "{} XML encoding {encoding}",
                    m.context
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
    m.cap(text.len(), o.max_input_bytes)?;
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

pub(super) enum Element<'a> {
    Start(&'a str, &'a Attributes),
    End(&'a str),
}
// Ancestors exclude the current element, including for empty-element End calls.
pub(super) fn scan_elements(
    text: &str,
    o: &XmlLimits,
    m: &mut Meter,
    mut callback: impl FnMut(Element<'_>, &[String], &mut Meter) -> Result<()>,
) -> Result<()> {
    // The borrowed reader does not build an owned XML tree or copy event bodies.
    m.spend(text.len().saturating_mul(16), 0)?;
    let mut reader = Reader::from_str(text);
    reader.config_mut().check_end_names = false;
    reader.config_mut().check_comments = true;
    let mut stack: Vec<String> = Vec::new();
    let (mut roots, mut elements, mut seen_any) = (0usize, 0usize, false);
    loop {
        let event = reader.read_event().map_err(|e| bad(e.to_string()))?;
        let empty = matches!(&event, Event::Empty(_));
        match event {
            Event::Start(e) | Event::Empty(e) => {
                elements = elements.checked_add(1).ok_or_else(|| m.limit())?;
                m.cap(elements, o.max_elements)?;
                m.cap(stack.len().saturating_add(1), o.max_depth)?;
                if stack.is_empty() {
                    roots += 1;
                    if roots != 1 {
                        return Err(bad("multiple XML roots"));
                    }
                }
                let tag = std::str::from_utf8(e.name().0).map_err(|_| bad("invalid XML tag"))?;
                name(tag)?;
                let a = attrs(&e, m)?;
                callback(Element::Start(tag, &a), &stack, m)?;
                if empty {
                    callback(Element::End(tag), &stack, m)?;
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
                callback(Element::End(tag), &stack, m)?;
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
                return Err(Error::Unsupported(format!(
                    "DTD and external entities in {} XML",
                    m.context
                )));
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
    Ok(())
}

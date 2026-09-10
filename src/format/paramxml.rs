// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenMS parameter XML (INI) transport, through version 1.8.0.
//! Reads the source scalar/list types, legacy restrictions and file tags.
//! See `docs/PARAMXML_SUPPORT.md` for deliberate checked adaptations.

use crate::data_structures::list::ListParse;
use crate::param::{Param, ParamBuilder, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub const VERSION: &str = "1.8.0";

/// Allocation/work bounds apply before parsing or writing the corresponding data.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub max_xml_bytes: usize,
    pub max_elements: usize,
    pub max_depth: usize,
    pub max_list_items: usize,
    pub max_path_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_xml_bytes: 64 * 1024 * 1024,
            max_elements: 1_000_000,
            max_depth: 128,
            max_list_items: 1_000_000,
            max_path_bytes: 65_536,
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

/// Read an owned parameter tree. No external entities or schemas are fetched.
pub fn read(input: impl Read) -> Result<Param> {
    read_with_limits(input, Limits::default())
}
pub fn read_with_limits(input: impl Read, limits: Limits) -> Result<Param> {
    parse(input, Param::new(), limits)
}
/// Source load semantics update existing keys and retain keys absent from the file.
/// A late parse or I/O failure leaves the caller's tree unchanged.
pub fn read_into(input: impl Read, target: &mut Param, limits: Limits) -> Result<()> {
    let draft = parse(input, target.checked_clone()?, limits)?;
    *target = draft;
    Ok(())
}
pub fn load(path: impl AsRef<Path>) -> Result<Param> {
    read(File::open(path)?)
}
pub fn load_into(path: impl AsRef<Path>, target: &mut Param) -> Result<()> {
    read_into(File::open(path)?, target, Limits::default())
}

fn document(input: impl Read, limit: usize) -> Result<String> {
    let count = u64::try_from(limit)
        .ok()
        .and_then(|v| v.checked_add(1))
        .ok_or_else(|| bad("parameter XML byte limit overflows"))?;
    let mut bytes = Vec::new();
    input.take(count).read_to_end(&mut bytes)?;
    if bytes.len() > limit {
        return Err(bad("parameter XML byte limit exceeded"));
    }
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
            let c = value.map_err(|_| bad("invalid UTF-16 parameter XML"))?;
            if decoded.len().saturating_add(c.len_utf8()) > limit {
                return Err(bad("decoded parameter XML byte limit exceeded"));
            }
            decoded.push(c);
        }
        decoded
    } else {
        let bom = bytes.starts_with(&[0xef, 0xbb, 0xbf]);
        let bytes = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(&bytes);
        let encoding = encoding(bytes)?;
        if bom
            && encoding
                .as_deref()
                .is_some_and(|s| !matches!(s, "utf-8" | "utf8"))
        {
            return Err(bad("XML declaration conflicts with UTF-8 byte order mark"));
        }
        match encoding.as_deref().unwrap_or("utf-8") {
            "utf-8" | "utf8" => std::str::from_utf8(bytes)
                .map_err(|_| bad("invalid UTF-8 parameter XML"))?
                .to_owned(),
            "us-ascii" | "ascii" if bytes.is_ascii() => String::from_utf8(bytes.to_vec()).unwrap(),
            "iso-8859-1" | "iso8859-1" | "latin1" => {
                let mut decoded = String::new();
                for &byte in bytes {
                    let c = char::from(byte);
                    if decoded.len().saturating_add(c.len_utf8()) > limit {
                        return Err(bad("decoded parameter XML byte limit exceeded"));
                    }
                    decoded.push(c);
                }
                decoded
            }
            value => return Err(unsupported(format!("parameter XML encoding {value}"))),
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
    xml_text(&text)?;
    // XML 1.0 line ending normalization happens before attribute normalization.
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
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

fn xml_text(text: &str) -> Result<()> {
    if text.chars().all(|c| {
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

fn valid_pi_target(bytes: &[u8]) -> bool {
    let Ok(name) = std::str::from_utf8(bytes) else {
        return false;
    };
    if name.eq_ignore_ascii_case("xml") {
        return false;
    }
    let mut chars = name.chars();
    chars.next().is_some_and(xml_name_start) && chars.all(|c| {
        xml_name_start(c)
            || c.is_ascii_digit()
            || matches!(c, '-' | '.' | '\u{b7}' | '\u{300}'..='\u{36f}' | '\u{203f}'..='\u{2040}')
    })
}
fn xml_name_start(c: char) -> bool {
    c.is_ascii_alphabetic()
        || matches!(c, ':' | '_'
        | '\u{c0}'..='\u{d6}' | '\u{d8}'..='\u{f6}' | '\u{f8}'..='\u{2ff}'
        | '\u{370}'..='\u{37d}' | '\u{37f}'..='\u{1fff}' | '\u{200c}'..='\u{200d}'
        | '\u{2070}'..='\u{218f}' | '\u{2c00}'..='\u{2fef}' | '\u{3001}'..='\u{d7ff}'
        | '\u{f900}'..='\u{fdcf}' | '\u{fdf0}'..='\u{fffd}' | '\u{10000}'..='\u{effff}')
}

type Attributes = BTreeMap<String, String>;
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
fn attributes(element: &BytesStart<'_>) -> Result<Attributes> {
    attribute_spacing(element)?;
    let allowed: &[&str] = match element.name().as_ref() {
        b"PARAMETERS" => &["version", "xmlns:xsi", "xsi:noNamespaceSchemaLocation"],
        b"NODE" => &["name", "description"],
        b"ITEM" | b"ITEMLIST" => &[
            "name",
            "type",
            "value",
            "description",
            "tags",
            "required",
            "advanced",
            "restrictions",
            "supported_formats",
            "short_description",
            "position",
            "default",
        ],
        b"LISTITEM" => &["value"],
        _ => return Err(unsupported("unknown parameter XML element")),
    };
    let mut result = Attributes::new();
    for item in element.attributes() {
        let item = item.map_err(|e| bad(e.to_string()))?;
        let key = std::str::from_utf8(item.key.as_ref())
            .map_err(|_| bad("invalid XML attribute name"))?;
        if !allowed.contains(&key) {
            return Err(unsupported(format!("parameter XML attribute {key}")));
        }
        let raw = std::str::from_utf8(&item.value).map_err(|_| bad("invalid XML attribute"))?;
        if raw.contains('<') {
            return Err(bad("unescaped less-than sign in XML attribute"));
        }
        let normalized = raw.replace(['\t', '\n', '\r'], " ");
        let value = quick_xml::escape::unescape(&normalized)
            .map_err(|e| bad(e.to_string()))?
            .into_owned();
        xml_text(&value)?;
        result.insert(key.into(), value);
    }
    Ok(result)
}
fn required<'a>(attrs: &'a Attributes, key: &str) -> Result<&'a str> {
    attrs
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| bad(format!("missing XML attribute {key}")))
}
fn optional<'a>(attrs: &'a Attributes, key: &str) -> &'a str {
    attrs.get(key).map_or("", String::as_str)
}
fn check_attrs(attrs: &Attributes, allowed: &[&str]) -> Result<()> {
    if let Some(key) = attrs.keys().find(|key| !allowed.contains(&key.as_str())) {
        Err(unsupported(format!("parameter XML attribute {key}")))
    } else {
        Ok(())
    }
}
fn split(value: &str, limits: Limits) -> Result<Vec<String>> {
    if value.is_empty() {
        Ok(Vec::new())
    } else {
        let count = value
            .bytes()
            .filter(|b| *b == b',')
            .count()
            .saturating_add(1);
        let bytes = count
            .checked_mul(std::mem::size_of::<String>())
            .and_then(|slots| slots.checked_add(value.len()))
            .ok_or_else(|| bad("parameter attribute list size overflows"))?;
        if count > limits.max_list_items
            || bytes > limits.max_xml_bytes.min(crate::param::MAX_PARAM_BYTES)
        {
            return Err(bad("parameter attribute list allocation limit exceeded"));
        }
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(count)
            .map_err(|_| bad("parameter attribute list allocation failed"))?;
        fields.extend(value.split(',').map(str::to_owned));
        Ok(fields)
    }
}
fn int(value: &str) -> Result<i32> {
    i32::from_list_item(value).map_err(|_| bad(format!("invalid parameter integer {value:?}")))
}
fn float(value: &str) -> Result<f64> {
    f64::from_list_item(value).map_err(|_| bad(format!("invalid parameter float {value:?}")))
}

struct Pending {
    key: String,
    attrs: Attributes,
    value: ParamValue,
    limits: Limits,
}
impl Pending {
    fn new(path: &str, attrs: Attributes, list: bool, limits: Limits) -> Result<Self> {
        if list && (attrs.contains_key("value") || attrs.contains_key("default")) {
            return Err(bad(
                "ITEMLIST cannot contain scalar value/default attributes",
            ));
        }
        check_attrs(
            &attrs,
            &[
                "name",
                "type",
                "value",
                "description",
                "tags",
                "required",
                "advanced",
                "restrictions",
                "supported_formats",
                "short_description",
                "position",
                "default",
            ],
        )?;
        let name = required(&attrs, "name")?;
        if name.is_empty() {
            return Err(bad("empty parameter name"));
        }
        if path.len().saturating_add(name.len()) > limits.max_path_bytes {
            return Err(bad("parameter path limit exceeded"));
        }
        let kind = required(&attrs, "type")?;
        let raw = if list { "" } else { required(&attrs, "value")? };
        let value = match (kind, list) {
            ("int", false) => ParamValue::Integer(i64::from(int(raw)?)),
            ("float" | "double", false) => ParamValue::Float(float(raw)?),
            ("string" | "bool" | "input-file" | "output-file" | "output-prefix", false) => {
                ParamValue::String(raw.into())
            }
            ("int", true) => ParamValue::IntegerList(Vec::new()),
            ("float" | "double", true) => ParamValue::FloatList(Vec::new()),
            ("string" | "input-file" | "output-file", true) => ParamValue::StringList(Vec::new()),
            _ => {
                return Err(unsupported(format!(
                    "parameter XML {} type {kind}",
                    if list { "list" } else { "scalar" }
                )));
            }
        };
        Ok(Self {
            key: format!("{path}{name}"),
            attrs,
            value,
            limits,
        })
    }
    fn push(&mut self, attrs: Attributes, limits: Limits) -> Result<()> {
        check_attrs(&attrs, &["value"])?;
        let raw = required(&attrs, "value")?;
        match &mut self.value {
            ParamValue::StringList(values) => {
                check_list(values.len(), limits)?;
                values.push(raw.into());
            }
            ParamValue::IntegerList(values) => {
                check_list(values.len(), limits)?;
                values.push(int(raw)?);
            }
            ParamValue::FloatList(values) => {
                check_list(values.len(), limits)?;
                values.push(float(raw)?);
            }
            _ => return Err(bad("LISTITEM outside ITEMLIST")),
        }
        Ok(())
    }
    fn store(self, param: &mut ParamBuilder) -> Result<()> {
        let kind = required(&self.attrs, "type")?;
        let mut tags = split(optional(&self.attrs, "tags"), self.limits)?;
        for tag in ["required", "advanced"] {
            if optional(&self.attrs, tag) == "true" {
                tags.push(tag.into());
            }
        }
        match kind {
            "input-file" => tags.push("input file".into()),
            "output-file" => tags.push("output file".into()),
            "output-prefix" => tags.push("output prefix".into()),
            _ => (),
        }
        let files = tags
            .iter()
            .any(|s| matches!(s.as_str(), "input file" | "output file" | "output prefix"));
        let key = &self.key;
        param.set_value(
            key,
            self.value,
            &optional(&self.attrs, "description").replace("#br#", "\n"),
            &tags,
        )?;
        if kind == "bool" {
            param.set_valid_strings(key, &["true".into(), "false".into()])?;
        } else {
            let restrictions = if files {
                self.attrs
                    .get("supported_formats")
                    .or_else(|| self.attrs.get("restrictions"))
            } else {
                self.attrs.get("restrictions")
            };
            if let Some(raw) = restrictions {
                match kind {
                    "string" | "input-file" | "output-file" | "output-prefix" => {
                        param.set_valid_strings(key, &split(raw, self.limits)?)?
                    }
                    "int" | "float" | "double" => {
                        let mut parts: Vec<_> = raw.split(':').collect();
                        if parts.len() != 2 {
                            parts = raw.split('-').collect();
                        }
                        // Source ignores malformed range shape (legacy empty attributes included).
                        if parts.len() == 2 {
                            if kind == "int" {
                                if !parts[0].is_empty() {
                                    param.set_min_int(key, int(parts[0])?)?;
                                }
                                if !parts[1].is_empty() {
                                    param.set_max_int(key, int(parts[1])?)?;
                                }
                            } else {
                                if !parts[0].is_empty() {
                                    param.set_min_float(key, float(parts[0])?)?;
                                }
                                if !parts[1].is_empty() {
                                    param.set_max_float(key, float(parts[1])?)?;
                                }
                            }
                        }
                    }
                    _ => (),
                }
            }
        }
        Ok(())
    }
}
fn check_list(count: usize, limits: Limits) -> Result<()> {
    if count >= limits.max_list_items {
        Err(bad("parameter list item limit exceeded"))
    } else {
        Ok(())
    }
}

fn parse(input: impl Read, param: Param, limits: Limits) -> Result<Param> {
    let mut param = ParamBuilder::new(param)?;
    let text = document(input, limits.max_xml_bytes)?;
    let mut reader = Reader::from_str(&text);
    reader.config_mut().check_comments = true;
    let mut stack: Vec<String> = Vec::new();
    let mut paths: Vec<usize> = Vec::new();
    let mut path = String::new();
    let mut pending: Option<Pending> = None;
    let mut root_seen = false;
    let mut declaration_allowed = true;
    let mut records = 0usize;
    loop {
        let event = reader.read_event().map_err(|e| bad(e.to_string()))?;
        match event {
            Event::Start(ref element) | Event::Empty(ref element) => {
                let empty = matches!(event, Event::Empty(_));
                if records >= limits.max_elements {
                    return Err(bad("parameter XML element limit exceeded"));
                }
                records += 1;
                if stack.len() >= limits.max_depth {
                    return Err(bad("parameter XML nesting limit exceeded"));
                }
                let name = std::str::from_utf8(element.name().as_ref())
                    .map_err(|_| bad("invalid XML name"))?
                    .to_owned();
                let attrs = attributes(element)?;
                let parent = stack.last().map(String::as_str);
                declaration_allowed = false;
                match name.as_str() {
                    "PARAMETERS" if parent.is_none() && !root_seen => {
                        root_seen = true;
                        check_attrs(
                            &attrs,
                            &["version", "xmlns:xsi", "xsi:noNamespaceSchemaLocation"],
                        )?;
                        if attrs
                            .get("xmlns:xsi")
                            .is_some_and(|v| v != "http://www.w3.org/2001/XMLSchema-instance")
                            || (attrs.contains_key("xsi:noNamespaceSchemaLocation")
                                && !attrs.contains_key("xmlns:xsi"))
                        {
                            return Err(bad(
                                "invalid or unbound parameter XML schema-instance namespace",
                            ));
                        }
                        let version = optional(&attrs, "version");
                        if !version.is_empty() {
                            let parts = version
                                .split('.')
                                .map(str::parse::<u32>)
                                .collect::<std::result::Result<Vec<_>, _>>()
                                .map_err(|_| bad("invalid parameter XML version"))?;
                            if !(2..=3).contains(&parts.len()) {
                                return Err(bad("invalid parameter XML version"));
                            }
                            if (parts[0], parts[1], *parts.get(2).unwrap_or(&0)) > (1, 8, 0) {
                                return Err(unsupported(format!(
                                    "parameter XML version {version}"
                                )));
                            }
                        }
                    }
                    "NODE" if matches!(parent, Some("PARAMETERS" | "NODE")) => {
                        check_attrs(&attrs, &["name", "description"])?;
                        let node = required(&attrs, "name")?;
                        if node.is_empty() {
                            return Err(bad("empty parameter section name"));
                        }
                        if path.len().saturating_add(node.len()).saturating_add(1)
                            > limits.max_path_bytes
                        {
                            return Err(bad("parameter path limit exceeded"));
                        }
                        paths.push(path.len());
                        path.push_str(node);
                        param.add_section(
                            &path,
                            &optional(&attrs, "description").replace("#br#", "\n"),
                        )?;
                        path.push(':');
                    }
                    "ITEM" | "ITEMLIST" if matches!(parent, Some("PARAMETERS" | "NODE")) => {
                        pending = Some(Pending::new(&path, attrs, name == "ITEMLIST", limits)?);
                    }
                    "LISTITEM" if parent == Some("ITEMLIST") => {
                        pending
                            .as_mut()
                            .ok_or_else(|| bad("LISTITEM has no list"))?
                            .push(attrs, limits)?;
                    }
                    _ => {
                        return Err(unsupported(format!(
                            "parameter XML element {name} in {parent:?}"
                        )));
                    }
                }
                stack.push(name);
                if empty {
                    close(&mut stack, &mut paths, &mut path, &mut pending, &mut param)?;
                }
            }
            Event::End(_) => close(&mut stack, &mut paths, &mut path, &mut pending, &mut param)?,
            Event::Text(value) if value.as_ref().iter().all(u8::is_ascii_whitespace) => {
                declaration_allowed = false;
            }
            Event::Comment(_) => {
                declaration_allowed = false;
            }
            Event::PI(value) if valid_pi_target(value.target()) => {
                declaration_allowed = false;
            }
            Event::Decl(_) if declaration_allowed => {
                declaration_allowed = false;
            }
            Event::Eof if root_seen && stack.is_empty() => return param.finish(),
            Event::Eof => return Err(bad("missing or unclosed parameter XML root")),
            Event::DocType(_) => return Err(unsupported("parameter XML DTD/entity declarations")),
            _ => return Err(bad("unexpected parameter XML content")),
        }
    }
}
fn close(
    stack: &mut Vec<String>,
    paths: &mut Vec<usize>,
    path: &mut String,
    pending: &mut Option<Pending>,
    param: &mut ParamBuilder,
) -> Result<()> {
    match stack.pop().as_deref() {
        Some("NODE") => path.truncate(paths.pop().ok_or_else(|| bad("unmatched parameter node"))?),
        Some("ITEM" | "ITEMLIST") => pending
            .take()
            .ok_or_else(|| bad("missing parameter entry"))?
            .store(param)?,
        Some(_) => (),
        None => return Err(bad("unmatched parameter XML end tag")),
    }
    Ok(())
}

/// Validate and serialize fully before touching the destination stream; flush errors propagate.
pub fn write(output: impl Write, param: &Param) -> Result<()> {
    write_with_limits(output, param, Limits::default())
}
pub fn write_with_limits(mut output: impl Write, param: &Param, limits: Limits) -> Result<()> {
    let mut xml = Output {
        bytes: Vec::new(),
        limits,
        elements: 0,
    };
    xml.add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<PARAMETERS version=\"1.8.0\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/Param_1_8_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n")?;
    xml.element(0)?;
    xml.node(param.root(), 1, 0)?;
    xml.add("</PARAMETERS>\n")?;
    output.write_all(&xml.bytes)?;
    output.flush()?;
    Ok(())
}
/// A filename of `-` writes to stdout, matching the source API.
pub fn store(path: impl AsRef<Path>, param: &Param) -> Result<()> {
    if path.as_ref() == Path::new("-") {
        return write(std::io::stdout().lock(), param);
    }
    let mut bytes = Vec::new();
    write(&mut bytes, param)?;
    let mut file = File::create(path)?;
    file.write_all(&bytes)?;
    file.flush()?;
    Ok(())
}

struct Output {
    bytes: Vec<u8>,
    limits: Limits,
    elements: usize,
}
impl Output {
    fn add(&mut self, text: &str) -> Result<()> {
        if self.bytes.len().saturating_add(text.len()) > self.limits.max_xml_bytes {
            return Err(bad("parameter XML output byte limit exceeded"));
        }
        self.bytes.extend_from_slice(text.as_bytes());
        Ok(())
    }
    fn element(&mut self, depth: usize) -> Result<()> {
        if self.elements >= self.limits.max_elements || depth >= self.limits.max_depth {
            return Err(bad("parameter XML output element/depth limit exceeded"));
        }
        self.elements += 1;
        Ok(())
    }
    fn attr(&mut self, name: &str, value: &str) -> Result<()> {
        self.add(" ")?;
        self.add(name)?;
        self.add("=\"")?;
        self.escaped(value, false)?;
        self.add("\"")
    }
    fn escaped(&mut self, value: &str, description: bool) -> Result<()> {
        xml_text(value)?;
        // Escape incrementally so a small byte limit cannot cause an oversized temporary.
        for c in value.chars() {
            match c {
                '&' => self.add("&amp;")?,
                '<' => self.add("&lt;")?,
                '>' => self.add("&gt;")?,
                '"' => self.add("&quot;")?,
                '\'' => self.add("&apos;")?,
                '\t' => self.add("&#x9;")?,
                '\n' if description => self.add("#br#")?,
                '\n' => self.add("&#xA;")?,
                '\r' => self.add("&#xD;")?,
                _ => {
                    let mut buffer = [0; 4];
                    self.add(c.encode_utf8(&mut buffer))?;
                }
            }
        }
        Ok(())
    }
    fn indent(&mut self, depth: usize) -> Result<()> {
        for _ in 0..depth {
            self.add("  ")?;
        }
        Ok(())
    }
    fn node(&mut self, node: &ParamNode, depth: usize, path_bytes: usize) -> Result<()> {
        if node.entries.len().saturating_add(node.nodes.len())
            > self.limits.max_elements.saturating_sub(self.elements)
        {
            return Err(bad("parameter XML output element limit exceeded"));
        }
        let mut seen = BTreeSet::new();
        for key in node
            .entries
            .iter()
            .map(|e| (false, e.name.as_str()))
            .chain(node.nodes.iter().map(|n| (true, n.name.as_str())))
        {
            if key.1.is_empty() || key.1.contains(':') || !seen.insert(key) {
                return Err(unsupported(
                    "parameter XML requires unique nonempty local names without colons",
                ));
            }
        }
        drop(seen);
        for entry in &node.entries {
            if path_bytes.saturating_add(entry.name.len()) > self.limits.max_path_bytes {
                return Err(bad("parameter path limit exceeded"));
            }
            self.entry(entry, depth)?;
        }
        for child in &node.nodes {
            let bytes = path_bytes
                .saturating_add(child.name.len())
                .saturating_add(1);
            if bytes > self.limits.max_path_bytes {
                return Err(bad("parameter path limit exceeded"));
            }
            self.element(depth)?;
            self.indent(depth)?;
            self.add("<NODE")?;
            self.attr("name", &child.name)?;
            self.description(&child.description)?;
            self.add(">\n")?;
            self.node(child, depth + 1, bytes)?;
            self.indent(depth)?;
            self.add("</NODE>\n")?;
        }
        Ok(())
    }
    fn description(&mut self, value: &str) -> Result<()> {
        if value.contains("#br#") {
            return Err(unsupported(
                "literal #br# in parameter description cannot round-trip through the source reader",
            ));
        }
        if value.len() > self.limits.max_xml_bytes {
            return Err(bad("parameter description byte limit exceeded"));
        }
        self.add(" description=\"")?;
        self.escaped(value, true)?;
        self.add("\"")
    }
    fn entry(&mut self, entry: &ParamEntry, depth: usize) -> Result<()> {
        if matches!(entry.value, ParamValue::Empty) {
            return Err(unsupported(
                "empty parameter value has no XML representation",
            ));
        }
        self.element(depth)?;
        let (list, mut kind) = match entry.value {
            ParamValue::Integer(_) => (false, "int"),
            ParamValue::Float(_) => (false, "double"),
            ParamValue::String(_) => (false, "string"),
            ParamValue::StringList(_) => (true, "string"),
            ParamValue::IntegerList(_) => (true, "int"),
            ParamValue::FloatList(_) => (true, "double"),
            ParamValue::Empty => unreachable!(),
        };
        let ints = matches!(
            entry.value,
            ParamValue::Integer(_) | ParamValue::IntegerList(_)
        );
        let floats = matches!(entry.value, ParamValue::Float(_) | ParamValue::FloatList(_));
        if (!ints && (entry.min_int != -i32::MAX || entry.max_int != i32::MAX))
            || (!floats && (entry.min_float != -f64::MAX || entry.max_float != f64::MAX))
            || ((ints || floats) && !entry.valid_strings.is_empty())
        {
            return Err(unsupported(
                "inactive parameter restrictions have no XML representation",
            ));
        }
        let mut type_tag = None;
        if kind == "string" {
            for (tag, xml_type) in [
                ("input file", "input-file"),
                ("output file", "output-file"),
                ("output prefix", "output-prefix"),
            ] {
                if tag == "output prefix" && list {
                    continue;
                }
                if entry.tags.contains(tag) {
                    kind = xml_type;
                    type_tag = Some(tag);
                    break;
                }
            }
            if !list
                && kind == "string"
                && entry.valid_strings == ["true", "false"]
                && matches!(&entry.value, ParamValue::String(v) if v == "false")
            {
                kind = "bool";
            }
        }
        self.indent(depth)?;
        self.add(if list { "<ITEMLIST" } else { "<ITEM" })?;
        self.attr("name", &entry.name)?;
        if !list {
            match &entry.value {
                ParamValue::String(v) => self.attr("value", v)?,
                ParamValue::Integer(v) => {
                    i32::try_from(*v).map_err(|_| {
                        unsupported("parameter XML source reader uses 32-bit integers")
                    })?;
                    self.attr("value", &v.to_string())?;
                }
                ParamValue::Float(v) => self.attr("value", &float_text(*v))?,
                _ => unreachable!(),
            }
        }
        self.attr("type", kind)?;
        self.description(&entry.description)?;
        for tag in ["required", "advanced"] {
            self.attr(
                tag,
                if entry.tags.contains(tag) {
                    "true"
                } else {
                    "false"
                },
            )?;
        }
        let remaining_tags = entry.tags.iter().filter(|s| {
            !matches!(s.as_str(), "required" | "advanced") && Some(s.as_str()) != type_tag
        });
        if remaining_tags.clone().next().is_some() {
            self.add(" tags=\"")?;
            for (index, tag) in remaining_tags.enumerate() {
                if tag.is_empty() || tag.contains(',') {
                    return Err(unsupported(
                        "empty or comma-containing parameter tag cannot round-trip",
                    ));
                }
                if index != 0 {
                    self.add(",")?;
                }
                self.escaped(tag, false)?;
            }
            self.add("\"")?;
        }
        if kind != "bool" {
            let restrictions = match &entry.value {
                ParamValue::Integer(_) | ParamValue::IntegerList(_) => {
                    bounds(entry.min_int, entry.max_int, -i32::MAX, i32::MAX, |v| {
                        v.to_string()
                    })
                }
                ParamValue::Float(_) | ParamValue::FloatList(_) => bounds(
                    entry.min_float,
                    entry.max_float,
                    -f64::MAX,
                    f64::MAX,
                    float_text,
                ),
                _ => {
                    if entry.valid_strings.iter().any(|v| v.contains(',')) {
                        return Err(unsupported(
                            "comma-containing valid string cannot round-trip",
                        ));
                    }
                    let size = entry
                        .valid_strings
                        .iter()
                        .try_fold(0usize, |n, s| n.checked_add(s.len().saturating_add(1)))
                        .ok_or_else(|| bad("restriction bytes overflow"))?;
                    if size > self.limits.max_xml_bytes {
                        return Err(bad("restriction byte limit exceeded"));
                    }
                    if entry.valid_strings == [""] {
                        return Err(unsupported("singleton empty restriction cannot round-trip"));
                    }
                    entry.valid_strings.join(",")
                }
            };
            if !restrictions.is_empty() {
                let file = entry
                    .tags
                    .iter()
                    .any(|s| matches!(s.as_str(), "input file" | "output file" | "output prefix"));
                self.attr(
                    if file {
                        "supported_formats"
                    } else {
                        "restrictions"
                    },
                    &restrictions,
                )?;
            }
        }
        if !list {
            return self.add(" />\n");
        }
        self.add(">\n")?;
        match &entry.value {
            ParamValue::StringList(values) => {
                if values.len() > self.limits.max_list_items {
                    return Err(bad("parameter list output limit exceeded"));
                }
                for value in values {
                    self.list_item(value, depth + 1)?;
                }
            }
            ParamValue::IntegerList(values) => {
                if values.len() > self.limits.max_list_items {
                    return Err(bad("parameter list output limit exceeded"));
                }
                for value in values {
                    self.list_item(&value.to_string(), depth + 1)?;
                }
            }
            ParamValue::FloatList(values) => {
                if values.len() > self.limits.max_list_items {
                    return Err(bad("parameter list output limit exceeded"));
                }
                for value in values {
                    self.list_item(&float_text(*value), depth + 1)?;
                }
            }
            _ => unreachable!(),
        }
        self.indent(depth)?;
        self.add("</ITEMLIST>\n")
    }
    fn list_item(&mut self, value: &str, depth: usize) -> Result<()> {
        self.element(depth)?;
        self.indent(depth)?;
        self.add("<LISTITEM")?;
        self.attr("value", value)?;
        self.add("/>\n")
    }
}
fn bounds<T: PartialEq>(
    min: T,
    max: T,
    default_min: T,
    default_max: T,
    format: impl Fn(T) -> String,
) -> String {
    if min == default_min && max == default_max {
        return String::new();
    }
    let low = if min == default_min {
        String::new()
    } else {
        format(min)
    };
    let high = if max == default_max {
        String::new()
    } else {
        format(max)
    };
    format!("{low}:{high}")
}
fn float_text(value: f64) -> String {
    if value.is_nan() {
        "NaN".into()
    } else if value == f64::INFINITY {
        "INF".into()
    } else if value == f64::NEG_INFINITY {
        "-INF".into()
    } else {
        value.to_string()
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A bounded SAX-style event stream for the attribute-driven identification
//! readers that the source builds on `Internal::XMLHandler`: `ProtXMLFile`,
//! `PepXMLFileMascot` and `XQuestResultXMLHandler`.
//!
//! The source hands each handler Xerces' `startElement` / `endElement` /
//! `characters` callbacks with the *qualified* element and attribute names, as
//! written. This module reproduces that contract over `quick-xml`:
//!
//! - element and attribute names are passed as written, prefix included, so a
//!   handler matching `protein` does not match `p:protein`, exactly as the
//!   source's `sm_.convert(qname) == "protein"` does not;
//! - attribute values are entity-expanded and normalised as XML 1.0 §3.3.3
//!   requires (line ends first, then every tab, line feed and carriage return
//!   becomes a space), which is what Xerces delivers;
//! - character data is delivered with predefined and numeric references
//!   expanded; a reference that needs a DTD is refused, and so is a DTD;
//! - the document must be one well-formed element: one root, closed, no
//!   character data outside it, matching end tags, no duplicate attributes.
//!
//! Every ceiling is checked before the work it bounds: the decoded size (with
//! the shared identification-XML decoder, which also handles a byte-order mark,
//! UTF-16 and Latin-1), the number of events, the nesting depth and the number
//! of attributes on one element.

use crate::format::identification_xml::decode_document;
use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::io::BufRead;

/// Resource ceilings of one parse.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SaxLimits {
    /// Maximum decoded document size in bytes.
    pub(crate) max_bytes: usize,
    /// Maximum XML events, including discarded ones.
    pub(crate) max_events: usize,
    /// Maximum element nesting depth.
    pub(crate) max_depth: usize,
    /// Maximum attributes on one element.
    pub(crate) max_attributes: usize,
}

/// Attributes of one element, in document order, names as written.
pub(crate) type Attributes = Vec<(String, String)>;

/// One callback of the event stream.
pub(crate) enum SaxEvent<'a> {
    /// `startElement(qname, attributes)`.
    Start {
        name: &'a str,
        attributes: &'a Attributes,
    },
    /// `endElement(qname)`.
    End { name: &'a str },
    /// `characters(chars)`: one piece of character data inside the root, with
    /// references expanded. Like Xerces, one run of text may arrive in pieces.
    Text(&'a str),
}

/// The value of attribute `name`, `XMLAttributes::value(name)`.
pub(crate) fn attribute<'a>(attributes: &'a Attributes, name: &str) -> Option<&'a str> {
    attributes
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value.as_str())
}

fn malformed(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// Parse `input` and call `handler` for every start tag, end tag and piece of
/// character data, in document order. The first handler error stops the parse
/// and is returned.
pub(crate) fn parse(
    input: impl BufRead,
    limits: &SaxLimits,
    mut handler: impl FnMut(SaxEvent<'_>) -> Result<()>,
) -> Result<()> {
    if limits.max_events == 0 || limits.max_depth == 0 || limits.max_attributes == 0 {
        return Err(Error::InvalidValue("invalid XML limits".into()));
    }
    let text = decode_document(input, limits.max_bytes, None)?;
    let mut reader = Reader::from_str(&text);
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().check_end_names = true;
    let mut events = 0usize;
    let mut depth = 0usize;
    let mut root_opened = false;
    loop {
        let event = reader
            .read_event()
            .map_err(|error| malformed(format!("XML is not well formed: {error}")))?;
        events += 1;
        if events > limits.max_events {
            return Err(Error::InvalidValue("XML event limit exceeded".into()));
        }
        match event {
            Event::Start(start) => {
                if depth == 0 {
                    if root_opened {
                        return Err(malformed("XML holds more than one root element"));
                    }
                    root_opened = true;
                }
                if depth >= limits.max_depth {
                    return Err(Error::InvalidValue(
                        "XML nesting depth limit exceeded".into(),
                    ));
                }
                depth += 1;
                let name = std::str::from_utf8(start.name().into_inner())
                    .map_err(|_| malformed("XML element name is not UTF-8"))?
                    .to_owned();
                let attributes = attributes(&start, limits)?;
                handler(SaxEvent::Start {
                    name: &name,
                    attributes: &attributes,
                })?;
            }
            Event::End(end) => {
                let name = std::str::from_utf8(end.name().into_inner())
                    .map_err(|_| malformed("XML element name is not UTF-8"))?
                    .to_owned();
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| malformed(format!("closing tag {name} has no opening tag")))?;
                handler(SaxEvent::End { name: &name })?;
            }
            Event::Text(body) => {
                let decoded = body
                    .xml_content()
                    .map_err(|error| malformed(format!("XML text: {error}")))?;
                if depth == 0 {
                    if decoded.trim_matches([' ', '\t', '\n', '\r']).is_empty() {
                        continue;
                    }
                    return Err(malformed(
                        "XML holds character data outside the root element",
                    ));
                }
                handler(SaxEvent::Text(&decoded))?;
            }
            Event::GeneralRef(reference) => {
                if depth == 0 {
                    return Err(malformed("XML entity reference outside the root element"));
                }
                let expanded = expand_reference(reference.as_ref())?;
                handler(SaxEvent::Text(&expanded))?;
            }
            Event::CData(data) => {
                if depth == 0 {
                    return Err(malformed("CDATA outside the root element"));
                }
                let bytes: &[u8] = &data;
                let decoded =
                    std::str::from_utf8(bytes).map_err(|_| malformed("CDATA is not UTF-8"))?;
                handler(SaxEvent::Text(decoded))?;
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported("XML DTDs are not supported".into()));
            }
            Event::Empty(_) => {
                // `expand_empty_elements` turns every empty element into a
                // start and an end event.
                return Err(malformed("unexpected empty-element event"));
            }
            Event::Eof => break,
            Event::Decl(_) | Event::PI(_) | Event::Comment(_) => {}
        }
    }
    if !root_opened {
        return Err(malformed("XML document has no root element"));
    }
    if depth != 0 {
        return Err(malformed("XML document ends inside an element"));
    }
    Ok(())
}

/// Attributes as Xerces delivers them: names as written, values normalised
/// (XML 1.0 §2.11 line ends, then §3.3.3 white space) and entity-expanded.
fn attributes(element: &BytesStart<'_>, limits: &SaxLimits) -> Result<Attributes> {
    let mut out: Attributes = Vec::new();
    for attribute in element.attributes().with_checks(true) {
        let attribute = attribute.map_err(|error| malformed(format!("XML attribute: {error}")))?;
        if out.len() >= limits.max_attributes {
            return Err(Error::InvalidValue(
                "XML attribute count limit exceeded".into(),
            ));
        }
        let key = std::str::from_utf8(attribute.key.as_ref())
            .map_err(|_| malformed("XML attribute name is not UTF-8"))?;
        let raw = std::str::from_utf8(&attribute.value)
            .map_err(|_| malformed("XML attribute value is not UTF-8"))?;
        if raw.contains('<') {
            return Err(malformed("raw '<' in XML attribute value"));
        }
        let normalized: String = raw
            .replace("\r\n", "\n")
            .chars()
            .map(|c| {
                if matches!(c, '\t' | '\n' | '\r') {
                    ' '
                } else {
                    c
                }
            })
            .collect();
        let value = quick_xml::escape::unescape(&normalized)
            .map_err(|error| malformed(format!("XML attribute {key}: {error}")))?;
        out.push((key.to_owned(), value.into_owned()));
    }
    Ok(out)
}

/// A predefined or numeric character reference; anything else needs a DTD.
fn expand_reference(raw: &[u8]) -> Result<String> {
    let name = std::str::from_utf8(raw).map_err(|_| malformed("XML reference is not UTF-8"))?;
    match name {
        "amp" => return Ok("&".into()),
        "lt" => return Ok("<".into()),
        "gt" => return Ok(">".into()),
        "quot" => return Ok("\"".into()),
        "apos" => return Ok("'".into()),
        _ => {}
    }
    let digits = name.strip_prefix('#').ok_or_else(|| {
        Error::Unsupported(format!(
            "XML entity reference `&{name};` needs a DTD, which is not supported"
        ))
    })?;
    let code = match digits
        .strip_prefix('x')
        .or_else(|| digits.strip_prefix('X'))
    {
        Some(hex) => u32::from_str_radix(hex, 16),
        None => digits.parse::<u32>(),
    }
    .map_err(|_| malformed("malformed numeric character reference"))?;
    char::from_u32(code)
        .map(String::from)
        .ok_or_else(|| malformed("numeric character reference is not a Unicode scalar value"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIMITS: SaxLimits = SaxLimits {
        max_bytes: 1 << 20,
        max_events: 1000,
        max_depth: 8,
        max_attributes: 8,
    };

    fn events(text: &str) -> Result<Vec<String>> {
        let mut out = Vec::new();
        parse(text.as_bytes(), &LIMITS, |event| {
            out.push(match event {
                SaxEvent::Start { name, attributes } => format!("<{name} {attributes:?}"),
                SaxEvent::End { name } => format!(">{name}"),
                SaxEvent::Text(text) => format!("'{text}"),
            });
            Ok(())
        })?;
        Ok(out)
    }

    #[test]
    fn names_are_qualified_and_values_normalised() {
        let got = events("<a:r x=\"1&amp;\t2\r\n3\"><b/>t&#65;</a:r>").unwrap();
        assert_eq!(
            got,
            [
                "<a:r [(\"x\", \"1& 2 3\")]",
                "<b []",
                ">b",
                "'t",
                "'A",
                ">a:r"
            ]
        );
    }

    #[test]
    fn malformed_documents_are_refused() {
        for text in [
            "",
            "<a></b>",
            "<a/><b/>",
            "text<a/>",
            "<a x=\"1\" x=\"2\"/>",
            "<a>&custom;</a>",
            "<!DOCTYPE a><a/>",
            "<a>",
        ] {
            assert!(events(text).is_err(), "{text:?}");
        }
        let deep = "<a>".repeat(9) + &"</a>".repeat(9);
        assert!(matches!(events(&deep), Err(Error::InvalidValue(_))));
        let wide = format!(
            "<a {}/>",
            (0..9)
                .map(|i| format!("x{i}=\"\""))
                .collect::<Vec<_>>()
                .join(" ")
        );
        assert!(matches!(events(&wide), Err(Error::InvalidValue(_))));
    }
}

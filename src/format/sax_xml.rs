// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Streaming, bounded SAX-style events for the knowledge-base XML readers.
//!
//! The source parses PTMXML, Unimod XML and UniProtKB XML through
//! `Internal::XMLFile::parse_`, which drives a Xerces SAX2 parser with
//! namespace processing switched off (`XMLFile.cpp:82-83`), so every handler
//! sees *qualified* names (`umod:mod`, not `mod`). This module is the part of
//! that machinery those three handlers observe, as a pull loop over
//! `quick-xml` that hands each handler the same three callbacks the source
//! overrides: `onStartElement`, `onEndElement` and `onCharacters`.
//!
//! Preserved Xerces conventions:
//!
//! - Element and attribute names are delivered with their prefix, unresolved.
//! - Attribute values are normalised as XML 1.0 requires: line ends first
//!   (`\r\n` and a lone `\r` become `\n`), then each literal tab, line feed or
//!   carriage return becomes a space. Character references are expanded
//!   afterwards, so `&#10;` stays a line feed.
//! - Character data is delivered with its predefined and numeric references
//!   expanded and XML 1.0 line ends normalised. A run of character data
//!   between two pieces of markup is one callback, as Xerces delivers it; a
//!   CDATA section is its own callback; a comment or processing instruction
//!   ends a run.
//! - White space inside the root element is character data (there is no DTD
//!   to declare it ignorable), so handlers see it.
//!
//! Native differences, each an explicit refusal rather than a guess:
//!
//! - A document type declaration is [`Error::Unsupported`]: an internal
//!   subset can declare entities, which Xerces would expand and this reader
//!   cannot. None of the three formats uses one.
//! - A declared encoding other than UTF-8 (or US-ASCII, its subset) is
//!   [`Error::Unsupported`]; Xerces transcodes several more.
//! - Input bytes, nesting depth, events, attributes per element and the bytes
//!   of one character run are bounded by [`SaxLimits`]; exceeding a bound is
//!   [`Error::InvalidValue`], never a truncated document.
//!
//! Every well-formedness error Xerces reports as `Exception::ParseError`
//! (mismatched or unclosed tags, a second root, text outside the root, an
//! undeclared entity, a duplicate attribute) is [`Error::Parse`].

use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::borrow::Cow;
use std::io::{self, BufRead, Read};

/// Resource bounds for one parse. Each bound is checked before the work it
/// limits is done, so exceeding it never yields a partial result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SaxLimits {
    /// Decompressed input bytes the parser may consume.
    pub max_input_bytes: u64,
    /// Open elements at any one time.
    pub max_depth: usize,
    /// Parser events (start, end, text, reference and so on) in the document.
    pub max_events: u64,
    /// Attributes on a single element.
    pub max_attributes: usize,
    /// Bytes of one run of character data, after expansion.
    pub max_text_bytes: usize,
}

/// Attributes of one start tag, in document order, names unresolved.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SaxAttributes(Vec<(String, String)>);

impl SaxAttributes {
    /// The value of the attribute whose qualified name is `name`, as the
    /// source's `XMLAttributes::value` returns it, or `None` when absent.
    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.0
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
    }
}

/// One callback of the source's `XMLHandler`.
#[derive(Debug)]
pub(crate) enum SaxEvent<'a> {
    /// `onStartElement(qname, attributes)`.
    Start(&'a str, &'a SaxAttributes),
    /// `onEndElement(qname)`.
    End(&'a str),
    /// `onCharacters(chars)`: one run of character data.
    Characters(&'a str),
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("XML {what} limit exceeded"))
}

const INPUT_LIMIT: &str = "XML input byte limit exceeded";

/// A reader that refuses, rather than truncates, input beyond its allowance.
struct Bounded<R> {
    inner: R,
    remaining: u64,
    exceeded: bool,
}

impl<R: BufRead> Read for Bounded<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let n = {
            let available = self.fill_buf()?;
            let n = available.len().min(out.len());
            out[..n].copy_from_slice(&available[..n]);
            n
        };
        self.consume(n);
        Ok(n)
    }
}

impl<R: BufRead> BufRead for Bounded<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let available = self.inner.fill_buf()?;
        if available.is_empty() {
            return Ok(available);
        }
        if self.remaining == 0 {
            self.exceeded = true;
            return Err(io::Error::other(INPUT_LIMIT));
        }
        let n = usize::try_from(self.remaining)
            .unwrap_or(usize::MAX)
            .min(available.len());
        Ok(&available[..n])
    }

    fn consume(&mut self, amount: usize) {
        self.remaining = self.remaining.saturating_sub(amount as u64);
        self.inner.consume(amount);
    }
}

/// Normalise a raw attribute value (XML 1.0 section 3.3.3), then expand its
/// references.
fn attribute_value(raw: &[u8]) -> Result<String> {
    let text =
        std::str::from_utf8(raw).map_err(|_| parse_error("XML attribute value is not UTF-8"))?;
    let mut normalized = String::with_capacity(text.len());
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\r' => {
                if chars.peek() == Some(&'\n') {
                    chars.next();
                }
                normalized.push(' ');
            }
            '\t' | '\n' => normalized.push(' '),
            _ => normalized.push(c),
        }
    }
    let value = quick_xml::escape::unescape(&normalized)
        .map_err(|error| parse_error(format!("XML attribute value: {error}")))?;
    Ok(value.into_owned())
}

fn attributes(start: &BytesStart<'_>, limits: &SaxLimits) -> Result<SaxAttributes> {
    let mut out = Vec::new();
    for attribute in start.attributes().with_checks(true) {
        if out.len() == limits.max_attributes {
            return Err(limit("attribute"));
        }
        let attribute =
            attribute.map_err(|error| parse_error(format!("XML attribute: {error}")))?;
        let key = std::str::from_utf8(attribute.key.as_ref())
            .map_err(|_| parse_error("XML attribute name is not UTF-8"))?
            .to_owned();
        let value = attribute_value(&attribute.value)?;
        out.push((key, value));
    }
    Ok(SaxAttributes(out))
}

fn expand_reference(name: &str) -> Result<Cow<'static, str>> {
    Ok(Cow::Borrowed(match name {
        "amp" => "&",
        "lt" => "<",
        "gt" => ">",
        "quot" => "\"",
        "apos" => "'",
        _ => {
            let digits = name
                .strip_prefix('#')
                .ok_or_else(|| parse_error(format!("undeclared XML entity `&{name};`")))?;
            let code = match digits.strip_prefix('x') {
                Some(hex) => u32::from_str_radix(hex, 16),
                None => digits.parse::<u32>(),
            }
            .map_err(|_| parse_error(format!("malformed character reference `&{name};`")))?;
            let character = char::from_u32(code)
                .filter(|&c| c != '\0')
                .ok_or_else(|| parse_error(format!("invalid character reference `&{name};`")))?;
            return Ok(Cow::Owned(character.to_string()));
        }
    }))
}

/// Strip a UTF-8 byte order mark, which Xerces consumes silently.
fn skip_bom(input: &mut impl BufRead) -> Result<()> {
    let bom = input.fill_buf()?.starts_with(&[0xEF, 0xBB, 0xBF]);
    if bom {
        input.consume(3);
    }
    Ok(())
}

/// Drive `handler` with the callbacks of one complete document.
///
/// The handler's own errors stop the parse and are returned unchanged.
pub(crate) fn parse(
    input: impl BufRead,
    limits: &SaxLimits,
    mut handler: impl FnMut(SaxEvent<'_>) -> Result<()>,
) -> Result<()> {
    let mut bounded = Bounded {
        inner: input,
        remaining: limits.max_input_bytes,
        exceeded: false,
    };
    skip_bom(&mut bounded).map_err(|error| input_error(error, &bounded))?;
    let mut reader = Reader::from_reader(bounded);
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().check_end_names = true;
    reader.config_mut().check_comments = true;
    let mut buffer = Vec::new();
    let mut text = String::new();
    let mut depth = 0usize;
    let mut root_seen = false;
    let mut events = 0u64;
    loop {
        buffer.clear();
        let event = match reader.read_event_into(&mut buffer) {
            Ok(event) => event,
            Err(error) => {
                if reader.get_ref().exceeded {
                    return Err(Error::InvalidValue(INPUT_LIMIT.into()));
                }
                return Err(match error {
                    quick_xml::Error::Io(error) => {
                        Error::Io(io::Error::new(error.kind(), error.to_string()))
                    }
                    other => parse_error(format!("XML is not well formed: {other}")),
                });
            }
        };
        events += 1;
        if events > limits.max_events {
            return Err(limit("event"));
        }
        let is_text = matches!(event, Event::Text(_) | Event::GeneralRef(_));
        if !is_text && !text.is_empty() {
            handler(SaxEvent::Characters(&text))?;
            text.clear();
        }
        match event {
            Event::Start(start) => {
                if depth == 0 {
                    if root_seen {
                        return Err(parse_error("XML document has more than one root element"));
                    }
                    root_seen = true;
                }
                if depth == limits.max_depth {
                    return Err(limit("nesting depth"));
                }
                depth += 1;
                let name = std::str::from_utf8(start.name().as_ref())
                    .map_err(|_| parse_error("XML element name is not UTF-8"))?
                    .to_owned();
                let attributes = attributes(&start, limits)?;
                handler(SaxEvent::Start(&name, &attributes))?;
            }
            Event::End(end) => {
                // check_end_names pairs every end tag with its start tag.
                depth = depth.saturating_sub(1);
                let name = std::str::from_utf8(end.name().as_ref())
                    .map_err(|_| parse_error("XML element name is not UTF-8"))?
                    .to_owned();
                handler(SaxEvent::End(&name))?;
            }
            Event::Empty(_) => {
                // expand_empty_elements turns every empty tag into Start + End.
                return Err(parse_error("unexpected empty-element event"));
            }
            Event::Text(body) => {
                let decoded = body
                    .xml10_content()
                    .map_err(|error| parse_error(format!("XML text: {error}")))?;
                if depth == 0 {
                    if decoded
                        .chars()
                        .all(|c| matches!(c, ' ' | '\t' | '\n' | '\r'))
                    {
                        continue;
                    }
                    return Err(parse_error("character data outside the XML root element"));
                }
                push_text(&mut text, &decoded, limits)?;
            }
            Event::GeneralRef(reference) => {
                if depth == 0 {
                    return Err(parse_error("XML reference outside the root element"));
                }
                let name = std::str::from_utf8(reference.as_ref())
                    .map_err(|_| parse_error("XML reference is not UTF-8"))?;
                let expanded = expand_reference(name)?;
                push_text(&mut text, &expanded, limits)?;
            }
            Event::CData(body) => {
                if depth == 0 {
                    return Err(parse_error("CDATA outside the XML root element"));
                }
                let decoded = body
                    .xml_content()
                    .map_err(|error| parse_error(format!("XML CDATA: {error}")))?;
                if decoded.len() > limits.max_text_bytes {
                    return Err(limit("character data"));
                }
                handler(SaxEvent::Characters(&decoded))?;
            }
            Event::Decl(declaration) => {
                if let Some(encoding) = declaration.encoding() {
                    let encoding = encoding
                        .map_err(|error| parse_error(format!("XML declaration: {error}")))?;
                    if !encoding.eq_ignore_ascii_case(b"UTF-8")
                        && !encoding.eq_ignore_ascii_case(b"US-ASCII")
                    {
                        return Err(Error::Unsupported(format!(
                            "XML encoding {} (only UTF-8 is read)",
                            String::from_utf8_lossy(&encoding)
                        )));
                    }
                }
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported(
                    "XML document type declarations are not supported".into(),
                ));
            }
            Event::Comment(_) | Event::PI(_) => {}
            Event::Eof => {
                if !root_seen {
                    return Err(parse_error("XML document has no root element"));
                }
                if depth != 0 {
                    return Err(parse_error("XML document ends inside an element"));
                }
                return Ok(());
            }
        }
    }
}

fn push_text(text: &mut String, piece: &str, limits: &SaxLimits) -> Result<()> {
    if text.len().saturating_add(piece.len()) > limits.max_text_bytes {
        return Err(limit("character data"));
    }
    text.push_str(piece);
    Ok(())
}

fn input_error<R>(error: Error, bounded: &Bounded<R>) -> Error {
    if bounded.exceeded {
        Error::InvalidValue(INPUT_LIMIT.into())
    } else {
        error
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIMITS: SaxLimits = SaxLimits {
        max_input_bytes: 1 << 20,
        max_depth: 16,
        max_events: 10_000,
        max_attributes: 8,
        max_text_bytes: 1 << 16,
    };

    fn collect(text: &str, limits: &SaxLimits) -> Result<Vec<String>> {
        let mut out = Vec::new();
        parse(text.as_bytes(), limits, |event| {
            out.push(match event {
                SaxEvent::Start(name, attributes) => format!("S {name} {:?}", attributes.0),
                SaxEvent::End(name) => format!("E {name}"),
                SaxEvent::Characters(text) => format!("C {text:?}"),
            });
            Ok(())
        })?;
        Ok(out)
    }

    #[test]
    fn delivers_qualified_names_and_coalesced_text() {
        let events = collect(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<a:r x=\"1&#10;2\t3\r\n4\"><b>x&amp;y&#65;</b><c/>\
             <!-- c -->t<![CDATA[<z>]]></a:r>\n",
            &LIMITS,
        )
        .unwrap();
        assert_eq!(
            events,
            [
                "S a:r [(\"x\", \"1\\n2 3 4\")]",
                "S b []",
                "C \"x&yA\"",
                "E b",
                "S c []",
                "E c",
                "C \"t\"",
                "C \"<z>\"",
                "E a:r",
            ]
        );
    }

    #[test]
    fn refuses_malformed_and_unsupported_documents() {
        for (text, parse_failure) in [
            ("<a><b></a>", true),
            ("<a/><b/>", true),
            ("<a/>x", true),
            ("<a>&foo;</a>", true),
            ("<a x=\"1\" x=\"2\"/>", true),
            ("", true),
            ("<a>", true),
            ("<!DOCTYPE a><a/>", false),
            ("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><a/>", false),
        ] {
            let error = collect(text, &LIMITS).unwrap_err();
            assert_eq!(
                matches!(error, Error::Parse { .. }),
                parse_failure,
                "{text}: {error}"
            );
        }
    }

    #[test]
    fn bounds_refuse_instead_of_truncating() {
        let document = "<a><b>text</b></a>";
        let tight = |f: fn(&mut SaxLimits)| {
            let mut limits = LIMITS;
            f(&mut limits);
            collect(document, &limits).unwrap_err()
        };
        for error in [
            tight(|l| l.max_input_bytes = 5),
            tight(|l| l.max_depth = 1),
            tight(|l| l.max_events = 3),
            tight(|l| l.max_text_bytes = 3),
        ] {
            assert!(matches!(error, Error::InvalidValue(_)), "{error}");
        }
        assert!(
            collect(document, &LIMITS).is_ok(),
            "the same document parses within the default bounds"
        );
        let mut exact = LIMITS;
        exact.max_input_bytes = document.len() as u64;
        assert!(collect(document, &exact).is_ok());
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded SAX-style event driver shared by the X! Tandem and OMSSA XML
//! readers.
//!
//! The source readers (`XTandemXMLFile`, `OMSSAXMLFile`) are Xerces SAX2
//! handlers derived from `Internal::XMLHandler` (`FORMAT/HANDLERS/XMLHandler.h`):
//! they keep state across `onStartElement`, `onEndElement` and `onCharacters`
//! callbacks, and several of their rules depend on *which* character events
//! arrive, including the white space between elements (OMSSA resets the
//! current modification on the white space inside `<MSHits_mods>`). This
//! driver therefore reports every text node, white space included, as one
//! [`SaxHandler::characters`] call, with the five predefined entities and
//! numeric character references expanded and CDATA sections appended as
//! text. Xerces may split one text node into several calls at buffer or
//! entity boundaries; the source handlers act on the first call only, so a
//! split node would lose its tail there, and delivering it whole is the
//! intended reading.
//!
//! Element and attribute names are reported as written, prefix included, as
//! Xerces reports the qualified name. A DTD is refused ([`Error::Unsupported`]),
//! because its entities would need expanding. Every input is bounded by
//! [`XmlReadLimits`], checked before the corresponding storage grows.

use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::Event;
use std::io::Read;

/// Resource ceilings for one search-engine XML document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct XmlReadLimits {
    /// Maximum input size in bytes.
    pub max_bytes: usize,
    /// Maximum XML events, text and markup alike.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum bytes of one text node after entity expansion.
    pub max_text_bytes: usize,
    /// Maximum attributes on one element.
    pub max_attributes: usize,
}

impl Default for XmlReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1024 * 1024 * 1024,
            max_events: 200_000_000,
            max_depth: 256,
            max_text_bytes: 64 * 1024 * 1024,
            max_attributes: 1024,
        }
    }
}

impl XmlReadLimits {
    fn validate(&self) -> Result<()> {
        if self.max_bytes == 0
            || self.max_events == 0
            || self.max_depth == 0
            || self.max_text_bytes == 0
            || self.max_attributes == 0
        {
            return Err(Error::InvalidValue(
                "search-engine XML limits must be positive".into(),
            ));
        }
        Ok(())
    }
}

/// How the document bytes become text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TextEncoding {
    /// UTF-8, with or without a byte-order mark; a declaration naming any
    /// other encoding is refused. Xerces honours the declaration, so this
    /// covers the documents the search engines write.
    Utf8,
    /// ISO-8859-1 whatever the declaration says, as the source's
    /// `enforceEncoding_("ISO-8859-1")` forces (X! Tandem). Every byte is
    /// one character, so no input is refused for its encoding.
    Latin1,
}

/// The attributes of one start tag, in document order.
#[derive(Clone, Debug, Default)]
pub(crate) struct Attributes {
    entries: Vec<(String, String)>,
}

impl Attributes {
    /// The value of the attribute with this qualified name, if present.
    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.entries
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
    }

    /// The value of a required attribute.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the source's `fatalError` message "Required
    /// attribute 'name' not present!" when it is missing.
    pub(crate) fn required(&self, name: &str) -> Result<&str> {
        self.get(name).ok_or_else(|| Error::Parse {
            line: 0,
            message: format!("Required attribute '{name}' not present!"),
        })
    }
}

/// The callbacks of one SAX-style handler.
///
/// An error returned with `line` 0 gets the line of the event that caused it.
pub(crate) trait SaxHandler {
    /// An element opened.
    fn start_element(&mut self, name: &str, attributes: &Attributes) -> Result<()>;
    /// An element closed.
    fn end_element(&mut self, name: &str) -> Result<()>;
    /// One text node inside the root element, white space included.
    fn characters(&mut self, text: &str) -> Result<()>;
}

/// Read a whole document as text, refusing more than `limits.max_bytes`.
pub(crate) fn read_document(
    input: impl Read,
    encoding: TextEncoding,
    limits: &XmlReadLimits,
) -> Result<String> {
    limits.validate()?;
    let count = u64::try_from(limits.max_bytes)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| Error::InvalidValue("XML byte limit overflows".into()))?;
    let mut bytes = Vec::new();
    input.take(count).read_to_end(&mut bytes)?;
    if bytes.len() > limits.max_bytes {
        return Err(Error::InvalidValue("XML byte limit exceeded".into()));
    }
    match encoding {
        TextEncoding::Latin1 => Ok(bytes.iter().map(|&byte| char::from(byte)).collect()),
        TextEncoding::Utf8 => {
            let body = bytes
                .strip_prefix(&[0xEF, 0xBB, 0xBF])
                .map_or(bytes.as_slice(), |rest| rest);
            String::from_utf8(body.to_vec()).map_err(|_| Error::Parse {
                line: 0,
                message: "XML input is not valid UTF-8".into(),
            })
        }
    }
}

/// One-based line of a byte offset.
fn line_at(text: &str, offset: usize) -> usize {
    let end = offset.min(text.len());
    text.as_bytes()[..end]
        .iter()
        .filter(|&&byte| byte == b'\n')
        .count()
        + 1
}

fn syntax(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// Give a line-less parse error the current line.
fn locate(error: Error, text: &str, offset: usize) -> Error {
    match error {
        Error::Parse { line: 0, message } => Error::Parse {
            line: line_at(text, offset),
            message,
        },
        other => other,
    }
}

/// Expand one predefined or numeric entity reference.
fn expand_reference(raw: &[u8]) -> Result<String> {
    let name = std::str::from_utf8(raw).map_err(|_| syntax("non-UTF-8 XML entity reference"))?;
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
    let code = match digits.strip_prefix('x') {
        Some(hex) => u32::from_str_radix(hex, 16),
        None => digits.parse::<u32>(),
    }
    .map_err(|_| syntax("malformed numeric character reference"))?;
    char::from_u32(code)
        .map(String::from)
        .ok_or_else(|| syntax("numeric character reference is not a Unicode scalar value"))
}

struct Pending {
    text: String,
    offset: usize,
}

impl Pending {
    fn push(&mut self, fragment: &str, offset: usize, limits: &XmlReadLimits) -> Result<()> {
        if self.text.is_empty() {
            self.offset = offset;
        }
        if self.text.len().saturating_add(fragment.len()) > limits.max_text_bytes {
            return Err(Error::InvalidValue("XML text node limit exceeded".into()));
        }
        self.text.push_str(fragment);
        Ok(())
    }

    fn flush(&mut self, handler: &mut impl SaxHandler, document: &str) -> Result<()> {
        if self.text.is_empty() {
            return Ok(());
        }
        let text = std::mem::take(&mut self.text);
        handler
            .characters(&text)
            .map_err(|error| locate(error, document, self.offset))
    }
}

/// Drive `handler` through one complete document.
///
/// # Errors
///
/// [`Error::Parse`] when the document is not well formed (unbalanced tags,
/// several or no root elements, text outside the root, a malformed reference
/// or attribute), [`Error::Unsupported`] for a DTD or an encoding other than
/// UTF-8 in [`TextEncoding::Utf8`] mode, [`Error::InvalidValue`] when a limit
/// is exceeded, and any error the handler returns, which ends the parse.
pub(crate) fn drive(
    text: &str,
    encoding: TextEncoding,
    handler: &mut impl SaxHandler,
    limits: &XmlReadLimits,
) -> Result<()> {
    limits.validate()?;
    let mut reader = Reader::from_str(text);
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().check_end_names = true;
    let mut open: Vec<String> = Vec::new();
    let mut events = 0usize;
    let mut root_seen = false;
    let mut pending = Pending {
        text: String::new(),
        offset: 0,
    };
    loop {
        let offset = usize::try_from(reader.buffer_position()).unwrap_or(usize::MAX);
        let event = reader.read_event().map_err(|error| Error::Parse {
            line: line_at(text, usize::try_from(reader.error_position()).unwrap_or(0)),
            message: format!("XML is not well formed: {error}"),
        })?;
        events += 1;
        if events > limits.max_events {
            return Err(Error::InvalidValue("XML event limit exceeded".into()));
        }
        match event {
            Event::Start(start) => {
                pending.flush(handler, text)?;
                if open.is_empty() {
                    if root_seen {
                        return Err(locate(
                            syntax("XML holds more than one root element"),
                            text,
                            offset,
                        ));
                    }
                    root_seen = true;
                }
                if open.len() >= limits.max_depth {
                    return Err(Error::InvalidValue(
                        "XML nesting depth limit exceeded".into(),
                    ));
                }
                let name = std::str::from_utf8(start.name().as_ref())
                    .map_err(|_| locate(syntax("XML element name is not UTF-8"), text, offset))?
                    .to_owned();
                let mut attributes = Attributes::default();
                for attribute in start.attributes() {
                    let attribute = attribute.map_err(|error| {
                        locate(syntax(format!("XML attribute: {error}")), text, offset)
                    })?;
                    if attributes.entries.len() >= limits.max_attributes {
                        return Err(Error::InvalidValue("XML attribute limit exceeded".into()));
                    }
                    let key = std::str::from_utf8(attribute.key.as_ref())
                        .map_err(|_| {
                            locate(syntax("XML attribute name is not UTF-8"), text, offset)
                        })?
                        .to_owned();
                    let value = attribute
                        .unescape_value()
                        .map_err(|error| {
                            locate(syntax(format!("XML attribute: {error}")), text, offset)
                        })?
                        .into_owned();
                    attributes.entries.push((key, value));
                }
                handler
                    .start_element(&name, &attributes)
                    .map_err(|error| locate(error, text, offset))?;
                open.push(name);
            }
            Event::End(end) => {
                pending.flush(handler, text)?;
                let name = std::str::from_utf8(end.name().as_ref())
                    .map_err(|_| locate(syntax("XML element name is not UTF-8"), text, offset))?
                    .to_owned();
                if open.pop().is_none() {
                    return Err(locate(
                        syntax(format!("closing tag {name} not matched by an opening tag")),
                        text,
                        offset,
                    ));
                }
                handler
                    .end_element(&name)
                    .map_err(|error| locate(error, text, offset))?;
            }
            Event::Text(body) => {
                let decoded = body
                    .xml_content()
                    .map_err(|error| locate(syntax(format!("XML text: {error}")), text, offset))?;
                if open.is_empty() {
                    if decoded.trim_matches([' ', '\t', '\r', '\n']).is_empty() {
                        continue;
                    }
                    return Err(locate(
                        syntax("XML holds character data outside the root element"),
                        text,
                        offset,
                    ));
                }
                pending.push(&decoded, offset, limits)?;
            }
            Event::GeneralRef(reference) => {
                if open.is_empty() {
                    return Err(locate(
                        syntax("XML entity reference outside the root element"),
                        text,
                        offset,
                    ));
                }
                let expanded = expand_reference(reference.as_ref())
                    .map_err(|error| locate(error, text, offset))?;
                pending.push(&expanded, offset, limits)?;
            }
            Event::CData(data) => {
                if open.is_empty() {
                    return Err(locate(
                        syntax("XML CDATA section outside the root element"),
                        text,
                        offset,
                    ));
                }
                let raw = std::str::from_utf8(data.as_ref())
                    .map_err(|_| locate(syntax("XML CDATA is not UTF-8"), text, offset))?;
                pending.push(raw, offset, limits)?;
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported("XML DTDs are not supported".into()));
            }
            Event::Decl(declaration) => {
                if encoding == TextEncoding::Utf8 {
                    if let Some(name) = declaration.encoding().transpose().map_err(|error| {
                        locate(syntax(format!("XML declaration: {error}")), text, offset)
                    })? {
                        if !name.eq_ignore_ascii_case(b"UTF-8") {
                            return Err(Error::Unsupported(
                                "this XML reader supports only UTF-8 documents".into(),
                            ));
                        }
                    }
                }
            }
            Event::Comment(_) | Event::PI(_) => pending.flush(handler, text)?,
            Event::Empty(_) => {
                return Err(syntax("unexpanded empty XML element"));
            }
            Event::Eof => {
                if !root_seen || !open.is_empty() {
                    return Err(Error::Parse {
                        line: line_at(text, text.len()),
                        message: "incomplete XML document".into(),
                    });
                }
                return Ok(());
            }
        }
    }
}

/// Xerces `XMLString::parseInt`, which the source's `attributeAsInt_` uses:
/// XML white space is trimmed, then an optionally signed decimal integer must
/// fill the rest and fit an `i32`.
pub(crate) fn xerces_parse_int(text: &str) -> Result<i32> {
    let trimmed = text.trim_matches([' ', '\t', '\r', '\n']);
    let digits = trimmed.strip_prefix(['+', '-']).unwrap_or(trimmed);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(syntax(format!("'{text}' is not an integer")));
    }
    trimmed
        .strip_prefix('+')
        .unwrap_or(trimmed)
        .parse::<i32>()
        .map_err(|_| syntax(format!("'{text}' is out of the integer range")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Recorder {
        events: Vec<String>,
    }

    impl SaxHandler for Recorder {
        fn start_element(&mut self, name: &str, attributes: &Attributes) -> Result<()> {
            self.events
                .push(format!("<{name} {}", attributes.get("a").unwrap_or("-")));
            Ok(())
        }
        fn end_element(&mut self, name: &str) -> Result<()> {
            self.events.push(format!("</{name}"));
            Ok(())
        }
        fn characters(&mut self, text: &str) -> Result<()> {
            self.events.push(format!("[{text}]"));
            Ok(())
        }
    }

    #[test]
    fn text_nodes_arrive_whole_with_white_space_and_references() {
        let mut recorder = Recorder::default();
        let text = "<?xml version=\"1.0\"?>\n<r a=\"1&amp;2\">\n <x>1&#46;5</x><y/></r>\n";
        drive(
            text,
            TextEncoding::Utf8,
            &mut recorder,
            &XmlReadLimits::default(),
        )
        .unwrap();
        assert_eq!(
            recorder.events,
            [
                "<r 1&2", "[\n ]", "<x -", "[1.5]", "</x", "<y -", "</y", "</r"
            ]
        );
    }

    #[test]
    fn malformed_documents_are_refused_with_a_line() {
        let mut recorder = Recorder::default();
        let error = drive(
            "<r>\n<x></y></r>",
            TextEncoding::Utf8,
            &mut recorder,
            &XmlReadLimits::default(),
        )
        .unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }), "{error:?}");
        for bad in ["", "<r/><s/>", "x<r/>", "<r>"] {
            assert!(
                drive(
                    bad,
                    TextEncoding::Utf8,
                    &mut Recorder::default(),
                    &XmlReadLimits::default()
                )
                .is_err(),
                "{bad}"
            );
        }
        assert!(matches!(
            drive(
                "<!DOCTYPE r><r/>",
                TextEncoding::Utf8,
                &mut Recorder::default(),
                &XmlReadLimits::default()
            ),
            Err(Error::Unsupported(_))
        ));
    }

    #[test]
    fn latin1_ignores_the_declared_encoding() {
        let text = read_document(
            &b"<?xml version=\"1.0\" encoding=\"UTF-16\"?><r>\xe9</r>"[..],
            TextEncoding::Latin1,
            &XmlReadLimits::default(),
        )
        .unwrap();
        let mut recorder = Recorder::default();
        drive(
            &text,
            TextEncoding::Latin1,
            &mut recorder,
            &XmlReadLimits::default(),
        )
        .unwrap();
        assert_eq!(recorder.events[1], "[\u{e9}]");
        assert!(
            drive(
                &text,
                TextEncoding::Utf8,
                &mut Recorder::default(),
                &XmlReadLimits::default()
            )
            .is_err()
        );
    }

    #[test]
    fn limits_are_enforced() {
        let limits = XmlReadLimits {
            max_depth: 1,
            ..XmlReadLimits::default()
        };
        assert!(matches!(
            drive(
                "<r><x/></r>",
                TextEncoding::Utf8,
                &mut Recorder::default(),
                &limits
            ),
            Err(Error::InvalidValue(_))
        ));
        let limits = XmlReadLimits {
            max_bytes: 3,
            ..XmlReadLimits::default()
        };
        assert!(read_document(&b"<r/>"[..], TextEncoding::Utf8, &limits).is_err());
    }

    #[test]
    fn xerces_integers() {
        assert_eq!(xerces_parse_int(" 12 ").unwrap(), 12);
        assert_eq!(xerces_parse_int("-3").unwrap(), -3);
        assert_eq!(xerces_parse_int("+3").unwrap(), 3);
        assert!(xerces_parse_int("3.0").is_err());
        assert!(xerces_parse_int("").is_err());
        assert!(xerces_parse_int("2147483648").is_err());
    }
}

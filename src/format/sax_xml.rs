// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Streaming, bounded SAX-style events with the callbacks Xerces delivers to an
//! OpenMS `Internal::XMLHandler`.
//!
//! The source parses its SAX-handled formats through `Internal::XMLFile::parse`
//! (`XMLFile.cpp:78-114`), a Xerces-C `SAX2XMLReader` with namespace processing
//! and namespace prefixes switched off, bridged to the handler's
//! `onStartElement`, `onEndElement` and `onCharacters`. This module is that
//! machinery as a pull loop over `quick-xml`, shared by every reader that ports
//! an `XMLHandler`: PTMXML, Unimod XML and UniProtKB XML (`kb-xml`), and the
//! X! Tandem, OMSSA and ProteinProphet protXML readers (`idxml`). Before wave-4
//! integration each of those lane groups carried its own driver, and all three
//! delivered character data differently from Xerces; this one follows the
//! Release build's Xerces, recorded by `../oracle/integ-F-xerces-sax` (see the
//! `release_callbacks` test).
//!
//! What a handler observes, as Xerces delivers it:
//!
//! - Element and attribute names are qualified names as written (`umod:mod`,
//!   not `mod`); nothing is resolved.
//! - Attribute values are normalised as XML 1.0 section 3.3.3 requires: a
//!   literal CR LF pair, CR, LF or tab becomes one space; references are
//!   expanded afterwards, so `&#10;` stays a line feed. A literal `<` in a
//!   value is a well-formedness error.
//! - Character data is split where Xerces splits it: the text before a
//!   reference is one callback, and the reference's replacement starts the
//!   next one, which runs to the next reference or markup. `1&#46;5` arrives as
//!   `1` and `.5`, `&lt;&gt;` as `<` and `>`. Markup, a comment or a processing
//!   instruction ends a run; a CDATA section is a callback of its own, even
//!   when empty. Line ends are normalised (XML 1.0 section 2.11).
//! - White space inside the root element is character data (there is no DTD to
//!   declare it ignorable); white space around the root element is not
//!   delivered.
//!
//! Every well-formedness error Xerces reports as a fatal error, which the
//! source turns into `Exception::ParseError` (mismatched or unclosed tags, a
//! second root, text outside the root, an undeclared entity, a character
//! reference to NUL or with an upper-case `X`, a duplicate attribute), is
//! [`Error::Parse`] with the line on which the parser stood, as Xerces reports
//! it. A handler error of kind [`Error::Parse`] with line 0 gets that line too,
//! the way `XMLHandler::fatalError` adds the locator's line.
//!
//! Native differences, each an explicit refusal rather than a guess:
//!
//! - A document type declaration is [`Error::Unsupported`]: its internal subset
//!   can declare entities, which Xerces expands and this reader cannot.
//! - In [`SaxEncoding::Utf8`] mode a declared encoding other than UTF-8 (or its
//!   subset US-ASCII) is [`Error::Unsupported`]; Xerces transcodes more.
//!   [`SaxEncoding::Declared`] honours ISO-8859-1 and UTF-16 as Xerces does, and
//!   [`SaxEncoding::Latin1`] reproduces the source's
//!   `XMLFile::enforceEncoding_("ISO-8859-1")`.
//! - Input bytes, nesting depth, events, attributes per element and the bytes
//!   of one run of character data are bounded by [`SaxLimits`]; exceeding a
//!   bound is [`Error::InvalidValue`], never a truncated document.
//! - Xerces may also split a run at an internal buffer boundary; a document
//!   held in memory is not split there (the oracle's 100,000-character run is
//!   one callback), and neither is a run here.

use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::io::{self, BufRead, Read};

/// Resource bounds for one parse. Each bound is checked before the work it
/// limits is done, so exceeding it never yields a partial result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SaxLimits {
    /// Input bytes the parser may consume (before decoding).
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

impl SaxLimits {
    fn validate(&self) -> Result<()> {
        if self.max_input_bytes == 0
            || self.max_depth == 0
            || self.max_events == 0
            || self.max_attributes == 0
            || self.max_text_bytes == 0
        {
            return Err(Error::InvalidValue("XML limits must be positive".into()));
        }
        Ok(())
    }
}

/// How the input bytes become text.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum SaxEncoding {
    /// UTF-8, with or without a byte-order mark; a declaration naming any
    /// encoding other than UTF-8 or US-ASCII is refused. Streams.
    #[default]
    Utf8,
    /// Every byte is one ISO-8859-1 character, whatever the declaration says:
    /// the source's `XMLFile::enforceEncoding_("ISO-8859-1")` (X! Tandem). The
    /// bounded input is decoded in memory first.
    #[cfg_attr(not(feature = "idxml"), allow(dead_code))]
    Latin1,
    /// The byte-order mark and declared encoding (UTF-8, UTF-16, ISO-8859-1),
    /// as Xerces honours them, through the identification-XML decoder. The
    /// bounded input is decoded in memory first.
    #[cfg(feature = "idxml")]
    Declared,
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

    /// The value of a required attribute (the source's
    /// `XMLHandler::attributeAsString_`).
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the source's `fatalError` message "Required
    /// attribute 'name' not present!" when it is missing; the driver adds the
    /// line.
    pub(crate) fn required(&self, name: &str) -> Result<&str> {
        self.get(name)
            .ok_or_else(|| parse_error(format!("Required attribute '{name}' not present!")))
    }
}

/// One callback of the source's `XMLHandler`.
#[derive(Debug)]
pub(crate) enum SaxEvent<'a> {
    /// `onStartElement(qname, attributes)`.
    Start(&'a str, &'a SaxAttributes),
    /// `onEndElement(qname)`.
    End(&'a str),
    /// `onCharacters(chars)`: one run of character data, split as Xerces
    /// splits it (see the module documentation).
    Characters(&'a str),
}

/// The callbacks of one handler, for readers that keep their state in a
/// struct rather than a closure.
#[cfg_attr(not(feature = "idxml"), allow(dead_code))]
pub(crate) trait SaxHandler {
    /// `onStartElement(qname, attributes)`.
    fn start_element(&mut self, name: &str, attributes: &SaxAttributes) -> Result<()>;
    /// `onEndElement(qname)`.
    fn end_element(&mut self, name: &str) -> Result<()>;
    /// `onCharacters(chars)`.
    fn characters(&mut self, text: &str) -> Result<()>;
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

/// A reader that refuses, rather than truncates, input beyond its allowance,
/// and counts the line feeds it has handed out, which is the line the parser
/// stands on.
struct Bounded<R> {
    inner: R,
    remaining: u64,
    exceeded: bool,
    newlines: usize,
}

impl<R> Bounded<R> {
    fn line(&self) -> usize {
        self.newlines.saturating_add(1)
    }
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
        if let Ok(available) = self.inner.fill_buf() {
            let counted = amount.min(available.len());
            self.newlines = self
                .newlines
                .saturating_add(available[..counted].iter().filter(|&&b| b == b'\n').count());
        }
        self.remaining = self.remaining.saturating_sub(amount as u64);
        self.inner.consume(amount);
    }
}

/// Whether `code` is an XML 1.0 `Char` (section 2.2); a character reference to
/// anything else is a well-formedness error, as Xerces reports for `&#0;`.
fn xml_char(code: u32) -> bool {
    matches!(code, 0x9 | 0xA | 0xD | 0x20..=0xD7FF | 0xE000..=0xFFFD | 0x10000..=0x10_FFFF)
}

/// Expand one reference (the name between `&` and `;`): the five predefined
/// entities and decimal or lower-case-`x` hexadecimal character references.
fn expand_reference(name: &str) -> Result<String> {
    let text = match name {
        "amp" => "&",
        "lt" => "<",
        "gt" => ">",
        "quot" => "\"",
        "apos" => "'",
        _ => {
            let Some(digits) = name.strip_prefix('#') else {
                return Err(parse_error(format!("entity '{name}' not found")));
            };
            let code = match digits.strip_prefix('x') {
                Some(hex) if !hex.is_empty() && hex.bytes().all(|b| b.is_ascii_hexdigit()) => {
                    u32::from_str_radix(hex, 16).ok()
                }
                Some(_) => None,
                None if digits.starts_with('X') => {
                    return Err(parse_error(
                        "hex radix character references must use 'x', not 'X'",
                    ));
                }
                None if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) => {
                    digits.parse::<u32>().ok()
                }
                None => None,
            };
            let character = code
                .filter(|&code| xml_char(code))
                .and_then(char::from_u32)
                .ok_or_else(|| parse_error("invalid character reference"))?;
            return Ok(character.to_string());
        }
    };
    Ok(text.to_owned())
}

/// Normalise a raw attribute value (XML 1.0 section 3.3.3), then expand its
/// references.
fn attribute_value(raw: &[u8], name: &str) -> Result<String> {
    let text =
        std::str::from_utf8(raw).map_err(|_| parse_error("XML attribute value is not UTF-8"))?;
    if text.contains('<') {
        return Err(parse_error(format!(
            "'<' character cannot be used in attribute value '{name}'; use &lt; instead"
        )));
    }
    let mut value = String::with_capacity(text.len());
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        match c {
            '\r' => {
                if chars.peek().map(|&(_, next)| next) == Some('\n') {
                    chars.next();
                }
                value.push(' ');
            }
            '\t' | '\n' => value.push(' '),
            '&' => {
                let rest = &text[at + 1..];
                let end = rest
                    .find(';')
                    .ok_or_else(|| parse_error("unterminated reference in XML attribute"))?;
                value.push_str(&expand_reference(&rest[..end])?);
                for _ in 0..=rest[..end].chars().count() {
                    chars.next();
                }
            }
            _ => value.push(c),
        }
    }
    Ok(value)
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
        let value = attribute_value(&attribute.value, &key)?;
        out.push((key, value));
    }
    Ok(SaxAttributes(out))
}

/// Strip a UTF-8 byte order mark, which Xerces consumes silently.
fn skip_bom(input: &mut impl BufRead) -> Result<()> {
    let bom = input.fill_buf()?.starts_with(&[0xEF, 0xBB, 0xBF]);
    if bom {
        input.consume(3);
    }
    Ok(())
}

/// Drive `handler` with the callbacks of one complete UTF-8 document.
///
/// # Errors
///
/// See [`parse_with`].
#[cfg_attr(not(feature = "kb-xml"), allow(dead_code))]
pub(crate) fn parse(
    input: impl BufRead,
    limits: &SaxLimits,
    handler: impl FnMut(SaxEvent<'_>) -> Result<()>,
) -> Result<()> {
    parse_with(input, SaxEncoding::Utf8, limits, handler)
}

/// Drive a [`SaxHandler`] with the callbacks of one complete document.
///
/// # Errors
///
/// See [`parse_with`].
#[cfg_attr(not(feature = "idxml"), allow(dead_code))]
pub(crate) fn drive(
    input: impl BufRead,
    encoding: SaxEncoding,
    limits: &SaxLimits,
    handler: &mut impl SaxHandler,
) -> Result<()> {
    parse_with(input, encoding, limits, |event| match event {
        SaxEvent::Start(name, attributes) => handler.start_element(name, attributes),
        SaxEvent::End(name) => handler.end_element(name),
        SaxEvent::Characters(text) => handler.characters(text),
    })
}

/// Drive `handler` with the callbacks of one complete document read in
/// `encoding`.
///
/// # Errors
///
/// [`Error::Parse`] (with the line) for a document that is not well formed;
/// [`Error::Unsupported`] for a DTD, or a declared encoding other than UTF-8 in
/// [`SaxEncoding::Utf8`] mode; [`Error::InvalidValue`] when a limit is zero or
/// exceeded; [`Error::Io`] when reading fails; and any error the handler
/// returns, which ends the parse (an [`Error::Parse`] with line 0 gets the
/// current line).
pub(crate) fn parse_with(
    input: impl BufRead,
    encoding: SaxEncoding,
    limits: &SaxLimits,
    handler: impl FnMut(SaxEvent<'_>) -> Result<()>,
) -> Result<()> {
    limits.validate()?;
    match encoding {
        SaxEncoding::Utf8 => run(input, limits, true, handler),
        SaxEncoding::Latin1 => {
            let text = read_latin1(input, limits.max_input_bytes)?;
            run(text.as_bytes(), &decoded(limits), false, handler)
        }
        #[cfg(feature = "idxml")]
        SaxEncoding::Declared => {
            let cap = usize::try_from(limits.max_input_bytes).unwrap_or(usize::MAX);
            let text = crate::format::identification_xml::decode_document(input, cap, None)?;
            run(text.as_bytes(), &decoded(limits), false, handler)
        }
    }
}

/// Limits for a document already decoded within `max_input_bytes`: its
/// decoded bytes (a Latin-1 byte can become two) are not counted again.
fn decoded(limits: &SaxLimits) -> SaxLimits {
    SaxLimits {
        max_input_bytes: u64::MAX,
        ..*limits
    }
}

fn read_latin1(input: impl Read, max_bytes: u64) -> Result<String> {
    let mut bytes = Vec::new();
    input
        .take(max_bytes.saturating_add(1))
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > max_bytes {
        return Err(Error::InvalidValue(INPUT_LIMIT.into()));
    }
    Ok(bytes.iter().map(|&byte| char::from(byte)).collect())
}

/// Give a line-less parse error the line the parser stands on.
fn locate(error: Error, line: usize) -> Error {
    match error {
        Error::Parse { line: 0, message } => Error::Parse { line, message },
        other => other,
    }
}

/// The current run of character data, delivered at the next reference or
/// markup.
struct Run {
    text: String,
}

impl Run {
    fn push(&mut self, piece: &str, limits: &SaxLimits) -> Result<()> {
        if self.text.len().saturating_add(piece.len()) > limits.max_text_bytes {
            return Err(limit("character data"));
        }
        self.text.push_str(piece);
        Ok(())
    }

    fn flush(&mut self, handler: &mut impl FnMut(SaxEvent<'_>) -> Result<()>) -> Result<()> {
        if self.text.is_empty() {
            return Ok(());
        }
        let result = handler(SaxEvent::Characters(&self.text));
        self.text.clear();
        result
    }
}

fn run(
    input: impl BufRead,
    limits: &SaxLimits,
    check_declared_encoding: bool,
    mut handler: impl FnMut(SaxEvent<'_>) -> Result<()>,
) -> Result<()> {
    let mut bounded = Bounded {
        inner: input,
        remaining: limits.max_input_bytes,
        exceeded: false,
        newlines: 0,
    };
    if let Err(error) = skip_bom(&mut bounded) {
        return Err(if bounded.exceeded {
            Error::InvalidValue(INPUT_LIMIT.into())
        } else {
            error
        });
    }
    let mut reader = Reader::from_reader(bounded);
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().check_end_names = true;
    reader.config_mut().check_comments = true;
    let mut buffer = Vec::new();
    let mut run = Run {
        text: String::new(),
    };
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
                // Xerces delivers the run before it scans the markup that
                // turns out to be malformed.
                let line = reader.get_ref().line();
                run.flush(&mut handler)
                    .map_err(|error| locate(error, line))?;
                return Err(match error {
                    quick_xml::Error::Io(error) => {
                        Error::Io(io::Error::new(error.kind(), error.to_string()))
                    }
                    other => Error::Parse {
                        line,
                        message: format!("XML is not well formed: {other}"),
                    },
                });
            }
        };
        let line = reader.get_ref().line();
        events += 1;
        if events > limits.max_events {
            return Err(limit("event"));
        }
        let step = match event {
            Event::Text(body) => (|| {
                let decoded = body
                    .xml10_content()
                    .map_err(|error| parse_error(format!("XML text: {error}")))?;
                if depth == 0 {
                    if decoded
                        .chars()
                        .all(|c| matches!(c, ' ' | '\t' | '\n' | '\r'))
                    {
                        return Ok(());
                    }
                    return Err(parse_error("character data outside the XML root element"));
                }
                run.push(&decoded, limits)
            })(),
            Event::GeneralRef(reference) => (|| {
                if depth == 0 {
                    return Err(parse_error("XML reference outside the root element"));
                }
                let name = std::str::from_utf8(reference.as_ref())
                    .map_err(|_| parse_error("XML reference is not UTF-8"))?;
                let expanded = expand_reference(name)?;
                // Xerces ends the run before a reference; the replacement
                // starts the next one.
                run.flush(&mut handler)?;
                run.push(&expanded, limits)
            })(),
            Event::CData(body) => (|| {
                if depth == 0 {
                    return Err(parse_error("CDATA outside the XML root element"));
                }
                run.flush(&mut handler)?;
                let decoded = body
                    .xml10_content()
                    .map_err(|error| parse_error(format!("XML CDATA: {error}")))?;
                if decoded.len() > limits.max_text_bytes {
                    return Err(limit("character data"));
                }
                handler(SaxEvent::Characters(&decoded))
            })(),
            Event::Start(start) => (|| {
                run.flush(&mut handler)?;
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
                handler(SaxEvent::Start(&name, &attributes))
            })(),
            Event::End(end) => (|| {
                run.flush(&mut handler)?;
                // check_end_names pairs every end tag with its start tag.
                depth = depth.saturating_sub(1);
                let name = std::str::from_utf8(end.name().as_ref())
                    .map_err(|_| parse_error("XML element name is not UTF-8"))?
                    .to_owned();
                handler(SaxEvent::End(&name))
            })(),
            Event::Empty(_) => {
                // expand_empty_elements turns every empty tag into Start + End.
                Err(parse_error("unexpected empty-element event"))
            }
            Event::Decl(declaration) => (|| {
                if !check_declared_encoding {
                    return Ok(());
                }
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
                Ok(())
            })(),
            Event::DocType(_) => Err(Error::Unsupported(
                "XML document type declarations are not supported".into(),
            )),
            Event::Comment(_) | Event::PI(_) => run.flush(&mut handler),
            Event::Eof => {
                if !root_seen {
                    return Err(Error::Parse {
                        line,
                        message: "XML document has no root element".into(),
                    });
                }
                if depth != 0 {
                    return Err(Error::Parse {
                        line,
                        message: "XML document ends inside an element".into(),
                    });
                }
                return Ok(());
            }
        };
        step.map_err(|error| locate(error, line))?;
    }
}

/// Xerces `XMLString::parseInt`, which the source's `attributeAsInt_` uses:
/// the text is trimmed of XML white space, then `strtol` in base 10 must
/// consume all of it; `ERANGE` (outside a 64-bit `long`) and unconsumed
/// characters are `NumberFormatException`s, which end the parse. Xerces then
/// casts the `long` to `int`, which wraps a value outside the 32-bit range;
/// that is refused here instead. The trimming and `strtol` rules are
/// `StringManager`'s (`format::string_manager`), which lane F7b replayed
/// against the Release build (`strtol` also skips a leading vertical tab or
/// form feed after the XML white-space trim).
///
/// # Errors
///
/// [`Error::Parse`] when the text is not such an integer or overflows a
/// `long`; [`Error::InvalidValue`] when it fits a `long` but not an `i32`.
#[cfg_attr(not(feature = "idxml"), allow(dead_code))]
pub(crate) fn xerces_parse_int(text: &str) -> Result<i32> {
    let units = crate::format::string_manager::StringManager::convert(text);
    match crate::format::string_manager::parse_long(&units) {
        Ok(wide) => i32::try_from(wide).map_err(|_| {
            Error::InvalidValue(format!(
                "'{text}' is outside the 32-bit range the source casts it to"
            ))
        }),
        Err(Error::InvalidRange(_)) => Err(parse_error(format!("'{text}' overflows a long"))),
        Err(_) => Err(parse_error(format!(
            "Could not convert '{text}' to an integer"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIMITS: SaxLimits = SaxLimits {
        max_input_bytes: 1 << 20,
        max_depth: 16,
        max_events: 1_000_000,
        max_attributes: 8,
        max_text_bytes: 1 << 20,
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

    // From lane F6's driver. Its expectation read `C "x&yA"`: the Release
    // build's Xerces delivers this very document's `x&amp;y&#65;` as three
    // callbacks (oracle case cdata_comment), so it now reads three.
    #[test]
    fn delivers_qualified_names_and_splits_at_references() {
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
                "C \"x\"",
                "C \"&y\"",
                "C \"A\"",
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
            ("x<a/>", true),
            ("<a>&foo;</a>", true),
            ("<a x=\"1\" x=\"2\"/>", true),
            ("", true),
            ("<a>", true),
            ("<a>&#0;</a>", true),
            ("<a>&#X41;</a>", true),
            ("<a x=\"a<b\"/>", true),
            ("<a x=\"&#0;\"/>", true),
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

    // From lane F3b's driver: a mismatched end tag on line 2 is reported there.
    #[test]
    fn malformed_documents_are_refused_with_a_line() {
        let error = collect("<r>\n<x></y></r>", &LIMITS).unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }), "{error:?}");
        let error = collect("<r>\n<x>\n</y></r>", &LIMITS).unwrap_err();
        assert!(matches!(error, Error::Parse { line: 3, .. }), "{error:?}");
    }

    #[test]
    fn handler_errors_get_the_line() {
        let text = "<r>\n\n<x/></r>";
        let error = parse(text.as_bytes(), &LIMITS, |event| match event {
            SaxEvent::Start("x", attributes) => attributes.required("id").map(|_| ()),
            _ => Ok(()),
        })
        .unwrap_err();
        match error {
            Error::Parse { line, message } => {
                assert_eq!(line, 3);
                assert_eq!(message, "Required attribute 'id' not present!");
            }
            other => panic!("{other:?}"),
        }
    }

    // From lane F3b's driver: X! Tandem's forced ISO-8859-1 ignores the
    // declaration, which UTF-8 mode refuses.
    #[test]
    fn latin1_ignores_the_declared_encoding() {
        let bytes = &b"<?xml version=\"1.0\" encoding=\"UTF-16\"?><r>\xe9</r>"[..];
        let mut out = Vec::new();
        parse_with(bytes, SaxEncoding::Latin1, &LIMITS, |event| {
            if let SaxEvent::Characters(text) = event {
                out.push(text.to_owned());
            }
            Ok(())
        })
        .unwrap();
        assert_eq!(out, ["\u{e9}"]);
        assert!(matches!(
            parse(bytes, &LIMITS, |_| Ok(())),
            Err(Error::Unsupported(_))
        ));
    }

    #[cfg(feature = "idxml")]
    #[test]
    fn declared_latin1_is_honoured() {
        // Oracle case latin1_declared: Xerces reads the byte 0xE9 as U+00E9.
        let bytes = &b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><r>\xe9</r>"[..];
        let mut out = Vec::new();
        parse_with(bytes, SaxEncoding::Declared, &LIMITS, |event| {
            if let SaxEvent::Characters(text) = event {
                out.push(text.to_owned());
            }
            Ok(())
        })
        .unwrap();
        assert_eq!(out, ["\u{e9}"]);
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
            tight(|l| l.max_attributes = 0),
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
        let mut latin1 = LIMITS;
        latin1.max_input_bytes = 3;
        assert!(matches!(
            parse_with(&b"<r/>"[..], SaxEncoding::Latin1, &latin1, |_| Ok(())),
            Err(Error::InvalidValue(_))
        ));
    }

    // From lane F3b's driver.
    #[test]
    fn xerces_integers() {
        assert_eq!(xerces_parse_int(" 12 ").unwrap(), 12);
        assert_eq!(xerces_parse_int("-3").unwrap(), -3);
        assert_eq!(xerces_parse_int("+3").unwrap(), 3);
        assert!(xerces_parse_int("3.0").is_err());
        assert!(xerces_parse_int("").is_err());
        assert!(xerces_parse_int("2147483648").is_err());
    }

    // From lane F3a's protXML reader (its own copy of parseInt, now this one).
    #[test]
    fn xerces_int_follows_strtol_on_the_trimmed_text() {
        assert_eq!(xerces_parse_int(" +2\n").unwrap(), 2);
        assert_eq!(xerces_parse_int("-3").unwrap(), -3);
        for text in ["", " ", "two", "2.0", "0x2", "+-2", "2 3"] {
            assert!(
                matches!(xerces_parse_int(text), Err(Error::Parse { .. })),
                "{text:?}"
            );
        }
        assert!(matches!(
            xerces_parse_int("4294967298"),
            Err(Error::InvalidValue(_))
        ));
        assert!(matches!(
            xerces_parse_int("99999999999999999999"),
            Err(Error::Parse { .. })
        ));
    }

    // From lane F3b's driver (`drive`, text as one string): the same document
    // is the oracle's refs_whitespace case, whose Release callbacks deliver
    // `1&#46;5` as `1` and `.5`; F3b expected one `[1.5]`.
    #[test]
    fn text_nodes_arrive_with_white_space_and_split_at_references() {
        #[derive(Default)]
        struct Recorder {
            events: Vec<String>,
        }
        impl SaxHandler for Recorder {
            fn start_element(&mut self, name: &str, attributes: &SaxAttributes) -> Result<()> {
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
        let mut recorder = Recorder::default();
        let text = "<?xml version=\"1.0\"?>\n<r a=\"1&amp;2\">\n <x>1&#46;5</x><y/></r>\n";
        drive(text.as_bytes(), SaxEncoding::Utf8, &LIMITS, &mut recorder).unwrap();
        assert_eq!(
            recorder.events,
            [
                "<r 1&2", "[\n ]", "<x -", "[1]", "[.5]", "</x", "<y -", "</y", "</r"
            ]
        );
    }

    // From lane F3b's driver.
    #[test]
    fn limits_are_enforced() {
        let limits = SaxLimits {
            max_depth: 1,
            ..LIMITS
        };
        assert!(matches!(
            collect("<r><x/></r>", &limits),
            Err(Error::InvalidValue(_))
        ));
        let limits = SaxLimits {
            max_input_bytes: 3,
            ..LIMITS
        };
        assert!(collect("<r/>", &limits).is_err());
    }

    // From lane F3a's driver, whose events were `'text` pieces: the Release
    // callbacks for this document (oracle case qualified_normalised) are the
    // same two, `t` and `A`.
    #[test]
    fn names_are_qualified_and_values_normalised() {
        let got = collect("<a:r x=\"1&amp;\t2\r\n3\"><b/>t&#65;</a:r>", &LIMITS).unwrap();
        assert_eq!(
            got,
            [
                "S a:r [(\"x\", \"1& 2 3\")]",
                "S b []",
                "E b",
                "C \"t\"",
                "C \"A\"",
                "E a:r"
            ]
        );
    }

    // From lane F3a's driver.
    #[test]
    fn malformed_documents_are_refused() {
        let limits = SaxLimits {
            max_depth: 8,
            max_attributes: 8,
            ..LIMITS
        };
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
            assert!(collect(text, &limits).is_err(), "{text:?}");
        }
        let deep = "<a>".repeat(9) + &"</a>".repeat(9);
        assert!(matches!(
            collect(&deep, &limits),
            Err(Error::InvalidValue(_))
        ));
        let wide = format!(
            "<a {}/>",
            (0..9)
                .map(|i| format!("x{i}=\"\""))
                .collect::<Vec<_>>()
                .join(" ")
        );
        assert!(matches!(
            collect(&wide, &limits),
            Err(Error::InvalidValue(_))
        ));
    }

    /// The Release build's Xerces callbacks for 28 documents
    /// (`../oracle/integ-F-xerces-sax`, ibminode06, run twice, identical),
    /// rendered as the oracle driver renders them.
    const RELEASE: &str =
        include_str!("../../tests/data/xerces_sax_release/cpp_release_records.txt");

    /// The oracle's documents, built exactly as `driver.cpp` builds them.
    fn release_documents() -> Vec<(&'static str, Vec<u8>)> {
        let big = "a".repeat(100_000);
        let big_refs = "ab&amp;".repeat(20_000);
        let text: Vec<(&'static str, String)> = vec![
            (
                "refs_whitespace",
                "<?xml version=\"1.0\"?>\n<r a=\"1&amp;2\">\n <x>1&#46;5</x><y/></r>\n".into(),
            ),
            (
                "qualified_normalised",
                "<a:r x=\"1&amp;\t2\r\n3\"><b/>t&#65;</a:r>".into(),
            ),
            (
                "cdata_comment",
                "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<a:r x=\"1&#10;2\t3\r\n4\"><b>x&amp;y&#65;</b><c/><!-- c -->t<![CDATA[<z>]]></a:r>\n".into(),
            ),
            ("text_comment_text", "<r>a<!--c-->b</r>".into()),
            ("text_pi_text", "<r>a<?pi x?>b</r>".into()),
            ("text_cdata_text", "<r>a<![CDATA[b]]>c</r>".into()),
            ("cdata_empty", "<r>a<![CDATA[]]>c</r>".into()),
            (
                "predefined_refs_only",
                "<r>&lt;&gt;&amp;&quot;&apos;</r>".into(),
            ),
            ("char_refs_mixed", "<r>x&#x41;y&#66;z</r>".into()),
            ("crlf_text", "<r>a\r\nb\rc</r>".into()),
            (
                "attr_ws",
                "<r x=\" a\tb\nc\r\nd \" y=\"&#9;&#10;&#13;\"/>".into(),
            ),
            ("big_text", format!("<r>{big}</r>")),
            ("big_text_refs", format!("<r>{big_refs}</r>")),
            ("mismatched_end_line2", "<r>\n<x></y></r>".into()),
            ("error_line3", "<r>\n<x>\n</y></r>".into()),
            ("two_roots", "<r/><s/>".into()),
            ("text_before_root", "x<r/>".into()),
            ("undeclared_entity", "<r>&foo;</r>".into()),
            ("duplicate_attr", "<r x=\"1\" x=\"2\"/>".into()),
            (
                "doctype_internal",
                "<!DOCTYPE r [<!ENTITY e \"EE\">]><r>&e;</r>".into(),
            ),
            ("nul_char_ref", "<r>&#0;</r>".into()),
            ("attr_lt", "<r x=\"a<b\"/>".into()),
            ("hex_upper_X", "<r>&#X41;</r>".into()),
            ("empty_doc", String::new()),
            ("unclosed", "<r>".into()),
            ("ws_after_root", "<r/>\n \t".into()),
        ];
        let mut documents: Vec<(&'static str, Vec<u8>)> = text
            .into_iter()
            .map(|(name, text)| (name, text.into_bytes()))
            .collect();
        documents.push((
            "latin1_declared",
            b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><r>\xe9</r>".to_vec(),
        ));
        documents.push(("bom_utf8", b"\xef\xbb\xbf<r>a</r>".to_vec()));
        documents
    }

    /// The oracle's escaping of text (control characters escaped).
    fn escape(text: &str) -> String {
        let mut out = String::new();
        for c in text.chars() {
            match c {
                '\n' => out.push_str("\\n"),
                '\r' => out.push_str("\\r"),
                '\t' => out.push_str("\\t"),
                '\\' => out.push_str("\\\\"),
                c if (c as u32) < 0x20 => out.push_str(&format!("\\x{:02x}", c as u32)),
                c => out.push(c),
            }
        }
        out
    }

    /// The oracle's rendering of a characters callback (UTF-16 length, long
    /// runs abbreviated).
    fn render_text(text: &str) -> String {
        let mut out = escape(text);
        let units = text.encode_utf16().count();
        if out.len() > 60 {
            out = format!("{}...<{units} units>", &out[..20]);
        }
        format!("  C len={units} \"{out}\"")
    }

    /// The oracle's records per case, with the locator positions (`@l:c`),
    /// which this driver does not expose, removed.
    fn release_records() -> Vec<(String, Vec<String>)> {
        let mut cases: Vec<(String, Vec<String>)> = Vec::new();
        for line in RELEASE.lines() {
            if let Some(name) = line.strip_prefix("CASE ") {
                cases.push((name.to_owned(), Vec::new()));
                continue;
            }
            let record = match line.rfind(" @") {
                Some(at) if !line.starts_with("  FATAL") => &line[..at],
                _ => line,
            };
            cases.last_mut().unwrap().1.push(record.to_owned());
        }
        cases
    }

    #[test]
    fn release_callbacks() {
        let documents = release_documents();
        let records = release_records();
        assert_eq!(records.len(), 28);
        assert_eq!(documents.len(), 28);
        for (name, expected) in records {
            let document = &documents
                .iter()
                .find(|(case, _)| *case == name)
                .unwrap_or_else(|| panic!("no document for {name}"))
                .1;
            // Xerces honours the ISO-8859-1 declaration; this driver does in
            // Declared mode (idxml) and refuses it in UTF-8 mode.
            let encoding = if name == "latin1_declared" {
                #[cfg(feature = "idxml")]
                {
                    SaxEncoding::Declared
                }
                #[cfg(not(feature = "idxml"))]
                {
                    SaxEncoding::Utf8
                }
            } else {
                SaxEncoding::Utf8
            };
            let mut got = Vec::new();
            let result = parse_with(document.as_slice(), encoding, &LIMITS, |event| {
                got.push(match event {
                    SaxEvent::Start(qname, attributes) => {
                        let rendered: Vec<String> = attributes
                            .0
                            .iter()
                            .map(|(key, value)| format!("{key}={}", escape(value)))
                            .collect();
                        format!("  S {qname} [{}]", rendered.join(", "))
                    }
                    SaxEvent::End(qname) => format!("  E {qname}"),
                    SaxEvent::Characters(text) => render_text(text),
                });
                Ok(())
            });
            let last = expected.last().unwrap().as_str();
            match last {
                "  OK" => {
                    if name == "doctype_internal"
                        || (name == "latin1_declared" && encoding == SaxEncoding::Utf8)
                    {
                        // Native refusal (module documentation).
                        assert!(
                            matches!(result, Err(Error::Unsupported(_))),
                            "{name}: {result:?}"
                        );
                        continue;
                    }
                    result.unwrap_or_else(|error| panic!("{name}: {error}"));
                    assert_eq!(got, expected[..expected.len() - 1], "{name}");
                }
                "  THROWN SAXParseException" => {
                    let fatal = &expected[expected.len() - 2];
                    let line: usize = fatal
                        .strip_prefix("  FATAL line=")
                        .and_then(|rest| rest.split(' ').next())
                        .and_then(|number| number.parse().ok())
                        .unwrap_or_else(|| panic!("{name}: {fatal}"));
                    match result {
                        Err(Error::Parse { line: ours, .. }) => {
                            assert_eq!(ours, line, "{name}: {fatal}")
                        }
                        other => panic!("{name}: expected a parse error, got {other:?}"),
                    }
                    // The callbacks before the fatal error match too.
                    assert_eq!(got, expected[..expected.len() - 2], "{name}");
                }
                other => panic!("{name}: unexpected oracle outcome {other}"),
            }
        }
    }
}

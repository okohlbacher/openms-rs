// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The reader of the legacy `.ttd` tool-description files: the native form of
//! core `FORMAT/ToolDescriptionFile.h` and `FORMAT/HANDLERS/ToolDescriptionHandler.h`
//! with their `.cpp` files (core `bc9cc12`).
//!
//! [`ToolHandler`](super::ToolHandler) reads these for its internal-tool
//! registry. The file format is a `<ttd>` element of `<tool>` elements: an
//! internal tool has a `<name>`, a `<category>` and `<type>` values; an
//! external one wraps a program with `<external>` details, file `<mappings>`
//! and an `<ini_param>` block of parameter-XML items.
//!
//! This lives with the command-line framework, which is its only consumer and
//! already has the XML reader (`paramxml`), rather than under `src/format`.

use crate::data_structures::{FileMapping, ToolDescription, ToolExternalDetails};
use crate::format::paramxml;
use crate::param::Param;
use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::Event;
use std::io::Read;
use std::path::Path;

/// Most bytes of one `.ttd` file.
pub const MAX_TTD_BYTES: u64 = 16 * 1024 * 1024;

/// Most XML elements in one `.ttd` file.
pub const MAX_TTD_ELEMENTS: usize = 100_000;

/// Deepest element nesting of one `.ttd` file.
pub const MAX_TTD_DEPTH: usize = 128;

/// The descriptions of one `.ttd` file and the non-fatal diagnostics its
/// reading raised.
#[derive(Clone, Debug, Default)]
pub struct LoadedToolDescriptions {
    /// One description per `<tool>` element, in document order.
    pub descriptions: Vec<ToolDescription>,
    /// The source handler's non-fatal `error(LOAD, …)` messages, each as the
    /// source prints it: `Non-fatal error while loading '<file>': <message>`.
    pub diagnostics: Vec<String>,
}

/// File adapter for tool-descriptor files (source `ToolDescriptionFile`).
///
/// Only [`load`](Self::load) is ported. The source's `store` throws
/// `Exception::NotImplemented` (its handler's `writeTo` is not implemented), so
/// there is nothing to port.
#[derive(Clone, Copy, Debug, Default)]
pub struct ToolDescriptionFile;

impl ToolDescriptionFile {
    /// The schema the source registers for the format.
    pub const SCHEMA_LOCATION: &'static str = "/SCHEMAS/ToolDescriptor_1_0.xsd";
    /// The schema version the source registers.
    pub const SCHEMA_VERSION: &'static str = "1.0.0";

    /// Load the tool descriptions of a `.ttd` file (source `load`).
    ///
    /// The element handling is `ToolDescriptionHandler`'s, including its
    /// leniency: an unknown element, an unknown `status` value of `<tool>` and
    /// text in an unexpected element are non-fatal diagnostics, and reading
    /// continues. Text of the one-value elements replaces the value, and each
    /// `<type>` adds one; text interrupted by an entity reference is joined.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read, and [`Error::Parse`] —
    /// `While loading '<file>': <reason>` — for XML that is not well formed,
    /// a `<tool>`, `<mapping>`, `<file_pre>` or `<file_post>` without a
    /// required attribute, a mapping `id` that is not an integer, a parameter
    /// block the parameter-XML reader refuses, or a file beyond
    /// [`MAX_TTD_BYTES`], [`MAX_TTD_ELEMENTS`] or [`MAX_TTD_DEPTH`]. The
    /// source's messages for malformed XML are Xerces's; these are
    /// quick-xml's.
    pub fn load(path: impl AsRef<Path>) -> Result<LoadedToolDescriptions> {
        let path = path.as_ref();
        let file = std::fs::File::open(path)?;
        let mut bytes = Vec::new();
        file.take(MAX_TTD_BYTES + 1).read_to_end(&mut bytes)?;
        let display = path.display().to_string();
        let fail = |reason: String| Error::Parse {
            line: 0,
            message: format!("While loading '{display}': {reason}"),
        };
        if bytes.len() as u64 > MAX_TTD_BYTES {
            return Err(fail(format!("the file exceeds {MAX_TTD_BYTES} bytes")));
        }
        let text = String::from_utf8(bytes).map_err(|_| fail("invalid UTF-8".into()))?;
        let mut handler = Handler::new(display.clone());
        let mut reader = Reader::from_str(&text);
        let mut depth = 0usize;
        let mut elements = 0usize;
        // Byte offset where the current <ini_param> content starts.
        let mut ini_start: Option<(usize, usize)> = None;
        loop {
            let before = usize::try_from(reader.buffer_position()).unwrap_or(usize::MAX);
            let event = reader.read_event().map_err(|e| fail(e.to_string()))?;
            match event {
                Event::Start(ref element) | Event::Empty(ref element) => {
                    let empty = matches!(event, Event::Empty(_));
                    elements += 1;
                    if elements > MAX_TTD_ELEMENTS || depth >= MAX_TTD_DEPTH {
                        return Err(fail(
                            "the tool description exceeds its element bounds".into(),
                        ));
                    }
                    let name = String::from_utf8_lossy(element.name().as_ref()).into_owned();
                    let mut attributes = Vec::new();
                    for attribute in element.attributes() {
                        let attribute = attribute.map_err(|e| fail(e.to_string()))?;
                        let key = String::from_utf8_lossy(attribute.key.as_ref()).into_owned();
                        // XML attribute-value normalisation, as the parameter
                        // reader applies it.
                        let raw = String::from_utf8_lossy(&attribute.value)
                            .replace(['\t', '\n', '\r'], " ");
                        let value = quick_xml::escape::unescape(&raw)
                            .map_err(|e| fail(e.to_string()))?
                            .into_owned();
                        attributes.push((key, value));
                    }
                    if ini_start.is_some() {
                        // Inside the parameter block: the parameter-XML reader
                        // sees the whole block once it ends.
                        if !empty {
                            depth += 1;
                        }
                        continue;
                    }
                    handler.start(&name, &attributes).map_err(fail)?;
                    if name == "ini_param" && handler.in_ini_section {
                        if empty {
                            handler.end("ini_param", Param::new());
                        } else {
                            let start =
                                usize::try_from(reader.buffer_position()).unwrap_or(usize::MAX);
                            ini_start = Some((start, depth));
                            depth += 1;
                        }
                        continue;
                    }
                    if empty {
                        handler.end(&name, Param::new());
                    } else {
                        depth += 1;
                    }
                }
                Event::End(ref element) => {
                    depth = depth.saturating_sub(1);
                    let name = String::from_utf8_lossy(element.name().as_ref()).into_owned();
                    if let Some((start, level)) = ini_start {
                        if name == "ini_param" && depth == level {
                            let inner = text.get(start..before).unwrap_or("");
                            let document = format!("<PARAMETERS>{inner}</PARAMETERS>");
                            let param = paramxml::read(document.as_bytes())
                                .map_err(|e| fail(e.to_string()))?;
                            ini_start = None;
                            handler.end("ini_param", param);
                        }
                        continue;
                    }
                    handler.end(&name, Param::new());
                }
                Event::Text(ref content) => {
                    if ini_start.is_none() {
                        let raw = content.xml_content().map_err(|e| fail(e.to_string()))?;
                        let value =
                            quick_xml::escape::unescape(&raw).map_err(|e| fail(e.to_string()))?;
                        handler.characters(&value);
                    }
                }
                Event::CData(ref content) => {
                    if ini_start.is_none() {
                        handler.characters(&String::from_utf8_lossy(content.as_ref()));
                    }
                }
                Event::GeneralRef(ref reference) => {
                    if ini_start.is_none() {
                        let name = reference.decode().map_err(|e| fail(e.to_string()))?;
                        let mut utf8 = [0u8; 4];
                        let value: &str = match reference
                            .resolve_char_ref()
                            .map_err(|e| fail(e.to_string()))?
                        {
                            Some(c) => c.encode_utf8(&mut utf8),
                            None => quick_xml::escape::resolve_xml_entity(&name)
                                .ok_or_else(|| fail(format!("unknown entity '{name}'")))?,
                        };
                        handler.characters(value);
                    }
                }
                Event::Eof => break,
                _ => {}
            }
        }
        if depth != 0 || ini_start.is_some() {
            return Err(fail(
                "input ended before all started tags were ended".into(),
            ));
        }
        handler.flush_text();
        Ok(LoadedToolDescriptions {
            descriptions: handler.descriptions,
            diagnostics: handler.diagnostics,
        })
    }
}

/// The source's `ToolDescriptionHandler` state outside the parameter block.
struct Handler {
    file: String,
    tool: ToolDescription,
    details: ToolExternalDetails,
    descriptions: Vec<ToolDescription>,
    diagnostics: Vec<String>,
    tag: String,
    open_tags: Vec<String>,
    in_ini_section: bool,
    /// Text of the current element, delivered when the element's text ends.
    text: String,
}

impl Handler {
    fn new(file: String) -> Self {
        Self {
            file,
            tool: ToolDescription::default(),
            details: ToolExternalDetails::default(),
            descriptions: Vec::new(),
            diagnostics: Vec::new(),
            tag: String::new(),
            open_tags: Vec::new(),
            in_ini_section: false,
            text: String::new(),
        }
    }

    /// Source `XMLHandler::error(LOAD, message)`.
    fn warn(&mut self, message: &str) {
        self.diagnostics.push(format!(
            "Non-fatal error while loading '{}': {message}",
            self.file
        ));
    }

    fn attribute<'a>(
        attributes: &'a [(String, String)],
        name: &str,
    ) -> std::result::Result<&'a str, String> {
        attributes
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
            .ok_or_else(|| format!("Required attribute '{name}' not present!"))
    }

    /// Source `onStartElement` (`ToolDescriptionHandler.cpp:30-112`).
    fn start(
        &mut self,
        name: &str,
        attributes: &[(String, String)],
    ) -> std::result::Result<(), String> {
        self.flush_text();
        self.tag = name.to_owned();
        self.open_tags.push(name.to_owned());
        match name {
            "tool" => {
                let status = Self::attribute(attributes, "status")?;
                match status {
                    "external" => self.tool.internal.is_internal = false,
                    "internal" => self.tool.internal.is_internal = true,
                    other => self.warn(&format!(
                        "ToolDescriptionHandler::startElement: Element 'status' if tag 'tool' has unknown value {other}'."
                    )),
                }
                return Ok(());
            }
            "mapping" => {
                let id = Self::attribute(attributes, "id")?;
                let id = super::context::to_int32(id)?;
                let command = Self::attribute(attributes, "cl")?.to_owned();
                self.details.tr_table.mapping.insert(id, command);
                return Ok(());
            }
            "file_post" | "file_pre" => {
                let mapping = FileMapping {
                    location: Self::attribute(attributes, "location")?.to_owned(),
                    target: Self::attribute(attributes, "target")?.to_owned(),
                };
                if name == "file_post" {
                    self.details.tr_table.post_moves.push(mapping);
                } else {
                    self.details.tr_table.pre_moves.push(mapping);
                }
                return Ok(());
            }
            "ini_param" => {
                self.in_ini_section = true;
                return Ok(());
            }
            "ttd" | "category" | "e_category" | "type" => return Ok(()),
            _ => {}
        }
        let allowed = if self.tool.internal.is_internal {
            name == "name"
        } else {
            matches!(
                name,
                "external"
                    | "cloptions"
                    | "path"
                    | "mappings"
                    | "mapping"
                    | "ini_param"
                    | "text"
                    | "onstartup"
                    | "onfail"
                    | "onfinish"
                    | "workingdirectory"
            )
        };
        if !allowed {
            self.warn(&format!(
                "ToolDescriptionHandler::startElement(): Unknown element found: '{name}', ignoring."
            ));
        }
        Ok(())
    }

    /// Collect text; it is delivered to the current tag by
    /// [`flush_text`](Self::flush_text).
    fn characters(&mut self, text: &str) {
        // Text outside the root element is not character data.
        if !self.open_tags.is_empty() {
            self.text.push_str(text);
        }
    }

    /// Source `onCharacters` (`ToolDescriptionHandler.cpp:114-172`) for the text
    /// collected since the last element boundary.
    fn flush_text(&mut self) {
        if self.text.is_empty() {
            return;
        }
        let text = std::mem::take(&mut self.text);
        match self.tag.as_str() {
            "ttd" | "tool" | "mappings" | "external" | "text" => {}
            "name" => self.tool.internal.name = text,
            "category" => self.tool.internal.category = text,
            "type" => self.tool.internal.types.push(text),
            "e_category" => self.details.category = text,
            "cloptions" => self.details.commandline = text,
            "path" => self.details.path = text,
            "onstartup" => self.details.text_startup = text,
            "onfail" => self.details.text_fail = text,
            "onfinish" => self.details.text_finish = text,
            "workingdirectory" => self.details.working_directory = text,
            other => {
                let message = format!(
                    "ToolDescriptionHandler::characters: Unknown character section found: '{other}', ignoring."
                );
                self.warn(&message);
            }
        }
    }

    /// Source `onEndElement` (`ToolDescriptionHandler.cpp:174-210`).
    fn end(&mut self, name: &str, param: Param) {
        self.flush_text();
        self.open_tags.pop();
        if let Some(last) = self.open_tags.last() {
            self.tag = last.clone();
        }
        match name {
            "ini_param" => {
                self.in_ini_section = false;
                self.details.param = param;
            }
            "external" => {
                self.tool
                    .external_details
                    .push(std::mem::take(&mut self.details));
            }
            "tool" => {
                self.descriptions.push(std::mem::take(&mut self.tool));
            }
            _ => {}
        }
    }
}

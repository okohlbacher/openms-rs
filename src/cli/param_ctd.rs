// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The Common Tool Description (CTD) writer: the native form of core
//! `FORMAT/ParamCTDFile.h` and `ParamCTDFile.cpp` (core `bc9cc12`).
//!
//! `TOPPBase` writes one CTD file per tool for `-write_ctd`, which workflow
//! systems such as KNIME and Galaxy consume. **The written bytes are the
//! product**, so this writer reproduces the source's text exactly, including
//! three defects a consumer can observe:
//!
//! 1. The source's character replacement skips the character after each
//!    replacement, so of two adjacent special characters only the first is
//!    escaped: `&&` is written `&amp;&`, and a second consecutive line break is
//!    left raw instead of becoming `#br#`.
//! 2. A tab in a string value becomes `&#x9;` before the value is escaped, so
//!    it is written `&amp;#x9;`, which reads back as the six characters
//!    `&#x9;`.
//! 3. The `<tool>` attributes and the description are not escaped at all.
//!
//! Each is logged in `OpenMS_CPP_ISSUES.md` (CPP-349). The writer is
//! byte-identical to the Release build for the eight ported tools and to the
//! retained class-test file `ParamCTDFile_test_writeCTDToStream.ctd`; see
//! `docs/TOPP_CLI_SUPPORT.md` (*Tool descriptions*).
//!
//! It lives with the command-line framework, its only consumer, rather than
//! under `src/format`.

use crate::data_structures::ToolInfo;
use crate::format::file_info::text_format::ostream_g;
use crate::param::{Param, ParamValue};
use crate::{Error, Result};
use std::io::Write;
use std::path::Path;

/// Significant digits of a list item: the source stream's precision,
/// `std::numeric_limits<double>::digits10`.
const LIST_ITEM_DIGITS: u32 = 15;

/// Most bytes one CTD document may take, checked before anything is written.
pub const MAX_CTD_BYTES: usize = 64 * 1024 * 1024;

/// The CTD writer (source class `ParamCTDFile`).
#[derive(Clone, Copy, Debug, Default)]
pub struct ParamCtdFile;

impl ParamCtdFile {
    /// The parameter-schema version written as `<PARAMETERS version=…>`
    /// (source `schema_version_`).
    pub const SCHEMA_VERSION: &'static str = "1.8.0";
    /// The schema path appended to the OpenMS share URL (source
    /// `schema_location_`).
    pub const SCHEMA_LOCATION: &'static str = "/SCHEMAS/Param_1_8_0.xsd";

    /// Write `param` and `tool_info` as a CTD file (source `store`).
    ///
    /// The source writes to standard output for the file name `-`; so does
    /// this.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] with the source's text, `Unable to create file: <path>`,
    /// when the file cannot be created, where the source throws
    /// `std::ios::failure` with that message; other failures as
    /// [`write_ctd_to_stream`](Self::write_ctd_to_stream). The document is
    /// built completely before the file is opened, so a document that cannot
    /// be built leaves no file behind.
    pub fn store(&self, path: impl AsRef<Path>, param: &Param, tool_info: &ToolInfo) -> Result<()> {
        let path = path.as_ref();
        let text = self.to_ctd_string(param, tool_info)?;
        if path == Path::new("-") {
            let mut out = std::io::stdout().lock();
            out.write_all(text.as_bytes())?;
            out.flush()?;
            return Ok(());
        }
        let mut file = std::fs::File::create(path).map_err(|error| {
            Error::Io(std::io::Error::new(
                error.kind(),
                format!("Unable to create file: {}", path.display()),
            ))
        })?;
        file.write_all(text.as_bytes())?;
        file.flush()?;
        Ok(())
    }

    /// Write the CTD document to a stream (source `writeCTDToStream`).
    ///
    /// # Errors
    ///
    /// As [`to_ctd_string`](Self::to_ctd_string), and [`Error::Io`] when the
    /// stream fails.
    pub fn write_ctd_to_stream(
        &self,
        out: &mut dyn Write,
        param: &Param,
        tool_info: &ToolInfo,
    ) -> Result<()> {
        out.write_all(self.to_ctd_string(param, tool_info)?.as_bytes())?;
        out.flush()?;
        Ok(())
    }

    /// The CTD document as the source writes it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the document would exceed
    /// [`MAX_CTD_BYTES`], and parameter-tree failures of
    /// [`Param::iter`].
    pub fn to_ctd_string(&self, param: &Param, tool_info: &ToolInfo) -> Result<String> {
        let mut os = Bounded::default();
        os.add("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
        os.add(&format!(
            "<tool ctdVersion=\"1.8\" version=\"{}\" name=\"{}\" docurl=\"{}\" category=\"{}\" >\n",
            tool_info.version, tool_info.name, tool_info.docurl, tool_info.category
        ))?;
        os.add(&format!(
            "<description><![CDATA[{}]]></description>\n",
            tool_info.description
        ))?;
        os.add(&format!(
            "<manual><![CDATA[{}]]></manual>\n",
            tool_info.description
        ))?;
        os.add("<citations>\n")?;
        for doi in &tool_info.citations {
            os.add(&format!("  <citation doi=\"{doi}\" url=\"\" />\n"))?;
        }
        os.add("</citations>\n")?;
        os.add(&format!(
            "<PARAMETERS version=\"{}\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS{}\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
            Self::SCHEMA_VERSION,
            Self::SCHEMA_LOCATION
        ))?;

        let mut indentation = 2usize;
        let items = param.iter()?;
        let end_trace = items.end_trace().to_vec();
        let mut any = false;
        for item in items {
            any = true;
            for trace in &item.trace {
                if trace.opened {
                    let description = replace(&trace.description, b'\n', "#br#");
                    os.add(&format!(
                        "{}<NODE name=\"{}\" description=\"{}\">\n",
                        " ".repeat(indentation),
                        escape_xml(&trace.name),
                        escape_xml(&description)
                    ))?;
                    indentation += 2;
                } else {
                    indentation = indentation.saturating_sub(2);
                    os.add(&format!("{}</NODE>\n", " ".repeat(indentation)))?;
                }
            }
            write_entry(&mut os, item.entry, indentation)?;
        }
        if any {
            for _ in &end_trace {
                indentation = indentation.saturating_sub(2);
                os.add(&format!("{}</NODE>\n", " ".repeat(indentation)))?;
            }
        }
        os.add("</PARAMETERS>\n")?;
        os.add("</tool>\n")?;
        Ok(os.text)
    }
}

/// A document under construction, bounded by [`MAX_CTD_BYTES`].
#[derive(Default)]
struct Bounded {
    text: String,
}

impl Bounded {
    fn add(&mut self, piece: &str) -> Result<()> {
        if self.text.len().saturating_add(piece.len()) > MAX_CTD_BYTES {
            return Err(Error::InvalidValue(format!(
                "the CTD document exceeds {MAX_CTD_BYTES} bytes"
            )));
        }
        self.text.push_str(piece);
        Ok(())
    }
}

/// One `<ITEM>` or `<ITEMLIST>` (`ParamCTDFile.cpp:89-325`).
fn write_entry(
    os: &mut Bounded,
    entry: &crate::param::ParamEntry,
    indentation: usize,
) -> Result<()> {
    let value = &entry.value;
    if matches!(value, ParamValue::Empty) {
        return Ok(());
    }
    let pad = " ".repeat(indentation);
    let scalar = matches!(
        value,
        ParamValue::String(_) | ParamValue::Integer(_) | ParamValue::Float(_)
    );
    let mut tags = entry.tags.clone();
    let mut is_flag = false;
    let mut line = String::new();
    if scalar {
        line.push_str(&format!(
            "{pad}<ITEM name=\"{}\" value=\"",
            escape_xml(&entry.name)
        ));
    } else {
        line.push_str(&format!(
            "{pad}<ITEMLIST name=\"{}",
            escape_xml(&entry.name)
        ));
    }
    let text = |v: &ParamValue| v.to_text(true);
    match value {
        ParamValue::Integer(_) => {
            line.push_str(&format!("{}\" type=\"int\"", text(value)?));
        }
        ParamValue::Float(_) => {
            line.push_str(&format!("{}\" type=\"double\"", text(value)?));
        }
        ParamValue::String(string) => {
            if tags.remove("input file") {
                line.push_str(&format!("{}\" type=\"input-file\"", escape_xml(string)));
            } else if tags.remove("output file") {
                line.push_str(&format!("{}\" type=\"output-file\"", escape_xml(string)));
            } else if tags.remove("output dir") {
                line.push_str(&format!("{}\" type=\"output-dir\"", escape_xml(string)));
            } else if tags.remove("output prefix") {
                line.push_str(&format!("{}\" type=\"output-prefix\"", escape_xml(string)));
            } else if entry.valid_strings.len() == 2
                && entry.valid_strings[0] == "true"
                && entry.valid_strings[1] == "false"
                && string == "false"
            {
                is_flag = true;
                line.push_str(&format!("{string}\" type=\"bool\""));
            } else {
                let tabbed = if string.contains('\t') {
                    replace(string, b'\t', "&#x9;")
                } else {
                    string.clone()
                };
                line.push_str(&format!("{}\" type=\"string\"", escape_xml(&tabbed)));
            }
        }
        ParamValue::StringList(_) => {
            if tags.remove("input file") {
                line.push_str("\" type=\"input-file\"");
            } else if tags.remove("output file") {
                line.push_str("\" type=\"output-file\"");
            } else {
                line.push_str("\" type=\"string\"");
            }
        }
        ParamValue::IntegerList(_) => line.push_str("\" type=\"int\""),
        ParamValue::FloatList(_) => line.push_str("\" type=\"double\""),
        ParamValue::Empty => {}
    }

    let description = replace(&entry.description, b'\n', "#br#");
    line.push_str(&format!(" description=\"{}\"", escape_xml(&description)));
    line.push_str(if tags.remove("required") {
        " required=\"true\""
    } else {
        " required=\"false\""
    });
    line.push_str(if tags.remove("advanced") {
        " advanced=\"true\""
    } else {
        " advanced=\"false\""
    });
    if !tags.is_empty() {
        let list: Vec<&str> = tags.iter().map(String::as_str).collect();
        line.push_str(&format!(" tags=\"{}\"", escape_xml(&list.join(","))));
    }

    if !is_flag {
        let mut restrictions = String::new();
        match value {
            ParamValue::Integer(_) | ParamValue::IntegerList(_) => {
                let min_set = entry.min_int != -i32::MAX;
                let max_set = entry.max_int != i32::MAX;
                if min_set || max_set {
                    if min_set {
                        restrictions.push_str(&entry.min_int.to_string());
                    }
                    restrictions.push(':');
                    if max_set {
                        restrictions.push_str(&entry.max_int.to_string());
                    }
                }
            }
            ParamValue::Float(_) | ParamValue::FloatList(_) => {
                let min_set = entry.min_float != -f64::MAX;
                let max_set = entry.max_float != f64::MAX;
                if min_set || max_set {
                    if min_set {
                        restrictions.push_str(&to_string_double(entry.min_float));
                    }
                    restrictions.push(':');
                    if max_set {
                        restrictions.push_str(&to_string_double(entry.max_float));
                    }
                }
            }
            ParamValue::String(_) | ParamValue::StringList(_) => {
                restrictions = entry.valid_strings.join(",");
            }
            ParamValue::Empty => {}
        }
        if !restrictions.is_empty() {
            let file_like = entry.tags.contains("input file")
                || entry.tags.contains("output file")
                || entry.tags.contains("output prefix");
            let attribute = if file_like {
                "supported_formats"
            } else {
                "restrictions"
            };
            line.push_str(&format!(" {attribute}=\"{}\"", escape_xml(&restrictions)));
        }
    }
    line.push_str(if scalar { " />\n" } else { " >\n" });
    os.add(&line)?;

    let inner = " ".repeat(indentation + 2);
    match value {
        ParamValue::StringList(items) => {
            for item in items {
                let tabbed = if item.contains('\t') {
                    replace(item, b'\t', "&#x9;")
                } else {
                    item.clone()
                };
                os.add(&format!(
                    "{inner}<LISTITEM value=\"{}\"/>\n",
                    escape_xml(&tabbed)
                ))?;
            }
        }
        ParamValue::IntegerList(items) => {
            for item in items {
                os.add(&format!("{inner}<LISTITEM value=\"{item}\"/>\n"))?;
            }
        }
        ParamValue::FloatList(items) => {
            for item in items {
                os.add(&format!(
                    "{inner}<LISTITEM value=\"{}\"/>\n",
                    ostream_g(*item, LIST_ITEM_DIGITS)
                ))?;
            }
        }
        _ => {}
    }
    if !scalar {
        os.add(&format!("{pad}</ITEMLIST>\n"))?;
    }
    Ok(())
}

/// `std::to_string(double)`: `printf("%f")`, with glibc's `nan`/`-nan`.
fn to_string_double(value: f64) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    format!("{value:.6}")
}

impl ParamCtdFile {
    /// Source `ParamCTDFile::escapeXML` (a private static there):
    /// [`escape_xml`] as an associated function.
    pub fn escape_xml(text: &str) -> String {
        escape_xml(text)
    }
    /// Source `ParamCTDFile::replace` (a private static there): [`replace`]
    /// as an associated function.
    pub fn replace(text: &str, from: u8, to: &str) -> String {
        replace(text, from, to)
    }
}

/// Source `ParamCTDFile::escapeXML`: `&`, `>`, `"`, `<` and `'` in that
/// order, each through [`replace`], so the source's skip after a replacement
/// applies to each pass.
pub fn escape_xml(text: &str) -> String {
    let mut copy = text.to_owned();
    for (from, to) in [
        (b'&', "&amp;"),
        (b'>', "&gt;"),
        (b'"', "&quot;"),
        (b'<', "&lt;"),
        (b'\'', "&apos;"),
    ] {
        if copy.as_bytes().contains(&from) {
            copy = replace(&copy, from, to);
        }
    }
    copy
}

/// Source `ParamCTDFile::replace` (`ParamCTDFile.cpp:358-368`): replace each
/// byte `from` with `to`, except that after a replacement the next byte is
/// skipped — the loop index is advanced past the inserted text and then
/// incremented once more. Two adjacent `from` bytes therefore leave the second
/// one in place.
///
/// `from` is always an ASCII byte here, so a replacement never splits a UTF-8
/// sequence and the skipped byte never starts one mid-character.
pub fn replace(text: &str, from: u8, to: &str) -> String {
    // Linear form of the source loop: a replaced byte is followed by the next
    // input byte copied unchecked.
    let input = text.as_bytes();
    let mut out = Vec::with_capacity(input.len());
    let mut j = 0usize;
    while j < input.len() {
        if input[j] == from {
            out.extend_from_slice(to.as_bytes());
            if let Some(next) = input.get(j + 1) {
                out.push(*next);
            }
            j += 2;
        } else {
            out.push(input[j]);
            j += 1;
        }
    }
    // Only an ASCII byte was exchanged for ASCII text, so this is valid UTF-8.
    String::from_utf8(out).unwrap_or_default()
}

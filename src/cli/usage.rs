// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Usage text, the native form of the source `TOPPBase::printUsage_`
//! (`TOPPBase.cpp:631-888`).
//!
//! The layout is the one the source writes when standard error is not a
//! terminal and `COLUMNS` is unset, which is how the product SDK oracle ran
//! (`../oracle/topp-cli-lifecycle/cli2/manifest.json`, cases `help_<tool>` and
//! `helphelp_<tool>`): no colours and no line shaping to a console width.
//! See `docs/TOPP_CLI_SUPPORT.md` for the remaining differences.

use super::parameter::{ParameterInformation, ParameterType};
use super::spec::ToolSpec;
use crate::data_structures::list::ListFormat;
use crate::param::{Param, ParamValue};
use std::collections::BTreeMap;
use std::io::{Result as IoResult, Write};

/// `TOPPBase::cite_openms`, as its `Citation::toString` renders it.
const CITE_OPENMS: &str = "Pfeuffer, J., Bielow, C., Wein, S. et al.. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. Nat Methods (2024). doi:10.1038/s41592-024-02197-7.";

/// Most lines one written item keeps, as the source `IndentedStream(cerr, 0, 10)`.
const MAX_LINES: usize = 10;

/// Source `IndentedStream`, with an unlimited line width.
///
/// Each item is split at its line breaks; every line after the first is
/// indented by the current indentation, an item ending in a line break leaves
/// an indented empty line, and an item of more than [`MAX_LINES`] lines keeps
/// its first eight lines, an indented `...` and its last line
/// (`ConsoleUtils::breakString_`).
struct Indented<'a> {
    out: &'a mut dyn Write,
    indentation: usize,
}

impl Indented<'_> {
    fn put(&mut self, item: &str) -> IoResult<()> {
        if item.is_empty() {
            return Ok(());
        }
        let prefix = " ".repeat(self.indentation);
        // A trailing line break yields an empty last part, which becomes the
        // indented empty line the source appends.
        let mut lines: Vec<String> = item
            .split('\n')
            .enumerate()
            .map(|(index, line)| {
                if index == 0 {
                    line.to_owned()
                } else {
                    format!("{prefix}{line}")
                }
            })
            .collect();
        if lines.len() > MAX_LINES {
            let last = lines.pop().unwrap_or_default();
            lines.truncate(MAX_LINES - 2);
            lines.push(format!("{prefix}..."));
            lines.push(last);
        }
        let mut first = true;
        for line in &lines {
            if !first {
                self.out.write_all(b"\n")?;
            }
            first = false;
            self.out.write_all(line.as_bytes())?;
        }
        Ok(())
    }
}

/// Source `StringUtils::fillRight`: pad with spaces to `width` bytes.
fn fill_right(text: &mut String, width: usize) {
    while text.len() < width {
        text.push(' ');
    }
}

/// Source `StringUtils::firstToUpper`: the first byte, if an ASCII letter,
/// becomes upper case.
fn first_to_upper(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(first) if first.is_ascii() => {
            let mut upper = String::with_capacity(text.len());
            upper.push(first.to_ascii_uppercase());
            upper.push_str(chars.as_str());
            upper
        }
        _ => text.to_owned(),
    }
}

/// Source `getSubsection_`: the name before its last `:`, or empty.
fn subsection_of(name: &str) -> &str {
    name.rfind(':').map_or("", |position| &name[..position])
}

/// Source `ParamValue::toString()` at full precision.
fn value_text(value: &ParamValue) -> String {
    value.to_text(true).unwrap_or_default()
}

/// The `(default: …)` addon of one parameter, as the source lists it.
fn default_addon(entry: &ParameterInformation) -> Option<String> {
    match entry.kind {
        ParameterType::String
        | ParameterType::Double
        | ParameterType::Int
        | ParameterType::StringList
        | ParameterType::IntList
        | ParameterType::DoubleList => {
            let text = value_text(&entry.default_value).replace(", ", " ");
            (!text.is_empty() && text != "[]").then(|| format!("default: '{text}'"))
        }
        _ => None,
    }
}

/// Source `StringUtils::toStr(double)`, as `std::string + double` appends it.
fn float_text(value: f64) -> String {
    value
        .to_list_text()
        .map(|text| text.into_owned())
        .unwrap_or_default()
}

/// The restriction addons of one parameter, as the source lists them.
fn restrictions(entry: &ParameterInformation) -> Vec<String> {
    let mut restrictions = Vec::new();
    match entry.kind {
        ParameterType::String | ParameterType::StringList => {
            if !entry.valid_strings.is_empty() {
                let quoted: Vec<String> = entry
                    .valid_strings
                    .iter()
                    .map(|value| format!("'{value}'"))
                    .collect();
                restrictions.push(format!("valid: {}", quoted.join(", ")));
            }
        }
        ParameterType::InputFile
        | ParameterType::OutputFile
        | ParameterType::OutputPrefix
        | ParameterType::OutputDir
        | ParameterType::InputFileList
        | ParameterType::OutputFileList => {
            let quoted: Vec<String> = entry
                .accepted_formats()
                .map(|value| format!("'{value}'"))
                .collect();
            if !quoted.is_empty() {
                restrictions.push(format!("valid formats: {}", quoted.join(", ")));
            }
        }
        ParameterType::Int | ParameterType::IntList => {
            if let Some(min) = entry.min_int {
                restrictions.push(format!("min: '{min}'"));
            }
            if let Some(max) = entry.max_int {
                restrictions.push(format!("max: '{max}'"));
            }
        }
        ParameterType::Double | ParameterType::DoubleList => {
            if let Some(min) = entry.min_float {
                restrictions.push(format!("min: '{}'", float_text(min)));
            }
            if let Some(max) = entry.max_float {
                restrictions.push(format!("max: '{}'", float_text(max)));
            }
        }
        _ => {}
    }
    restrictions
}

/// The documentation URL of a tool, as `getDocumentationURL` builds it for a
/// release core version.
pub(crate) fn documentation_url(name: &str) -> String {
    format!(
        "http://www.openms.de/doxygen/release/{}/html/TOPP_{name}.html",
        crate::CORE_SDK_VERSION
    )
}

/// Print the usage block, as `printUsage_`.
///
/// `version` is the verbose version line, `subsection_defaults` the tool's
/// subsection parameters as `getSubsectionDefaults_` returns them, and
/// `verbose` the `--helphelp` request, which lists advanced parameters and
/// every subsection parameter instead of the subsection summary.
pub(crate) fn print(
    out: &mut dyn Write,
    name: &str,
    description: &str,
    version: &str,
    spec: &ToolSpec,
    subsection_defaults: &Param,
    verbose: bool,
) -> IoResult<()> {
    let url = documentation_url(name);
    let mut stream = Indented {
        out,
        indentation: 0,
    };
    stream.put(&format!(
        "\n{name} -- {description}\nFull documentation: {url}\nVersion: {version}\nTo cite OpenMS:\n + "
    ))?;
    stream.indentation = 3;
    stream.put(CITE_OPENMS)?;
    stream.indentation = 0;
    stream.put("\n\nUsage:\n")?;
    stream.put(&format!("  {name} <options>\n\n"))?;

    if !spec.subsections().is_empty() && !verbose {
        stream.put("This tool has algorithm parameters that are not shown here! Please check the ini file for a detailed description or use the --helphelp option\n\n")?;
    }

    // Under --helphelp the subsection parameters join the registered ones, as
    // the source's registerFullParam_; a failure leaves what was registered.
    let mut full = spec.clone();
    if verbose {
        let _ = full.register_full_param(subsection_defaults);
    }
    let topp_subsections: BTreeMap<&str, &str> = full
        .topp_subsections()
        .iter()
        .map(|(section, text)| (section.as_str(), text.as_str()))
        .collect();

    stream.put("Options (mandatory options marked with '*'):\n")?;
    let shown = |entry: &&ParameterInformation| !entry.advanced || verbose;
    let widest = full
        .parameters()
        .iter()
        .filter(shown)
        .map(|entry| entry.name.len() + entry.argument.len() + usize::from(entry.required))
        .max()
        .unwrap_or(0);
    let offset = 6 + widest;

    let mut current_subsection = "";
    for entry in full.parameters().iter().filter(shown) {
        let subsection = subsection_of(&entry.name);
        if !subsection.is_empty() && current_subsection != subsection {
            current_subsection = subsection;
            let text = topp_subsections
                .get(subsection)
                .copied()
                .filter(|text| !text.is_empty())
                .unwrap_or(subsection);
            stream.put("\n")?;
            stream.put(&format!("{text}:\n"))?;
        } else if subsection.is_empty() && !current_subsection.is_empty() {
            current_subsection = "";
            stream.put("\n")?;
        }

        let mut left = format!("  -{} {}", entry.name, entry.argument);
        if entry.required {
            left.push('*');
        }
        if entry.kind == ParameterType::Newline {
            left.clear();
        }
        fill_right(&mut left, offset);
        let text = first_to_upper(&entry.description);
        if entry.kind == ParameterType::Text {
            stream.put(&text)?;
        } else {
            let addons: Vec<String> = default_addon(entry).into_iter().collect();
            let restrictions = restrictions(entry);
            stream.indentation = offset;
            stream.put(&left)?;
            stream.put(&text)?;
            if !addons.is_empty() {
                stream.put(&format!(" ({})", addons.join(" ")))?;
            }
            if !restrictions.is_empty() {
                stream.put(&format!(" ({})", restrictions.join(" ")))?;
            }
            stream.indentation = 0;
        }
        stream.put("\n")?;
    }

    if !spec.subsections().is_empty() && !verbose {
        let sections: BTreeMap<&str, &str> = spec
            .subsections()
            .iter()
            .map(|(section, text)| (section.as_str(), text.as_str()))
            .collect();
        let indent = sections
            .keys()
            .map(|section| section.len())
            .max()
            .unwrap_or(0)
            + 6;
        stream.put("\nThe following configuration subsections are valid:\n")?;
        for (section, text) in sections {
            let mut line = format!(" - {section}");
            fill_right(&mut line, indent);
            stream.indentation = indent;
            stream.put(&format!("{line}{text}"))?;
            stream.indentation = 0;
            stream.put("\n")?;
        }
        stream.put(&format!(
            "\nYou can write an example INI file using the '-write_ini' option.\nDocumentation of subsection parameters can be found in the doxygen documentation or the INIFileEditor.\nFor more information, please consult the online documentation for this tool:\n  - {url}\n"
        ))?;
    }
    stream.put("\n")?;
    stream.out.flush()
}

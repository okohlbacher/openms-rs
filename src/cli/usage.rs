// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Usage text, the native form of the source `TOPPBase::printUsage_`
//! (`TOPPBase.cpp:631-888`).
//!
//! The text is written item by item as the source inserts it into its
//! `IndentedStream`, so that a console width shapes it as the source's does
//! ([`Console`]); written through explicit streams it is the layout the
//! source writes when standard error is not a terminal and `COLUMNS` is
//! unset: no colours and no line shaping. See `docs/TOPP_CLI_SUPPORT.md`.

use super::console::{Console, break_string};
use super::defs::{CITE_OPENMS, Citation};
use super::parameter::{ParameterInformation, ParameterType};
use super::spec::ToolSpec;
use crate::concept::log_stream::LogColor;
use crate::data_structures::list::ListFormat;
use crate::param::{Param, ParamValue};
use std::collections::BTreeMap;
use std::io::{Result as IoResult, Write};

/// Most lines one written item keeps, as the source `IndentedStream(cerr, 0, 10)`.
const MAX_LINES: usize = 10;

/// Source `IndentedStream` (`IndentedStream.h:57-88`, `IndentedStream.cpp`)
/// with the source's `Colorizer` insertions.
///
/// Each item is broken by `ConsoleUtils::breakStringList` at the current
/// indentation, starting at the current column; the column afterwards is the
/// length of the last line written. A coloured item is broken the same way
/// and wrapped in its colour's codes when the console is coloured, also when
/// it is empty.
struct Indented<'a> {
    out: &'a mut dyn Write,
    indentation: usize,
    column: usize,
    console: Console,
}

impl Indented<'_> {
    /// Break `item` and advance the column, returning the bytes to write.
    fn shape(&mut self, item: &[u8]) -> Vec<u8> {
        let lines = break_string(
            item,
            self.indentation,
            MAX_LINES,
            self.column,
            self.console.width,
        );
        let Some(last) = lines.last() else {
            return Vec::new();
        };
        if lines.len() == 1 {
            self.column += last.len();
        } else {
            self.column = last.len();
        }
        lines.join(&b'\n')
    }

    /// `is << item`.
    fn put(&mut self, item: &str) -> IoResult<()> {
        self.put_bytes(item.as_bytes())
    }

    /// `is << item` for an item that may hold a partial UTF-8 sequence, as a
    /// line the source broke inside one.
    fn put_bytes(&mut self, item: &[u8]) -> IoResult<()> {
        let bytes = self.shape(item);
        self.out.write_all(&bytes)
    }

    /// `is << colour(item)`: the item in a colour, then back.
    fn coloured(&mut self, colour: LogColor, item: &str) -> IoResult<()> {
        let bytes = self.shape(item.as_bytes());
        self.start(colour)?;
        self.out.write_all(&bytes)?;
        self.end(colour)
    }

    /// `is << colour()`: switch a colour on.
    fn start(&mut self, colour: LogColor) -> IoResult<()> {
        if self.console.colour {
            self.out.write_all(colour.enable())?;
        }
        Ok(())
    }

    /// `is << colour().undo()`: switch it off.
    fn end(&mut self, colour: LogColor) -> IoResult<()> {
        if self.console.colour {
            self.out.write_all(colour.disable())?;
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

/// Print the usage block, as `printUsage_`, laid out for `console`.
///
/// `version` is the verbose version line, `citations` the tool's own
/// citations, printed as `To cite <tool>:` after the OpenMS citation,
/// `subsection_defaults` the tool's subsection parameters as
/// `getSubsectionDefaults_` returns them, and `verbose` the `--helphelp`
/// request, which lists advanced parameters and every subsection parameter
/// instead of the subsection summary.
#[allow(clippy::too_many_arguments)]
pub(crate) fn print(
    out: &mut dyn Write,
    console: Console,
    name: &str,
    description: &str,
    version: &str,
    citations: &[Citation],
    spec: &ToolSpec,
    subsection_defaults: &Param,
    verbose: bool,
) -> IoResult<()> {
    let url = documentation_url(name);
    let mut stream = Indented {
        out,
        indentation: 0,
        column: 0,
        console,
    };
    // TOPPBase.cpp:638-645, one insertion at a time.
    stream.put("\n")?;
    stream.coloured(LogColor::Invert, name)?;
    stream.put(" -- ")?;
    stream.put(description)?;
    stream.put("\n")?;
    stream.coloured(LogColor::Bright, "Full documentation: ")?;
    stream.coloured(LogColor::Underline, &url)?;
    stream.put("\n")?;
    stream.coloured(LogColor::Bright, "Version: ")?;
    stream.put(version)?;
    stream.put("\n")?;
    stream.coloured(LogColor::Bright, "To cite OpenMS:\n")?;
    stream.put(" + ")?;
    stream.indentation = 3;
    stream.put(&CITE_OPENMS.to_source_string())?;
    stream.indentation = 0;
    stream.put("\n")?;
    // The tool's own citations (TOPPBase.cpp:646-651).
    if !citations.is_empty() {
        stream.start(LogColor::Bright)?;
        stream.put("To cite ")?;
        stream.put(name)?;
        stream.put(":")?;
        stream.end(LogColor::Bright)?;
        stream.indentation = 0;
        stream.put("\n")?;
        for citation in citations {
            stream.put(" + ")?;
            stream.indentation = 3;
            stream.put(&citation.to_source_string())?;
            stream.indentation = 0;
            stream.put("\n")?;
        }
    }
    stream.indentation = 0;
    stream.put("\n")?;
    stream.coloured(LogColor::Invert, "Usage:")?;
    stream.put("\n")?;
    stream.put("  ")?;
    stream.coloured(LogColor::Bright, name)?;
    stream.put(" <options>")?;
    stream.put("\n")?;
    stream.put("\n")?;

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

    stream.coloured(LogColor::Bright, "Options")?;
    stream.put(" (")?;
    stream.coloured(LogColor::Green, "mandatory options marked with '*'")?;
    stream.put("):\n")?;
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
            stream.put(text)?;
            stream.put(":\n")?;
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
            // The source's name column is empty for a text entry.
            stream.put(&text)?;
        } else {
            let addons: Vec<String> = default_addon(entry).into_iter().collect();
            let restrictions = restrictions(entry);
            let addons = if addons.is_empty() {
                String::new()
            } else {
                format!(" ({})", addons.join(" "))
            };
            let restrictions = if restrictions.is_empty() {
                String::new()
            } else {
                format!(" ({})", restrictions.join(" "))
            };
            stream.indentation = offset;
            if entry.required {
                stream.coloured(LogColor::Green, &left)?;
            } else {
                stream.put(&left)?;
            }
            stream.put(&text)?;
            stream.coloured(LogColor::Cyan, &addons)?;
            stream.coloured(LogColor::Magenta, &restrictions)?;
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
        stream.put("\n")?;
        stream.put("The following configuration subsections are valid:\n")?;
        for (section, text) in sections {
            let mut line = format!(" - {section}");
            fill_right(&mut line, indent);
            line.push_str(text);
            // `ConsoleUtils::breakString(tmp + description, indent, 10)`,
            // which the stream then breaks again at its own indentation.
            let broken =
                break_string(line.as_bytes(), indent, MAX_LINES, 0, console.width).join(&b'\n');
            stream.put_bytes(&broken)?;
            stream.put("\n")?;
        }
        stream.put("\n")?;
        stream.put("You can write an example INI file using the '-write_ini' option.\n")?;
        stream.put("Documentation of subsection parameters can be found in the doxygen documentation or the INIFileEditor.\n")?;
        stream.put(
            "For more information, please consult the online documentation for this tool:\n",
        )?;
        stream.put("  - ")?;
        stream.coloured(LogColor::Underline, &url)?;
        stream.put("\n")?;
    }
    // `is << endl`, which bypasses the shaping.
    stream.out.write_all(b"\n")?;
    stream.out.flush()
}

#[cfg(test)]
mod tests {
    //! The usage text on a terminal against the Release build
    //! (`../oracle/toppbase-completion/console.sh`, cases `tty_*`, retained in
    //! `tests/data/topp_cli_console`). A test process has no terminal, so the
    //! layout the executable would pick is passed in: `stty cols 70` gives
    //! width 69, `COLUMNS=50` width 49, a 0x0 terminal no shaping. That the
    //! executable picks it is compared in a pseudo-terminal by
    //! `../oracle/toppbase-completion/compare_tty.py`.
    use super::*;
    use crate::cli::console::UNSHAPED;
    use crate::cli::tools::BaselineFilter;
    use crate::cli::{
        Tool, ToolHandler, ToolRegistrySources, product_versions, subsection_defaults, tool_spec,
    };

    fn oracle(case: &str, file: &str) -> Vec<u8> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/topp_cli_console")
            .join(case)
            .join(file);
        std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
    }

    /// What reached the terminal, less the reset the Release build writes to
    /// it at exit, which the executable writes after the usage text.
    fn terminal(case: &str) -> Vec<u8> {
        let bytes = oracle(case, "tty.bin");
        bytes
            .strip_suffix(b"\x1b[0m")
            .unwrap_or_else(|| panic!("{case}: the Release build resets the terminal at exit"))
            .to_vec()
    }

    fn usage<T: Tool>(console: Console, verbose: bool) -> Vec<u8> {
        let spec = tool_spec::<T>().unwrap();
        let subsections = subsection_defaults::<T>(&spec).unwrap();
        let mut sources = ToolRegistrySources::only_prefixes(Vec::new()).unwrap();
        sources.builtin_manifest = true;
        let (_, version) = product_versions(&ToolHandler::new(sources), T::NAME).unwrap();
        let mut out = Vec::new();
        print(
            &mut out,
            console,
            T::NAME,
            T::DESCRIPTION,
            &version,
            T::CITATIONS,
            &spec,
            &subsections,
            verbose,
        )
        .unwrap();
        out
    }

    fn coloured(width: i64) -> Console {
        Console {
            width,
            colour: true,
        }
    }

    #[test]
    fn usage_on_a_terminal_matches_the_release_build() {
        assert_eq!(
            usage::<BaselineFilter>(coloured(69), false),
            terminal("tty_help_BaselineFilter")
        );
        assert_eq!(
            usage::<BaselineFilter>(coloured(UNSHAPED), false),
            terminal("tty_help_nosize")
        );
        assert_eq!(
            usage::<BaselineFilter>(coloured(49), false),
            terminal("tty_help_columns")
        );
        // Standard error redirected, standard input on the terminal: shaped
        // to the terminal, not coloured, and nothing reaches the terminal.
        assert!(oracle("tty_help_errfile", "tty.bin").is_empty());
        assert_eq!(
            usage::<BaselineFilter>(
                Console {
                    width: 69,
                    colour: false
                },
                false
            ),
            oracle("tty_help_errfile", "stderr.txt")
        );
    }

    #[cfg(feature = "mzml")]
    #[test]
    fn subsections_on_a_terminal_match_the_release_build() {
        use crate::cli::tools::PeakPickerHiRes;
        assert_eq!(
            usage::<PeakPickerHiRes>(coloured(69), false),
            terminal("tty_help_PeakPickerHiRes")
        );
        assert_eq!(
            usage::<PeakPickerHiRes>(coloured(69), true),
            terminal("tty_helphelp_PeakPickerHiRes")
        );
    }

    #[cfg(all(feature = "mzml", feature = "featurexml"))]
    #[test]
    fn file_info_usage_on_a_terminal_matches_the_release_build() {
        use crate::cli::tools::FileInfo;
        assert_eq!(
            usage::<FileInfo>(coloured(69), false),
            terminal("tty_help_FileInfo")
        );
    }
}

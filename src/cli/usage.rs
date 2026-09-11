// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Usage text, the native form of the source `printUsage_`.

use super::parameter::ParameterType;
use super::spec::ToolSpec;
use crate::param::ParamValue;
use std::io::{Result as IoResult, Write};

/// Column at which descriptions start, as in the source layout.
const DESCRIPTION_COLUMN: usize = 34;

/// Source rendering of a default value, using the parameter text conversion
/// rather than Rust's own float formatting.
fn default_text(value: &ParamValue) -> String {
    value.to_text(false).unwrap_or_default()
}

/// Wrap `text` to `width`, indenting every line after the first by `indent`.
fn wrap(text: &str, width: usize, indent: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();
    for word in text.split_whitespace() {
        if !current.is_empty() && current.len() + 1 + word.len() > width {
            lines.push(std::mem::take(&mut current));
        }
        if !current.is_empty() {
            current.push(' ');
        }
        current.push_str(word);
    }
    if !current.is_empty() {
        lines.push(current);
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    for line in lines.iter_mut().skip(1) {
        *line = format!("{}{line}", " ".repeat(indent));
    }
    lines
}

/// Print the usage block. `advanced` includes parameters registered as advanced,
/// which is what the source `--helphelp` flag selects.
pub fn print(
    out: &mut dyn Write,
    name: &str,
    description: &str,
    spec: &ToolSpec,
    advanced: bool,
) -> IoResult<()> {
    writeln!(out)?;
    writeln!(out, "{name} -- {description}")?;
    writeln!(out, "Version: {}", crate::CORE_SDK_VERSION)?;
    writeln!(out)?;
    writeln!(out, "Usage:")?;
    writeln!(out, "  {name} <options>")?;
    writeln!(out)?;
    writeln!(out, "Options (mandatory options marked with '*'):")?;

    let width = 96usize.saturating_sub(DESCRIPTION_COLUMN);
    let mut hidden = 0usize;
    for entry in spec.parameters() {
        match entry.kind {
            ParameterType::Newline => {
                writeln!(out)?;
                continue;
            }
            ParameterType::Text => {
                writeln!(out, "{}", entry.description)?;
                continue;
            }
            _ => {}
        }
        if entry.advanced && !advanced {
            hidden += 1;
            continue;
        }
        let mut left = format!(
            "  {}{}",
            entry.token(),
            if entry.argument.is_empty() {
                String::new()
            } else {
                format!(" {}", entry.argument)
            }
        );
        if entry.required {
            left.push('*');
        }
        let mut right = entry.description.clone();
        let default = default_text(&entry.default_value);
        if !default.is_empty() && !matches!(entry.kind, ParameterType::Flag) {
            right.push_str(&format!(" (default: '{default}')"));
        }
        if !entry.valid_strings.is_empty() {
            right.push_str(&format!(" (valid: {})", entry.valid_strings.join(", ")));
        }
        if !entry.valid_formats.is_empty() {
            right.push_str(&format!(" (formats: {})", entry.valid_formats.join(", ")));
        }
        let lines = wrap(&right, width, DESCRIPTION_COLUMN);
        if left.len() >= DESCRIPTION_COLUMN {
            writeln!(out, "{left}")?;
            for line in &lines {
                writeln!(
                    out,
                    "{}{}",
                    " ".repeat(DESCRIPTION_COLUMN),
                    line.trim_start()
                )?;
            }
        } else {
            writeln!(
                out,
                "{left}{}{}",
                " ".repeat(DESCRIPTION_COLUMN - left.len()),
                lines[0]
            )?;
            for line in &lines[1..] {
                writeln!(out, "{line}")?;
            }
        }
    }
    if hidden > 0 {
        writeln!(out)?;
        writeln!(
            out,
            "You can list {hidden} more advanced parameter(s) with '--helphelp'."
        )?;
    }
    writeln!(out)?;
    Ok(())
}

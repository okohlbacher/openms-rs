// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Dedicated bounded readers; source row failures are reported and skipped.

use super::{
    MAX_BYTES, MAX_RECORDS, MAX_WORK, RibonucleotideDiagnostic, RibonucleotideEntry,
    RibonucleotideLoadReport, add, invalid,
};
use crate::chemistry::ribonucleotide::{
    MAX_RIBONUCLEOTIDE_CODE_BYTES, MAX_RIBONUCLEOTIDE_TEXT_BYTES, Ribonucleotide,
    RibonucleotideRecord, RibonucleotideTermSpecificity,
};
use crate::{Error, Result};
use std::sync::Arc;

const MAX_INPUT_BYTES: usize = 16 * 1024 * 1024;
const MAX_ROW_BYTES: usize = 1024 * 1024;
const HEADER: &str = "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass";

#[derive(Default)]
struct Limits {
    bytes: usize,
    work: usize,
    rows: usize,
}
impl Limits {
    fn new(text: &str) -> Result<Self> {
        if text.len() > MAX_INPUT_BYTES {
            return Err(invalid("RNA provider input exceeds 16 MiB"));
        }
        Ok(Self {
            work: text.len() * 2,
            ..Self::default()
        })
    }
    fn row(&mut self) -> Result<()> {
        add(&mut self.rows, 1, MAX_RECORDS, "RNA provider row count")
    }
    fn fields(&mut self, text: &[&str], code: &str, formulas: &[&str]) -> Result<()> {
        if code.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES
            || text.iter().any(|s| s.len() > MAX_RIBONUCLEOTIDE_TEXT_BYTES)
        {
            return Err(invalid("RNA provider text/code field limit exceeded"));
        }
        for formula in formulas {
            if formula.len() > MAX_RIBONUCLEOTIDE_TEXT_BYTES {
                return Err(invalid("RNA formula text limit exceeded"));
            }
            add(
                &mut self.work,
                formula.len().saturating_add(1).saturating_mul(256),
                MAX_WORK,
                "RNA formula work",
            )?;
        }
        Ok(())
    }
    fn output(
        &mut self,
        report: &mut RibonucleotideLoadReport,
        index: usize,
        entry: Result<Option<RibonucleotideEntry>>,
    ) -> Result<()> {
        match entry {
            Ok(Some(entry)) => {
                let extra = entry
                    .alternatives
                    .as_ref()
                    .map_or(0, |a| a[0].len() + a[1].len());
                add(
                    &mut self.bytes,
                    entry
                        .ribonucleotide
                        .payload_bytes()?
                        .saturating_add(extra)
                        .saturating_add(96),
                    MAX_BYTES,
                    "RNA provider output",
                )?;
                report.entries.push(entry);
            }
            other => {
                let message = match other {
                    Err(error) => error.to_string(),
                    _ => "empty ribonucleotide code omitted".into(),
                };
                add(
                    &mut self.bytes,
                    message.len() + 64,
                    MAX_BYTES,
                    "RNA provider diagnostics",
                )?;
                report.skipped += 1;
                report
                    .diagnostics
                    .push(RibonucleotideDiagnostic { index, message });
            }
        }
        Ok(())
    }
}

/// Source TSV schema, including zero-versus-missing mass behavior and original
/// branch order. Input is bounded to 16 MiB, 100k lines and 1 MiB per data row;
/// record/output/formula-work limits are fatal. Ordinary malformed rows are skipped.
pub fn read_tsv(text: &str) -> Result<RibonucleotideLoadReport> {
    let mut limits = Limits::new(text)?;
    let mut rows = text.lines().enumerate();
    let header = loop {
        let Some((index, row)) = rows.next() else {
            return Err(invalid("RNA TSV header is missing"));
        };
        limits.row()?;
        if !row.starts_with('#') {
            break (index, row);
        }
    };
    if !header.1.starts_with(HEADER) {
        return Err(Error::Parse {
            line: header.0 + 1,
            message: "unexpected RNA TSV header".into(),
        });
    }
    let mut report = RibonucleotideLoadReport::default();
    for (index, row) in rows {
        limits.row()?;
        if row.len() > MAX_ROW_BYTES {
            return Err(invalid("RNA TSV row byte limit exceeded"));
        }
        // Replace PRIME throughout the row before interpreting any field.
        let normalized = row.replace('\u{2032}', "'");
        let fields: Vec<&str> = normalized.split('\t').take(10).collect();
        if fields.len() >= 9 {
            let code = fields[1]
                .strip_suffix("tRNA")
                .filter(|_| fields[1].ends_with("QtRNA"))
                .unwrap_or(fields[1]);
            limits.fields(&[fields[0], fields[2], fields[5]], code, &[fields[6]])?;
            if let Some(alternatives) = fields.get(9).filter(|_| {
                !fields[2].ends_with('N') && !fields[1].starts_with('d') && fields[1].ends_with('?')
            }) {
                // Only first/last space-delimited pieces are used by source.
                if let Some((first, _)) = alternatives.split_once(' ') {
                    let last = alternatives.rsplit_once(' ').unwrap().1;
                    if first.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES
                        || last.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES
                    {
                        return Err(invalid("RNA alternative code byte limit exceeded"));
                    }
                }
            }
        }
        limits.output(&mut report, index + 1, parse_tsv_row(&fields))?;
    }
    Ok(report)
}

fn parse_tsv_row(fields: &[&str]) -> Result<Option<RibonucleotideEntry>> {
    if fields.len() < 9 {
        return Err(invalid("RNA TSV requires at least nine fields"));
    }
    if fields[1].is_empty() {
        return Err(invalid("RNA TSV code is empty"));
    }
    let mut record = RibonucleotideRecord {
        name: fields[0].into(),
        code: if fields[1].ends_with("QtRNA") {
            fields[1][..fields[1].len() - 4].into()
        } else {
            fields[1].into()
        },
        new_code: fields[2].into(),
        html_code: fields[5].into(),
        ..Default::default()
    };
    if fields[3] == "preQ0base" {
        record.origin = 'G';
    } else if fields[3].len() == 1 {
        record.origin = fields[3].as_bytes()[0] as char;
    }
    if !fields[6].is_empty() && fields[6] != "-" {
        record.formula = fields[6].parse()?;
    }
    if !fields[7].is_empty() && fields[7] != "None" {
        record.mono_mass = number(fields[7])?;
        if record.mono_mass == 0.0 && !record.formula.is_empty() {
            record.mono_mass = record.formula.mono_mass();
        }
    }
    if !fields[8].is_empty() && fields[8] != "None" {
        record.average_mass = number(fields[8])?;
        if record.average_mass == 0.0 && !record.formula.is_empty() {
            record.average_mass = record.formula.average_mass();
        }
    }
    let mut alternatives = None;
    if fields[2].ends_with('N') {
        if fields[2].contains("55") || fields[2] == "N" {
            record.term_specificity = RibonucleotideTermSpecificity::FivePrime;
        } else if fields[2].contains("33") {
            record.term_specificity = RibonucleotideTermSpecificity::ThreePrime;
        }
    } else if fields[1].starts_with('d') {
        record.baseloss_formula = "C5H10O4".parse()?;
    } else if fields[1].ends_with('m') || fields[1].ends_with("m*") {
        record.baseloss_formula = "C6H12O5".parse()?;
    } else if fields[1].ends_with('?') {
        let value = fields
            .get(9)
            .ok_or_else(|| invalid("RNA ambiguity requires a tenth field"))?;
        let (first, _) = value
            .split_once(' ')
            .ok_or_else(|| invalid("RNA ambiguity requires space-separated alternatives"))?;
        let last = value.rsplit_once(' ').unwrap().1;
        alternatives = Some([first.into(), last.into()]);
    } else if matches!(fields[1], "Ar(p)" | "Gr(p)") {
        record.baseloss_formula = "C10H19O21P".parse()?;
    }
    Ok(Some(RibonucleotideEntry {
        ribonucleotide: Arc::new(Ribonucleotide::from_record(record)?),
        alternatives,
    }))
}
fn number(text: &str) -> Result<f64> {
    let value = text
        .trim_ascii()
        .parse::<f64>()
        .map_err(|_| invalid("invalid RNA mass value"))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("RNA mass must be finite"))
    }
}

/// MODOMICS object keys iterate lexically; JSON arrays retain input order.
/// Bounds before JSON allocation: 16 MiB input, depth 64, 250k structural tokens
/// and 393,216 raw bytes per string. No renormalization of declared masses occurs.
#[cfg(feature = "rna-json")]
pub fn read_modomics_json(text: &str) -> Result<RibonucleotideLoadReport> {
    let mut limits = Limits::new(text)?;
    json_preflight(text)?;
    let document: serde_json::Value = serde_json::from_str(text).map_err(|error| Error::Parse {
        line: error.line(),
        message: format!("RNA JSON: {error}"),
    })?;
    let elements: Box<dyn Iterator<Item = &serde_json::Value> + '_> = match &document {
        serde_json::Value::Object(object) => Box::new(object.values()),
        serde_json::Value::Array(array) => Box::new(array.iter()),
        serde_json::Value::Null => Box::new(std::iter::empty()),
        other => Box::new(std::iter::once(other)),
    };
    let mut report = RibonucleotideLoadReport::default();
    for (index, element) in elements.enumerate() {
        limits.row()?;
        // Bound consumed string fields before converting into owned records.
        if let Some(object) = element.as_object() {
            let text = ["name", "short_name", "abbrev"].map(|key| {
                object
                    .get(key)
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("")
            });
            let formulas = ["formula", "baseloss_formula"].map(|key| {
                object
                    .get(key)
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("")
            });
            limits.fields(&text, text[1], &formulas)?;
            if let Some(alternatives) = object
                .get("alternatives")
                .and_then(serde_json::Value::as_array)
                .filter(|_| text[1].ends_with('?') || text[1].ends_with("?*"))
            {
                if alternatives
                    .iter()
                    .take(2)
                    .filter_map(serde_json::Value::as_str)
                    .any(|code| code.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES)
                {
                    return Err(invalid("RNA alternative code byte limit exceeded"));
                }
            }
        }
        limits.output(&mut report, index + 1, parse_json_row(element))?;
    }
    Ok(report)
}

#[cfg(feature = "rna-json")]
fn parse_json_row(value: &serde_json::Value) -> Result<Option<RibonucleotideEntry>> {
    use serde_json::Value;
    let object = value
        .as_object()
        .ok_or_else(|| invalid("RNA JSON entry must be an object"))?;
    let string = |key| {
        object
            .get(key)
            .and_then(Value::as_str)
            .ok_or_else(|| invalid("RNA JSON required string is missing or has wrong type"))
    };
    let code = string("short_name")?;
    let mut record = RibonucleotideRecord {
        name: string("name")?.into(),
        code: code.into(),
        new_code: code.into(),
        formula: string("formula")?.parse()?,
        ..Default::default()
    };
    let moieties = object
        .get("reference_moiety")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid("RNA reference_moiety must be an array"))?;
    if moieties.len() == 1 && moieties[0].as_str().is_some_and(|m| m.len() == 1) {
        record.origin = moieties[0].as_str().unwrap().as_bytes()[0] as char;
    } else if moieties.len() == 4 {
        record.origin = 'X';
        if code.ends_with("pN") {
            record.term_specificity = RibonucleotideTermSpecificity::FivePrime;
        } else if code.starts_with('N') && code.ends_with('p') {
            record.term_specificity = RibonucleotideTermSpecificity::ThreePrime;
        }
    } else {
        return Err(invalid("unsupported RNA reference moieties"));
    }
    if object.contains_key("abbrev") {
        record.html_code = string("abbrev")?.into();
    }
    if let Some(mass) = object.get("mass_avg").filter(|m| !m.is_null()) {
        record.average_mass = mass
            .as_f64()
            .ok_or_else(|| invalid("RNA average mass must be numeric"))?;
    }
    record.mono_mass = if let Some(mass) = object.get("mass_monoiso").filter(|m| !m.is_null()) {
        mass.as_f64()
            .ok_or_else(|| invalid("RNA mono mass must be numeric"))?
    } else {
        record.formula.mono_mass()
    };
    record.baseloss_formula =
        if let Some(value) = object.get("baseloss_formula").filter(|m| !m.is_null()) {
            value
                .as_str()
                .ok_or_else(|| invalid("RNA base-loss formula must be a string"))?
                .parse()?
        } else if code.starts_with('d') {
            "C5H10O4".parse()?
        } else if code.ends_with('m') || code.ends_with("m*") {
            "C6H12O5".parse()?
        } else if code.ends_with("Ar(p)") || code.ends_with("Gr(p)") {
            "C10H19O21P".parse()?
        } else {
            "C5H10O5".parse()?
        };
    let alternatives = if code.ends_with('?') || code.ends_with("?*") {
        let values = object
            .get("alternatives")
            .and_then(Value::as_array)
            .filter(|a| a.len() >= 2)
            .ok_or_else(|| invalid("RNA ambiguity requires two alternative strings"))?;
        let first = values[0]
            .as_str()
            .ok_or_else(|| invalid("RNA alternative must be a string"))?;
        let second = values[1]
            .as_str()
            .ok_or_else(|| invalid("RNA alternative must be a string"))?;
        // The JSON provider copies alternatives only when its first string is
        // nonempty; the TSV provider writes its two fields directly instead.
        (!first.is_empty()).then(|| [first.into(), second.into()])
    } else {
        None
    };
    if code.is_empty() {
        return Ok(None);
    }
    Ok(Some(RibonucleotideEntry {
        ribonucleotide: Arc::new(Ribonucleotide::from_record(record)?),
        alternatives,
    }))
}

#[cfg(feature = "rna-json")]
fn json_preflight(text: &str) -> Result<()> {
    let (mut quoted, mut escaped, mut depth, mut tokens, mut string_bytes) =
        (false, false, 0usize, 0usize, 0usize);
    for byte in text.bytes() {
        if quoted {
            string_bytes += 1;
            if string_bytes > 6 * MAX_RIBONUCLEOTIDE_TEXT_BYTES {
                return Err(invalid("RNA JSON raw string byte limit exceeded"));
            }
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                quoted = false;
            }
        } else {
            match byte {
                b'"' => {
                    quoted = true;
                    string_bytes = 0;
                }
                b'[' | b'{' => {
                    depth += 1;
                    tokens += 1;
                }
                b']' | b'}' => {
                    depth = depth.saturating_sub(1);
                }
                b',' | b':' => {
                    tokens += 1;
                }
                _ => {}
            }
            if depth > 64 || tokens > 250_000 {
                return Err(invalid("RNA JSON depth/token limit exceeded"));
            }
        }
    }
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenMS version-one escaped definition records.
use super::{ModificationRecord, NeutralLoss, ResidueModification, TermSpecificity};
use crate::chemistry::EmpiricalFormula;
use crate::data_structures::list::ListParse;
use crate::{Error, Result};

pub(crate) const MAX_DEFINITION_BYTES: usize = 16 * 1024 * 1024;
pub(crate) const MAX_DEFINITION_RECORDS: usize = 100_000;
fn invalid(text: &str) -> Error {
    Error::InvalidValue(text.into())
}
fn bounded(text: &str) -> Result<()> {
    if text.len() > MAX_DEFINITION_BYTES {
        Err(invalid("modification definition byte limit exceeded"))
    } else {
        Ok(())
    }
}

impl ResidueModification {
    /// Source version-one field projection. Fields absent from that format (for
    /// example absolute composition or aliases) are not included. Interchange
    /// adapters use the stricter lossless `modification_definitions::encode`.
    pub fn to_definition_string(&self) -> Result<String> {
        if self.name().is_empty() {
            return Err(invalid("definition requires a nonempty name"));
        }
        if self.payload_bytes()? > MAX_DEFINITION_BYTES / 2 {
            return Err(invalid("definition payload limit exceeded"));
        }
        let losses = self
            .neutral_losses()
            .iter()
            .map(|loss| loss.formula().to_string())
            .collect::<Vec<_>>()
            .join(",");
        let fields = [
            "1".to_owned(),
            self.name().into(),
            self.full_id().into(),
            self.full_name().into(),
            self.origin().unwrap_or('X').to_string(),
            if self.term_specificity() == TermSpecificity::Anywhere {
                "none".into()
            } else {
                self.term_specificity().name().into()
            },
            if self.diff_formula().atoms.is_empty() {
                String::new()
            } else {
                self.diff_formula().to_string()
            },
            mass_text(self.diff_mono_mass()),
            mass_text(self.diff_average_mass()),
            losses,
        ];
        let mut output = String::new();
        for (index, field) in fields.iter().enumerate() {
            if index > 0 {
                output.push('|');
            }
            for c in field.chars() {
                if matches!(c, '\\' | '|' | ';') {
                    output.push('\\');
                }
                output.push(c);
            }
        }
        bounded(&output)?;
        Ok(output)
    }
    /// Decode an owned definition. Explicit delta masses are retained even when
    /// a formula exists. Empty numeric fields mean zero; loss masses derive from
    /// their formulas. Unknown extra fields are rejected rather than discarded.
    pub fn from_definition_string(text: &str) -> Result<Self> {
        bounded(text)?;
        let mut fields = vec![String::new()];
        let mut chars = text.chars();
        while let Some(c) = chars.next() {
            match c {
                '\\' => fields
                    .last_mut()
                    .unwrap()
                    .push(chars.next().unwrap_or('\\')),
                '|' => {
                    if fields.len() == 10 {
                        return Err(invalid("too many modification definition fields"));
                    }
                    fields.push(String::new());
                }
                _ => fields.last_mut().unwrap().push(c),
            }
        }
        if fields.len() < 9 || fields[0] != "1" || fields[1].is_empty() {
            return Err(invalid("invalid modification definition version or fields"));
        }
        let mut origin = fields[4].chars();
        let origin = origin
            .next()
            .filter(|c| c.is_ascii_uppercase())
            .filter(|_| origin.next().is_none())
            .ok_or_else(|| invalid("invalid definition origin"))?;
        let term = match fields[5].as_str() {
            "none" => TermSpecificity::Anywhere,
            "N-term" => TermSpecificity::NTerm,
            "C-term" => TermSpecificity::CTerm,
            "Protein N-term" => TermSpecificity::ProteinNTerm,
            "Protein C-term" => TermSpecificity::ProteinCTerm,
            _ => return Err(invalid("invalid definition specificity")),
        };
        let mass = |s: &str| -> Result<f64> {
            let value = if s.is_empty() {
                0.0
            } else {
                f64::from_list_item(s)?
            };
            if value.is_finite() {
                Ok(value)
            } else {
                Err(invalid("definition masses must be finite"))
            }
        };
        let mut losses = Vec::new();
        if let Some(text) = fields.get(9) {
            for value in text.split(',').filter(|s| !s.is_empty()) {
                if losses.len() == MAX_DEFINITION_RECORDS {
                    return Err(invalid("definition neutral loss limit exceeded"));
                }
                let formula = EmpiricalFormula::parse(value)?;
                losses.push(NeutralLoss::new(
                    formula.clone(),
                    formula.mono_mass(),
                    formula.average_mass(),
                )?);
            }
        }
        Self::from_record(ModificationRecord {
            name: fields[1].clone(),
            full_id: fields[2].clone(),
            full_name: fields[3].clone(),
            origin: Some(origin),
            term_specificity: term,
            diff_formula: EmpiricalFormula::parse(&fields[6])?,
            diff_mono_mass: mass(&fields[7])?,
            diff_average_mass: mass(&fields[8])?,
            neutral_losses: losses,
            ..Default::default()
        })
    }
    /// Split semicolon-separated records without stripping field escapes.
    pub fn split_definition_records(text: &str) -> Result<Vec<&str>> {
        bounded(text)?;
        let mut records = Vec::new();
        let mut start = 0;
        let mut escaped = false;
        for (index, c) in text.char_indices() {
            if escaped {
                escaped = false;
            } else if c == '\\' {
                escaped = true;
            } else if c == ';' {
                if index > start {
                    if records.len() == MAX_DEFINITION_RECORDS {
                        return Err(invalid("definition record limit exceeded"));
                    }
                    records.push(&text[start..index]);
                }
                start = index + 1;
            }
        }
        if start < text.len() {
            if records.len() == MAX_DEFINITION_RECORDS {
                return Err(invalid("definition record limit exceeded"));
            }
            records.push(&text[start..]);
        }
        Ok(records)
    }
}

fn mass_text(value: f64) -> String {
    let fixed = value.to_string();
    let raw = format!("{value:e}");
    let (mantissa, exponent) = raw.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    let scientific = format!(
        "{mantissa}e{}{:02}",
        if exponent < 0 { '-' } else { '+' },
        exponent.unsigned_abs()
    );
    if scientific.len() < fixed.len() {
        scientific
    } else {
        fixed
    }
}

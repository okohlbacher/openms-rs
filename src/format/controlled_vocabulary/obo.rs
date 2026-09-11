// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::*;
use std::ops::Bound::{Excluded, Unbounded};

pub(super) const EMBEDDED: &[(&str, &[u8], OboEncoding)] = &[
    (
        "MS",
        include_bytes!("../../../resources/cv/psi-ms.obo"),
        OboEncoding::Utf8,
    ),
    (
        "PATO",
        include_bytes!("../../../resources/cv/quality.obo"),
        OboEncoding::Utf8,
    ),
    (
        "UO",
        include_bytes!("../../../resources/cv/unit.obo"),
        OboEncoding::Utf8,
    ),
    (
        "BTO",
        include_bytes!("../../../resources/cv/brenda.obo"),
        OboEncoding::Windows1252,
    ),
    (
        "GO",
        include_bytes!("../../../resources/cv/goslim_goa.obo"),
        OboEncoding::Utf8,
    ),
];

fn copy(value: &str, m: &mut Meter<'_>) -> Result<String> {
    m.text(value.len())?;
    let mut result = String::new();
    result.try_reserve_exact(value.len()).map_err(|_| limit())?;
    result.push_str(value);
    Ok(result)
}
fn trim(value: &str) -> &str {
    value.trim_matches([' ', '\t', '\n', '\r'])
}
fn after(value: &str, character: char) -> &str {
    value.split_once(character).map_or(value, |(_, tail)| tail)
}
fn suffix(value: &str, character: char) -> &str {
    value.rsplit_once(character).map_or(value, |(_, tail)| tail)
}
fn prefix(value: &str, character: char) -> &str {
    value.split_once(character).map_or(value, |(head, _)| head)
}
fn take(entries: &mut usize) -> Result<()> {
    *entries = entries.checked_sub(1).ok_or_else(limit)?;
    Ok(())
}

pub(super) fn load(
    cv: &mut ControlledVocabulary,
    name: &str,
    mut reader: impl BufRead,
    encoding: OboEncoding,
    m: &mut Meter<'_>,
) -> Result<OboLoadReport> {
    cv.name = copy(name, m)?;
    let limits = cv.limits;
    let mut input_left = limits.max_input_bytes;
    let mut entries = limits.max_entries;
    let mut raw = Vec::new();
    let mut current = CVTermDefinition::default();
    let mut report = OboLoadReport::default();
    let mut in_term = false;
    let mut line_number = 0usize;
    while read_line(
        &mut reader,
        &mut raw,
        &mut input_left,
        limits.max_line_bytes,
        m,
    )? {
        line_number = add(line_number, 1)?;
        // Repeated lexical scans have a fixed bounded count in the projection.
        m.spend(mul(raw.len(), 24)?)?;
        let decoded = decode(&raw, encoding, line_number, m)?;
        let line = trim(&decoded);
        if line.is_empty() {
            continue;
        }
        m.text(line.len())?;
        let mut compact = String::new();
        compact.try_reserve_exact(line.len()).map_err(|_| limit())?;
        compact.extend(
            line.chars()
                .filter(|c| !matches!(c, ' ' | '\t' | '\n' | '\r')),
        );
        if compact.starts_with("data-version:") {
            cv.version = copy(trim(after(line, ':')), m)?;
        }
        if compact.starts_with("default-namespace:") {
            cv.label = copy(trim(after(line, ':')), m)?;
        }
        if compact.starts_with("remark:URL:") {
            if let Some(position) = line.find("http://").or_else(|| line.find("https://")) {
                cv.url = copy(trim(&line[position..]), m)?;
            } else {
                diagnostic(
                    &mut report,
                    line_number,
                    "No URL found in the line.",
                    m,
                    &mut entries,
                )?;
            }
        }
        if compact.starts_with('[') {
            if compact.eq_ignore_ascii_case("[term]") {
                in_term = true;
                commit(cv, &mut current, &mut report, m, &mut entries)?;
                current = CVTermDefinition::default();
            } else {
                in_term = false;
            }
            continue;
        }
        if !in_term {
            continue;
        }
        take(&mut entries)?;
        if compact.starts_with("id:") {
            current.id = copy(trim(after(line, ':')), m)?;
        } else if compact.starts_with("name:") {
            current.name = copy(trim(after(line, ':')), m)?;
        } else if compact.starts_with("is_a:") {
            let id = trim(prefix(after(line, ':'), '!'));
            check_name(
                cv,
                id,
                line,
                &current.id,
                "parent",
                line_number,
                &mut report,
                m,
                &mut entries,
            )?;
            insert_set(&mut current.parents, id, m)?;
        } else if (compact.starts_with("relationship:DRV")
            || compact.starts_with("relationship:part_of"))
            && name == "brenda"
        {
            let relation = if compact.starts_with("relationship:DRV") {
                "DRV"
            } else {
                "part_of"
            };
            let id = relationship_id(line, relation, m)?;
            check_name(
                cv,
                &id,
                line,
                &current.id,
                relation,
                line_number,
                &mut report,
                m,
                &mut entries,
            )?;
            insert_set(&mut current.parents, &id, m)?;
        } else if compact.starts_with("relationship:has_units") {
            let id = relationship_id(line, "has_units", m)?;
            check_name(
                cv,
                &id,
                line,
                &current.id,
                "has_units",
                line_number,
                &mut report,
                m,
                &mut entries,
            )?;
            insert_set(&mut current.units, &id, m)?;
        } else if compact.starts_with("def:") {
            current.description = copy(quoted(line), m)?;
        } else if compact.starts_with("synonym:") {
            push(&mut current.synonyms, quoted(line), m)?;
        } else if compact == "is_obsolete:true" {
            current.obsolete = true;
        } else if compact.starts_with("xref:value-type")
            || compact.starts_with("xref_analog:value-type")
        {
            compact.retain(|c| c != '\\');
            if let Some(kind) = value_type(&compact, false) {
                current.xref_type = kind;
            } else {
                diagnostic(
                    &mut report,
                    line_number,
                    "unknown xsd value type, ignoring",
                    m,
                    &mut entries,
                )?;
            }
        } else if compact.starts_with("relationship:has_value_type") {
            if let Some(kind) = value_type(&compact, true) {
                current.xref_type = kind;
            } else {
                diagnostic(
                    &mut report,
                    line_number,
                    "unknown xsd value type, ignoring",
                    m,
                    &mut entries,
                )?;
            }
        } else if compact.starts_with("xref:binary-data-type")
            || compact.starts_with("xref_analog:binary-data-type")
        {
            compact.retain(|c| c != '\\');
            // Correct source's fixed22 for the29-byte xref_analog prefix.
            let length = if compact.starts_with("xref_analog:") {
                29
            } else {
                22
            };
            let value = prefix(&compact, '"');
            let value = value
                .get(length.min(value.len())..)
                .ok_or_else(|| Error::Parse {
                    line: line_number,
                    message: "binary-data-type byte offset splits UTF-8".into(),
                })?;
            push(&mut current.xref_binary, trim(value), m)?;
        } else {
            push(&mut current.unparsed, line, m)?;
        }
    }
    commit(cv, &mut current, &mut report, m, &mut entries)?;
    indexes(cv, m, &mut entries)?;
    Ok(report)
}

fn read_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    input_left: &mut usize,
    max_line: usize,
    m: &mut Meter<'_>,
) -> Result<bool> {
    line.clear();
    loop {
        let data = reader.fill_buf()?;
        if data.is_empty() {
            return Ok(!line.is_empty());
        }
        let allowed = (*input_left).min(max_line.saturating_sub(line.len()));
        // Inspect at most the remaining bound plus a one-byte excess sentinel.
        let inspected = data.len().min(allowed.saturating_add(1));
        // Charge each examined byte before touching it. Charging the entire
        // fill_buf slice per short line would penalize a borrowed whole file
        // quadratically compared with an ordinary buffered file reader.
        let mut end = None;
        for (i, byte) in data[..inspected].iter().enumerate() {
            m.spend(1)?;
            if *byte == b'\n' {
                end = Some(i);
                break;
            }
        }
        let count = end.map_or(inspected, |i| i + 1);
        if count > allowed {
            return Err(limit());
        }
        *input_left -= count;
        if line.capacity() - line.len() < count {
            let capacity = add(line.len(), count)?
                .max(line.capacity().saturating_mul(2))
                .max(64)
                .min(max_line);
            m.text(capacity)?;
            line.try_reserve_exact(capacity - line.len())
                .map_err(|_| limit())?;
        }
        m.spend(count)?;
        line.extend_from_slice(&data[..count]);
        reader.consume(count);
        if end.is_some() {
            return Ok(true);
        }
    }
}

fn decode(bytes: &[u8], encoding: OboEncoding, line: usize, m: &mut Meter<'_>) -> Result<String> {
    if encoding == OboEncoding::Utf8 {
        return copy(
            std::str::from_utf8(bytes).map_err(|_| Error::Parse {
                line,
                message: "invalid UTF-8 OBO text".into(),
            })?,
            m,
        );
    }
    // Undefined Windows-1252 bytes are explicit errors; never replacement text.
    const SPECIAL: [u32; 32] = [
        0x20ac, 0, 0x201a, 0x0192, 0x201e, 0x2026, 0x2020, 0x2021, 0x02c6, 0x2030, 0x0160, 0x2039,
        0x0152, 0, 0x017d, 0, 0, 0x2018, 0x2019, 0x201c, 0x201d, 0x2022, 0x2013, 0x2014, 0x02dc,
        0x2122, 0x0161, 0x203a, 0x0153, 0, 0x017e, 0x0178,
    ];
    let capacity = mul(bytes.len(), 3)?;
    m.text(capacity)?;
    let mut result = String::new();
    result.try_reserve_exact(capacity).map_err(|_| limit())?;
    for byte in bytes {
        let code = if (0x80..0xa0).contains(byte) {
            let code = SPECIAL[usize::from(*byte - 0x80)];
            if code == 0 {
                return Err(Error::Parse {
                    line,
                    message: "undefined Windows-1252 OBO byte".into(),
                });
            }
            code
        } else {
            u32::from(*byte)
        };
        result.push(char::from_u32(code).unwrap());
    }
    Ok(result)
}
fn quoted(line: &str) -> &str {
    trim(prefix(trim(after(line, '"')), '"'))
}
fn relationship_id(line: &str, relation: &str, m: &mut Meter<'_>) -> Result<String> {
    // Source size_t npos arithmetic and StringUtils::substr clamping are
    // defined even for oddly spaced or truncated relationship keywords.
    let start = line
        .find(relation)
        .unwrap_or(usize::MAX)
        .wrapping_add(relation.len() + 1)
        .min(line.len());
    let head = prefix(
        line.get(start..)
            .ok_or_else(|| invalid("OBO relationship byte offset splits UTF-8"))?,
        ':',
    );
    let tail = trim(prefix(suffix(line, ':'), '!'));
    let capacity = add(add(head.len(), 1)?, tail.len())?;
    m.text(capacity)?;
    let mut result = String::new();
    result.try_reserve_exact(capacity).map_err(|_| limit())?;
    result.push_str(head);
    result.push(':');
    result.push_str(tail);
    Ok(result)
}
fn value_type(text: &str, relationship: bool) -> Option<XRefType> {
    let needle = |kind: &str| {
        if relationship {
            text.contains(kind)
        } else {
            text.match_indices("value-type:")
                .any(|(i, _)| text[i + 11..].starts_with(kind))
        }
    };
    for (names, kind) in [
        (&["xsd:string"][..], XRefType::String),
        (&["xsd:integer", "xsd:int"][..], XRefType::Integer),
        (
            &["xsd:decimal", "xsd:float", "xsd:double"][..],
            XRefType::Decimal,
        ),
        (&["xsd:negativeInteger"][..], XRefType::NegativeInteger),
        (&["xsd:positiveInteger"][..], XRefType::PositiveInteger),
        (
            &["xsd:nonNegativeInteger"][..],
            XRefType::NonNegativeInteger,
        ),
        (
            &["xsd:nonPositiveInteger"][..],
            XRefType::NonPositiveInteger,
        ),
        (&["xsd:boolean", "xsd:bool"][..], XRefType::Boolean),
        (&["xsd:date"][..], XRefType::Date),
        (&["xsd:anyURI"][..], XRefType::AnyUri),
    ] {
        if names.iter().any(|n| needle(n)) {
            return Some(kind);
        }
    }
    if relationship
        && ["MS:1002711", "MS:1002712", "MS:1002713"]
            .iter()
            .any(|n| text.contains(n))
    {
        Some(XRefType::String)
    } else {
        None
    }
}
fn diagnostic(
    report: &mut OboLoadReport,
    line: usize,
    message: &str,
    m: &mut Meter<'_>,
    entries: &mut usize,
) -> Result<()> {
    take(entries)?;
    reserve(&mut report.diagnostics, m)?;
    let message = copy(message, m)?;
    report.diagnostics.push(OboDiagnostic { line, message });
    Ok(())
}
#[allow(clippy::too_many_arguments)]
fn check_name(
    cv: &ControlledVocabulary,
    id: &str,
    line: &str,
    term: &str,
    relation: &str,
    number: usize,
    report: &mut OboLoadReport,
    m: &mut Meter<'_>,
    entries: &mut usize,
) -> Result<()> {
    if !line.contains('!') {
        return Ok(());
    }
    m.lookup(cv.terms.len(), cv.max_id_bytes, id.len())?;
    if let Some(known) = cv.terms.get(id) {
        let supplied = trim(suffix(line, '!'));
        m.spend(add(known.name.len(), supplied.len())?)?;
        if !known.name.eq_ignore_ascii_case(supplied) {
            let capacity = add(add(term.len(), relation.len())?, 64)?;
            m.text(capacity)?;
            let message =
                format!("term '{term}': {relation} relationship name and identifier differ");
            diagnostic(report, number, &message, m, entries)?;
        }
    }
    Ok(())
}
fn reserve<T>(values: &mut Vec<T>, m: &mut Meter<'_>) -> Result<()> {
    m.spend(1)?;
    if values.len() == values.capacity() {
        let capacity = values.capacity().saturating_mul(2).max(4);
        m.slots::<T>(capacity)?;
        m.spend(values.len())?;
        values
            .try_reserve_exact(capacity - values.len())
            .map_err(|_| limit())?;
    }
    Ok(())
}
fn push(values: &mut Vec<String>, value: &str, m: &mut Meter<'_>) -> Result<()> {
    reserve(values, m)?;
    let value = copy(value, m)?;
    values.push(value);
    Ok(())
}
fn insert_set(set: &mut BTreeSet<String>, value: &str, m: &mut Meter<'_>) -> Result<()> {
    // Every key comparison is bounded by the query's own length.
    m.lookup(set.len(), usize::MAX, value.len())?;
    m.tree::<String>(1)?;
    set.insert(copy(value, m)?);
    Ok(())
}
fn commit(
    cv: &mut ControlledVocabulary,
    term: &mut CVTermDefinition,
    report: &mut OboLoadReport,
    m: &mut Meter<'_>,
    entries: &mut usize,
) -> Result<()> {
    if term.id.is_empty() {
        return Ok(());
    }
    take(entries)?;
    m.lookup(cv.terms.len(), cv.max_id_bytes, term.id.len())?;
    if cv.terms.len() >= cv.limits.max_terms && !cv.terms.contains_key(&term.id) {
        return Err(limit());
    }
    m.tree::<(String, CVTermDefinition)>(1)?;
    let key = copy(&term.id, m)?;
    cv.max_id_bytes = cv.max_id_bytes.max(key.len());
    cv.terms.insert(key, std::mem::take(term));
    report.definitions = add(report.definitions, 1)?;
    Ok(())
}
fn indexes(cv: &mut ControlledVocabulary, m: &mut Meter<'_>, entries: &mut usize) -> Result<()> {
    let mut current = match cv.terms.first_key_value() {
        None => return Ok(()),
        Some((id, _)) => copy(id, m)?,
    };
    loop {
        m.lookup(cv.terms.len(), cv.max_id_bytes, current.len())?;
        let term = cv.get_term(&current)?;
        // Mutating the map can insert source placeholders; own only this small
        // descriptor snapshot, not a full cloned term or registry.
        m.slots::<String>(term.parents.len())?;
        let mut parents = Vec::new();
        parents
            .try_reserve_exact(term.parents.len())
            .map_err(|_| limit())?;
        for parent in &term.parents {
            parents.push(copy(parent, m)?);
        }
        let name = copy(&term.name, m)?;
        let description = copy(&term.description, m)?;
        for parent in parents {
            take(entries)?;
            m.lookup(cv.terms.len(), cv.max_id_bytes, parent.len())?;
            if !cv.terms.contains_key(&parent) {
                if cv.terms.len() >= cv.limits.max_terms {
                    return Err(limit());
                }
                m.tree::<(String, CVTermDefinition)>(1)?;
                cv.max_id_bytes = cv.max_id_bytes.max(parent.len());
                m.lookup(cv.terms.len(), cv.max_id_bytes, parent.len())?;
                cv.terms
                    .insert(copy(&parent, m)?, CVTermDefinition::default());
            }
            m.lookup(cv.terms.len(), cv.max_id_bytes, parent.len())?;
            insert_set(
                &mut cv.terms.get_mut(&parent).unwrap().children,
                &current,
                m,
            )?;
        }
        take(entries)?;
        m.lookup(cv.names.len(), cv.max_name_bytes, name.len())?;
        let key = if cv.names.contains_key(&name) {
            let n = add(name.len(), description.len())?;
            m.text(n)?;
            let mut joined = String::new();
            joined.try_reserve_exact(n).map_err(|_| limit())?;
            joined.push_str(&name);
            joined.push_str(&description);
            joined
        } else {
            name
        };
        m.lookup(cv.names.len(), cv.max_name_bytes, key.len())?;
        m.tree::<(String, String)>(1)?;
        cv.max_name_bytes = cv.max_name_bytes.max(key.len());
        cv.names.entry(key).or_insert(copy(&current, m)?);
        // Source ++it occurs after insertions, so newly inserted later keys
        // participate in this pass, while earlier placeholders do not.
        m.lookup(cv.terms.len(), cv.max_id_bytes, current.len())?;
        let next = cv
            .terms
            .range::<str, _>((Excluded(current.as_str()), Unbounded))
            .next();
        match next {
            None => break,
            Some((id, _)) => current = copy(id, m)?,
        }
    }
    Ok(())
}

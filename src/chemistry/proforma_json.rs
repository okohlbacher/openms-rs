// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::*;
use serde_json::{Map, Number, Value};
use std::{io, mem::size_of};

/// Maximum input or output JSON text bytes per operation.
pub const MAX_PROFORMA_JSON_TEXT_BYTES: usize = 4 * 1024 * 1024;
/// Maximum cumulative input structural tokens and consumed collection items.
pub const MAX_PROFORMA_JSON_ITEMS: usize = 1_000_000;
/// Shared scanning, comparison, copying and conversion work per operation.
pub const MAX_PROFORMA_JSON_WORK: usize = 50_000_000;
/// Conservative cumulative requested-capacity allowance, not physical memory.
pub const MAX_PROFORMA_JSON_BYTES: usize = 256 * 1024 * 1024;
/// Maximum JSON container depth, including ignored fields.
pub const MAX_PROFORMA_JSON_DEPTH: usize = 64;

impl Peptidoform {
    /// Write the source tagged JSON schema; resolved chemistry is omitted.
    pub fn to_json(&self) -> Result<String> {
        let mut work = Work::new();
        let value = write_chain(self, &mut work)?;
        encode(&value, &mut work)
    }
    /// Read the source JSON schema without resolving chemistry or text grammar.
    pub fn from_json(input: &str) -> Result<Self> {
        let mut work = Work::new();
        let value = decode(input, &mut work)?;
        read_chain(&value, &mut work)
    }
}
impl PeptidoformIon {
    /// Write all source fields, including ion name and per-chain charges.
    pub fn to_json(&self) -> Result<String> {
        let mut work = Work::new();
        let value = write_ion(self, &mut work)?;
        encode(&value, &mut work)
    }
    /// Read an ion atomically; missing is_chimeric defaults to false.
    pub fn from_json(input: &str) -> Result<Self> {
        let mut work = Work::new();
        let value = decode(input, &mut work)?;
        read_ion(&value, &mut work)
    }
}
fn bad(message: &'static str) -> Error {
    invalid(message)
}
struct Work {
    remaining: usize,
    bytes: usize,
    items: usize,
}
impl Work {
    fn new() -> Self {
        Self {
            remaining: MAX_PROFORMA_JSON_WORK,
            bytes: MAX_PROFORMA_JSON_BYTES,
            items: MAX_PROFORMA_JSON_ITEMS,
        }
    }
    fn consume(&mut self, n: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or_else(|| bad("ProForma JSON work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, n: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(n)
            .ok_or_else(|| bad("ProForma JSON capacity limit exceeded"))?;
        Ok(())
    }
    fn items(&mut self, n: usize) -> Result<()> {
        self.items = self
            .items
            .checked_sub(n)
            .ok_or_else(|| bad("ProForma JSON item limit exceeded"))?;
        self.consume(
            n.checked_mul(64)
                .ok_or_else(|| bad("ProForma JSON size overflow"))?,
        )
    }
    fn vector<T>(&mut self, n: usize) -> Result<Vec<T>> {
        self.items(n)?;
        self.allocate(
            n.checked_mul(size_of::<T>())
                .ok_or_else(|| bad("ProForma JSON size overflow"))?,
        )?;
        let mut v = Vec::new();
        v.try_reserve_exact(n)
            .map_err(|_| bad("Cannot allocate ProForma JSON collection"))?;
        Ok(v)
    }
    fn string(&mut self, s: &str) -> Result<String> {
        self.consume(s.len())?;
        self.allocate(s.len())?;
        let mut out = String::new();
        out.try_reserve_exact(s.len())
            .map_err(|_| bad("Cannot allocate ProForma JSON string"))?;
        out.push_str(s);
        Ok(out)
    }
}
fn decode(input: &str, w: &mut Work) -> Result<Value> {
    if input.len() > MAX_PROFORMA_JSON_TEXT_BYTES {
        return Err(bad("ProForma JSON input limit exceeded"));
    }
    w.consume(input.len().saturating_mul(4))?;
    let (mut quoted, mut escaped, mut depth, mut items) = (false, false, 0usize, 1usize);
    for byte in input.bytes() {
        if quoted {
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
                    items += 1;
                }
                b'[' | b'{' => {
                    depth += 1;
                    items += 1;
                }
                b']' | b'}' => depth = depth.saturating_sub(1),
                b',' | b':' => items += 1,
                _ => (),
            }
            if depth > MAX_PROFORMA_JSON_DEPTH || items > MAX_PROFORMA_JSON_ITEMS {
                return Err(bad("ProForma JSON depth/item limit exceeded"));
            }
        }
    }
    w.items(items)?;
    // Covers sorted object-key comparisons and serde's DOM/string scratch before
    // allocation, including duplicate keys and ignored fields. This is a
    // conservative bound over raw text, not a second JSON syntax parser.
    let levels = (usize::BITS - items.leading_zeros()) as usize;
    w.consume(input.len().saturating_mul(levels + 2))?;
    w.allocate(
        items
            .saturating_mul(512)
            .saturating_add(input.len().saturating_mul(4)),
    )?;
    let normalized = normalize_integer_zero(input, w)?;
    let input = normalized.as_deref().unwrap_or(input);
    // nlohmann accepts a leading UTF-8 BOM; serde_json does not.
    let input = input.strip_prefix('\u{feff}').unwrap_or(input);
    serde_json::from_str(input).map_err(|e| Error::Parse {
        line: e.line(),
        message: format!("ProForma JSON: {e}"),
    })
}
// Preserve the source integer/float distinction: nlohmann parses bare -0 as
// integer zero while serde_json promotes it to float -0.0. This small lexical
// normalization touches no strings or decimal/exponent tokens and keeps offsets.
fn normalize_integer_zero(input: &str, w: &mut Work) -> Result<Option<String>> {
    let mut output = None;
    let (mut quoted, mut escaped) = (false, false);
    w.consume(input.len())?;
    for (index, byte) in input.bytes().enumerate() {
        if quoted {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                quoted = false;
            }
        } else if byte == b'"' {
            quoted = true;
        } else if byte == b'-'
            && (index == 0
                || (index == 3 && input.starts_with('\u{feff}'))
                || matches!(
                    input.as_bytes()[index - 1],
                    b' ' | b'\n' | b'\r' | b'\t' | b':' | b',' | b'['
                ))
            && input.as_bytes().get(index + 1) == Some(&b'0')
            && input
                .as_bytes()
                .get(index + 2)
                .is_none_or(|b| matches!(b, b' ' | b'\n' | b'\r' | b'\t' | b',' | b']' | b'}'))
        {
            if output.is_none() {
                output = Some(w.string(input)?.into_bytes());
            }
            output.as_mut().unwrap()[index] = b' ';
        }
    }
    output
        .map(|v| String::from_utf8(v).map_err(|_| bad("Invalid normalized ProForma JSON")))
        .transpose()
}
struct Output<'a> {
    bytes: Vec<u8>,
    work: &'a mut Work,
    failure: Option<Error>,
}
impl io::Write for Output<'_> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let result = (|| {
            self.work.consume(data.len().saturating_add(1))?;
            let wanted = self
                .bytes
                .len()
                .checked_add(data.len())
                .ok_or_else(|| bad("ProForma JSON output size overflow"))?;
            if wanted > MAX_PROFORMA_JSON_TEXT_BYTES {
                return Err(bad("ProForma JSON output limit exceeded"));
            }
            if wanted > self.bytes.capacity() {
                let cap = wanted
                    .max(self.bytes.capacity().saturating_mul(2))
                    .clamp(128, MAX_PROFORMA_JSON_TEXT_BYTES);
                self.work.consume(self.bytes.len())?;
                self.work.allocate(cap)?;
                self.bytes
                    .try_reserve_exact(cap - self.bytes.len())
                    .map_err(|_| bad("Cannot allocate ProForma JSON output"))?;
            }
            self.bytes.extend_from_slice(data);
            Ok(data.len())
        })();
        result.map_err(|error| {
            self.failure = Some(error);
            io::Error::other("ProForma JSON output limit")
        })
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn encode(value: &Value, w: &mut Work) -> Result<String> {
    let mut output = Output {
        bytes: Vec::new(),
        work: w,
        failure: None,
    };
    if serde_json::to_writer(&mut output, value).is_err() {
        return Err(output
            .failure
            .unwrap_or_else(|| bad("Cannot encode ProForma JSON")));
    }
    String::from_utf8(output.bytes).map_err(|_| bad("Invalid ProForma JSON output UTF-8"))
}
fn required<'a>(v: &'a Value, key: &str) -> Result<&'a Value> {
    v.get(key)
        .ok_or_else(|| bad("Missing required ProForma JSON field or invalid object"))
}
fn optional<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    v.get(key).filter(|v| !v.is_null())
}
fn text(v: &Value) -> Result<&str> {
    v.as_str()
        .ok_or_else(|| bad("ProForma JSON field must be a string"))
}
fn read_text(v: &Value, w: &mut Work) -> Result<String> {
    w.string(text(v)?)
}
fn bool_value(v: &Value) -> Result<bool> {
    v.as_bool()
        .ok_or_else(|| bad("ProForma JSON field must be boolean"))
}
fn number(v: &Value) -> Result<f64> {
    v.as_f64()
        .filter(|x| x.is_finite())
        .ok_or_else(|| bad("ProForma JSON field must be a finite number"))
}
fn integer(v: &Value) -> Result<i32> {
    // Source generic int overload accepts bool; its double overload does not.
    if let Some(x) = v.as_bool() {
        return Ok(i32::from(x));
    }
    if let Some(x) = v.as_i64() {
        return i32::try_from(x).map_err(|_| bad("ProForma JSON integer exceeds i32"));
    }
    if let Some(x) = v.as_u64() {
        return i32::try_from(x).map_err(|_| bad("ProForma JSON integer exceeds i32"));
    }
    let x = number(v)?.trunc();
    if x < f64::from(i32::MIN) || x > f64::from(i32::MAX) {
        return Err(bad("ProForma JSON integer exceeds i32"));
    }
    Ok(x as i32)
}
fn variant<'a>(v: &'a Value, w: &mut Work) -> Result<(&'a str, &'a Value)> {
    let kind = text(required(v, "type")?)?;
    w.consume(kind.len().saturating_mul(8))?;
    Ok((kind, required(v, "value")?))
}
fn read_values<T>(
    v: &Value,
    strict: bool,
    mut f: impl FnMut(&Value, &mut Work) -> Result<T>,
    w: &mut Work,
) -> Result<Vec<T>> {
    if strict && !v.is_array() {
        return Err(bad("ProForma JSON field must be an array"));
    }
    let n = match v {
        Value::Array(a) => a.len(),
        Value::Object(a) => a.len(),
        Value::Null => 0,
        _ => 1,
    };
    let mut out = w.vector(n)?;
    match v {
        Value::Array(a) => {
            for x in a {
                out.push(f(x, w)?);
            }
        }
        Value::Object(a) => {
            for x in a.values() {
                out.push(f(x, w)?);
            }
        }
        Value::Null => (),
        x => out.push(f(x, w)?),
    }
    Ok(out)
}
fn read_cv(v: &Value, w: &mut Work) -> Result<CvDatabase> {
    let s = text(v)?;
    w.consume(s.len().saturating_mul(5))?;
    match s {
        "UNIMOD" => Ok(CvDatabase::Unimod),
        "MOD" => Ok(CvDatabase::Mod),
        "RESID" => Ok(CvDatabase::Resid),
        "XLMOD" => Ok(CvDatabase::Xlmod),
        "GNO" => Ok(CvDatabase::Gno),
        _ => Err(bad("Unknown ProForma JSON CV database")),
    }
}
fn read_formula(v: &Value, w: &mut Work) -> Result<FormulaTag> {
    Ok(FormulaTag {
        formula_string: read_text(required(v, "formula_string")?, w)?,
        charge: optional(v, "charge").map(integer).transpose()?,
    })
}
fn read_label(v: &Value, w: &mut Work) -> Result<Label> {
    let kind = text(required(v, "type")?)?;
    w.consume(kind.len().saturating_mul(3))?;
    let label_type = match kind {
        "CROSSLINK" => LabelType::Crosslink,
        "BRANCH" => LabelType::Branch,
        "AMBIGUOUS" => LabelType::Ambiguous,
        _ => return Err(bad("Unknown ProForma JSON label type")),
    };
    Ok(Label {
        label_type,
        identifier: read_text(required(v, "identifier")?, w)?,
        score: optional(v, "score").map(number).transpose()?,
    })
}
fn read_tag(v: &Value, w: &mut Work) -> Result<ModificationTag> {
    let (kind, v) = variant(v, w)?;
    Ok(match kind {
        "cv_accession" => ModificationTag::CvAccession(CvAccession {
            database: read_cv(required(v, "database")?, w)?,
            accession: read_text(required(v, "accession")?, w)?,
        }),
        "named_mod" => ModificationTag::NamedMod(NamedMod {
            name: read_text(required(v, "name")?, w)?,
            cv_hint: optional(v, "cv_hint").map(|x| read_cv(x, w)).transpose()?,
        }),
        "mass_delta" => {
            let s = text(required(v, "source")?)?;
            w.consume(s.len().saturating_mul(7))?;
            let source = match s {
                "NONE" => MassDeltaSource::None,
                "OBS" => MassDeltaSource::Obs,
                "U" => MassDeltaSource::U,
                "M" => MassDeltaSource::M,
                "R" => MassDeltaSource::R,
                "X" => MassDeltaSource::X,
                "G" => MassDeltaSource::G,
                _ => return Err(bad("Unknown ProForma JSON mass source")),
            };
            ModificationTag::MassDelta(MassDelta {
                source,
                mass: number(required(v, "mass")?)?,
                original_text: read_text(required(v, "original_text")?, w)?,
            })
        }
        "formula" => ModificationTag::FormulaTag(read_formula(v, w)?),
        "glycan" => ModificationTag::GlycanComposition(GlycanComposition {
            components: read_values(
                v,
                false,
                |item, w| {
                    let (kind, value) = variant(required(item, "monosaccharide")?, w)?;
                    let mono = match kind {
                        "name" => GlycanComponent::Name(read_text(value, w)?),
                        "formula" => GlycanComponent::Formula(read_formula(value, w)?),
                        _ => return Err(bad("Unknown ProForma JSON glycan type")),
                    };
                    Ok((mono, integer(required(item, "count")?)?))
                },
                w,
            )?,
        }),
        "info" => ModificationTag::InfoTag(InfoTag {
            text: read_text(required(v, "text")?, w)?,
        }),
        "position" => {
            let s = text(required(v, "residues")?)?;
            w.consume(s.len())?;
            if !s.is_ascii() {
                return Err(bad("ProForma JSON residue characters must be ASCII bytes"));
            }
            let mut residues = w.vector(s.len())?;
            residues.extend(s.chars());
            ModificationTag::PositionConstraint(PositionConstraint {
                residues,
                n_term: v
                    .get("n_term")
                    .map(bool_value)
                    .transpose()?
                    .unwrap_or(false),
                c_term: v
                    .get("c_term")
                    .map(bool_value)
                    .transpose()?
                    .unwrap_or(false),
            })
        }
        _ => return Err(bad("Unknown ProForma JSON modification tag")),
    })
}
fn read_modification(v: &Value, w: &mut Work) -> Result<Modification> {
    Ok(Modification {
        alternatives: read_values(
            v,
            false,
            |alt, w| {
                Ok((
                    read_tag(required(alt, "tag")?, w)?,
                    optional(alt, "label")
                        .map(|v| read_label(v, w))
                        .transpose()?,
                ))
            },
            w,
        )?,
        resolved_mod: None,
    })
}
fn read_element(v: &Value, w: &mut Work) -> Result<SequenceElement> {
    let aa = text(required(v, "amino_acid")?)?;
    if aa.len() != 1 {
        return Err(bad("ProForma JSON amino_acid must be exactly one byte"));
    }
    Ok(SequenceElement {
        amino_acid: char::from(aa.as_bytes()[0]),
        modifications: read_values(required(v, "modifications")?, true, read_modification, w)?,
    })
}
fn read_section(v: &Value, w: &mut Work) -> Result<SequenceSection> {
    let (kind, v) = variant(v, w)?;
    Ok(match kind {
        "element" => SequenceSection::Element(read_element(v, w)?),
        "ambiguous_region" => SequenceSection::AmbiguousRegion(AmbiguousRegion {
            elements: read_values(required(v, "elements")?, true, read_element, w)?,
        }),
        "modified_range" => SequenceSection::ModifiedRange(ModifiedRange {
            elements: read_values(required(v, "elements")?, true, read_element, w)?,
            modifications: read_values(required(v, "modifications")?, true, read_modification, w)?,
        }),
        _ => return Err(bad("Unknown ProForma JSON sequence section")),
    })
}
fn read_global(v: &Value, w: &mut Work) -> Result<GlobalModEntry> {
    let (kind, v) = variant(v, w)?;
    Ok(match kind {
        "isotope_replacement" => GlobalModEntry::IsotopeReplacement(IsotopeReplacement {
            isotope: read_text(required(v, "isotope")?, w)?,
        }),
        "global_modification" => GlobalModEntry::GlobalModification(GlobalModification {
            modification: read_modification(required(v, "modification")?, w)?,
            locations: read_values(required(v, "locations")?, false, read_text, w)?,
        }),
        _ => return Err(bad("Unknown ProForma JSON global entry")),
    })
}
fn read_charge(v: &Value, w: &mut Work) -> Result<ChargeState> {
    let (kind, v) = variant(v, w)?;
    Ok(match kind {
        "simple" => ChargeState::Simple(integer(v)?),
        "adducts" => ChargeState::Adducts(read_values(
            v,
            true,
            |v, w| {
                Ok(AdductIon {
                    formula: read_text(required(v, "formula")?, w)?,
                    charge: integer(required(v, "charge")?)?,
                    occurrence: optional(v, "occurrence").map(integer).transpose()?,
                })
            },
            w,
        )?),
        _ => return Err(bad("Unknown ProForma JSON charge type")),
    })
}
fn read_chain(v: &Value, w: &mut Work) -> Result<Peptidoform> {
    Ok(Peptidoform {
        global_mods: v
            .get("global_mods")
            .map(|v| read_values(v, true, read_global, w))
            .transpose()?
            .unwrap_or_default(),
        unlocalised_mods: read_values(
            required(v, "unlocalised_mods")?,
            true,
            |v, w| {
                Ok(UnlocalisedMod {
                    modifications: read_values(
                        required(v, "modifications")?,
                        true,
                        read_modification,
                        w,
                    )?,
                    occurrence: optional(v, "occurrence").map(integer).transpose()?,
                })
            },
            w,
        )?,
        labile_mods: read_values(
            required(v, "labile_mods")?,
            true,
            |v, w| {
                Ok(LabileModification {
                    modification: read_modification(required(v, "modification")?, w)?,
                })
            },
            w,
        )?,
        n_term_mods: read_values(required(v, "n_term_mods")?, true, read_modification, w)?,
        sequence: read_values(required(v, "sequence")?, true, read_section, w)?,
        c_term_mods: read_values(required(v, "c_term_mods")?, true, read_modification, w)?,
        name: optional(v, "name").map(|v| read_text(v, w)).transpose()?,
        charge: optional(v, "charge")
            .map(|v| read_charge(v, w))
            .transpose()?,
    })
}
fn read_ion(v: &Value, w: &mut Work) -> Result<PeptidoformIon> {
    Ok(PeptidoformIon {
        chains: read_values(required(v, "chains")?, true, read_chain, w)?,
        name: optional(v, "name").map(|v| read_text(v, w)).transpose()?,
        charge: optional(v, "charge")
            .map(|v| read_charge(v, w))
            .transpose()?,
        is_chimeric: v
            .get("is_chimeric")
            .map(bool_value)
            .transpose()?
            .unwrap_or(false),
    })
}

fn string(v: &str, w: &mut Work) -> Result<Value> {
    Ok(Value::String(w.string(v)?))
}
fn real(v: f64) -> Result<Value> {
    Number::from_f64(v)
        .map(Value::Number)
        .ok_or_else(|| bad("ProForma JSON stored floating value must be finite"))
}
// Charge a whole sparse BTree node per inserted field: 32 key/value slots
// plus 1KiB for node headers/edges and allocator headroom. This deliberately
// exceeds a typical leaf even when the object contains only one field.
const MAP_ENTRY_BYTES: usize = 1024 + 32 * size_of::<(String, Value)>();
fn object<const N: usize>(fields: [(&str, Value); N], w: &mut Work) -> Result<Value> {
    w.items(N + 1)?;
    w.allocate(N.saturating_mul(MAP_ENTRY_BYTES))?;
    let mut map = Map::new();
    for (key, value) in fields {
        map.insert(w.string(key)?, value);
    }
    Ok(Value::Object(map))
}
fn field(object: &mut Value, key: &str, value: Value, w: &mut Work) -> Result<()> {
    w.items(1)?;
    w.allocate(MAP_ENTRY_BYTES)?;
    object
        .as_object_mut()
        .unwrap()
        .insert(w.string(key)?, value);
    Ok(())
}
fn array<T>(
    values: &[T],
    mut f: impl FnMut(&T, &mut Work) -> Result<Value>,
    w: &mut Work,
) -> Result<Value> {
    let mut array = w.vector(values.len())?;
    for value in values {
        array.push(f(value, w)?);
    }
    Ok(Value::Array(array))
}
fn tagged(kind: &str, value: Value, w: &mut Work) -> Result<Value> {
    object([("type", string(kind, w)?), ("value", value)], w)
}
fn cv_name(v: CvDatabase) -> &'static str {
    match v {
        CvDatabase::Unimod => "UNIMOD",
        CvDatabase::Mod => "MOD",
        CvDatabase::Resid => "RESID",
        CvDatabase::Xlmod => "XLMOD",
        CvDatabase::Gno => "GNO",
    }
}
fn write_formula(v: &FormulaTag, w: &mut Work) -> Result<Value> {
    let mut out = object([("formula_string", string(&v.formula_string, w)?)], w)?;
    if let Some(charge) = v.charge {
        field(&mut out, "charge", charge.into(), w)?;
    }
    Ok(out)
}
fn write_label(v: &Label, w: &mut Work) -> Result<Value> {
    let kind = match v.label_type {
        LabelType::Crosslink => "CROSSLINK",
        LabelType::Branch => "BRANCH",
        LabelType::Ambiguous => "AMBIGUOUS",
    };
    let mut out = object(
        [
            ("type", string(kind, w)?),
            ("identifier", string(&v.identifier, w)?),
        ],
        w,
    )?;
    if let Some(score) = v.score {
        field(&mut out, "score", real(score)?, w)?;
    }
    Ok(out)
}
fn write_tag(v: &ModificationTag, w: &mut Work) -> Result<Value> {
    let (kind, value) = match v {
        ModificationTag::CvAccession(v) => (
            "cv_accession",
            object(
                [
                    ("database", string(cv_name(v.database), w)?),
                    ("accession", string(&v.accession, w)?),
                ],
                w,
            )?,
        ),
        ModificationTag::NamedMod(v) => {
            let mut out = object([("name", string(&v.name, w)?)], w)?;
            if let Some(cv) = v.cv_hint {
                field(&mut out, "cv_hint", string(cv_name(cv), w)?, w)?;
            }
            ("named_mod", out)
        }
        ModificationTag::MassDelta(v) => {
            let source = match v.source {
                MassDeltaSource::None => "NONE",
                MassDeltaSource::Obs => "OBS",
                MassDeltaSource::U => "U",
                MassDeltaSource::M => "M",
                MassDeltaSource::R => "R",
                MassDeltaSource::X => "X",
                MassDeltaSource::G => "G",
            };
            (
                "mass_delta",
                object(
                    [
                        ("source", string(source, w)?),
                        ("mass", real(v.mass)?),
                        ("original_text", string(&v.original_text, w)?),
                    ],
                    w,
                )?,
            )
        }
        ModificationTag::FormulaTag(v) => ("formula", write_formula(v, w)?),
        ModificationTag::GlycanComposition(v) => (
            "glycan",
            array(
                &v.components,
                |(mono, count), w| {
                    let mono = match mono {
                        GlycanComponent::Name(name) => tagged("name", string(name, w)?, w)?,
                        GlycanComponent::Formula(formula) => {
                            tagged("formula", write_formula(formula, w)?, w)?
                        }
                    };
                    object([("monosaccharide", mono), ("count", (*count).into())], w)
                },
                w,
            )?,
        ),
        ModificationTag::InfoTag(v) => ("info", object([("text", string(&v.text, w)?)], w)?),
        ModificationTag::PositionConstraint(v) => {
            w.items(v.residues.len())?;
            w.allocate(v.residues.len())?;
            let mut residues = String::new();
            residues
                .try_reserve_exact(v.residues.len())
                .map_err(|_| bad("Cannot allocate ProForma JSON residues"))?;
            for &ch in &v.residues {
                if !ch.is_ascii() {
                    return Err(bad("ProForma JSON residue characters must be ASCII bytes"));
                }
                residues.push(ch);
            }
            (
                "position",
                object(
                    [
                        ("residues", Value::String(residues)),
                        ("n_term", v.n_term.into()),
                        ("c_term", v.c_term.into()),
                    ],
                    w,
                )?,
            )
        }
    };
    tagged(kind, value, w)
}
fn write_modification(v: &Modification, w: &mut Work) -> Result<Value> {
    // Source JSON deliberately omits the external resolved pointer. Do not
    // traverse or clone its potentially large owned native Arc payload.
    array(
        &v.alternatives,
        |(tag, label), w| {
            let mut out = object([("tag", write_tag(tag, w)?)], w)?;
            if let Some(label) = label {
                field(&mut out, "label", write_label(label, w)?, w)?;
            }
            Ok(out)
        },
        w,
    )
}
fn write_element(v: &SequenceElement, w: &mut Work) -> Result<Value> {
    if !v.amino_acid.is_ascii() {
        return Err(bad("ProForma JSON amino_acid must be one ASCII byte"));
    }
    let mut byte = [0];
    let aa = v.amino_acid.encode_utf8(&mut byte);
    object(
        [
            ("amino_acid", string(aa, w)?),
            (
                "modifications",
                array(&v.modifications, write_modification, w)?,
            ),
        ],
        w,
    )
}
fn write_section(v: &SequenceSection, w: &mut Work) -> Result<Value> {
    let (kind, value) = match v {
        SequenceSection::Element(v) => ("element", write_element(v, w)?),
        SequenceSection::AmbiguousRegion(v) => (
            "ambiguous_region",
            object([("elements", array(&v.elements, write_element, w)?)], w)?,
        ),
        SequenceSection::ModifiedRange(v) => (
            "modified_range",
            object(
                [
                    ("elements", array(&v.elements, write_element, w)?),
                    (
                        "modifications",
                        array(&v.modifications, write_modification, w)?,
                    ),
                ],
                w,
            )?,
        ),
    };
    tagged(kind, value, w)
}
fn write_global(v: &GlobalModEntry, w: &mut Work) -> Result<Value> {
    let (kind, value) = match v {
        GlobalModEntry::IsotopeReplacement(v) => (
            "isotope_replacement",
            object([("isotope", string(&v.isotope, w)?)], w)?,
        ),
        GlobalModEntry::GlobalModification(v) => (
            "global_modification",
            object(
                [
                    ("modification", write_modification(&v.modification, w)?),
                    ("locations", array(&v.locations, |v, w| string(v, w), w)?),
                ],
                w,
            )?,
        ),
    };
    tagged(kind, value, w)
}
fn write_charge(v: &ChargeState, w: &mut Work) -> Result<Value> {
    let (kind, value) = match v {
        ChargeState::Simple(v) => ("simple", (*v).into()),
        ChargeState::Adducts(v) => (
            "adducts",
            array(
                v,
                |v, w| {
                    let mut out = object(
                        [
                            ("formula", string(&v.formula, w)?),
                            ("charge", v.charge.into()),
                        ],
                        w,
                    )?;
                    if let Some(count) = v.occurrence {
                        field(&mut out, "occurrence", count.into(), w)?;
                    }
                    Ok(out)
                },
                w,
            )?,
        ),
    };
    tagged(kind, value, w)
}
fn write_chain(v: &Peptidoform, w: &mut Work) -> Result<Value> {
    let mut out = object(
        [
            ("global_mods", array(&v.global_mods, write_global, w)?),
            (
                "unlocalised_mods",
                array(
                    &v.unlocalised_mods,
                    |v, w| {
                        let mut out = object(
                            [(
                                "modifications",
                                array(&v.modifications, write_modification, w)?,
                            )],
                            w,
                        )?;
                        if let Some(count) = v.occurrence {
                            field(&mut out, "occurrence", count.into(), w)?;
                        }
                        Ok(out)
                    },
                    w,
                )?,
            ),
            (
                "labile_mods",
                array(
                    &v.labile_mods,
                    |v, w| {
                        object(
                            [("modification", write_modification(&v.modification, w)?)],
                            w,
                        )
                    },
                    w,
                )?,
            ),
            ("n_term_mods", array(&v.n_term_mods, write_modification, w)?),
            ("sequence", array(&v.sequence, write_section, w)?),
            ("c_term_mods", array(&v.c_term_mods, write_modification, w)?),
        ],
        w,
    )?;
    if let Some(name) = &v.name {
        field(&mut out, "name", string(name, w)?, w)?;
    }
    if let Some(charge) = &v.charge {
        field(&mut out, "charge", write_charge(charge, w)?, w)?;
    }
    Ok(out)
}
fn write_ion(v: &PeptidoformIon, w: &mut Work) -> Result<Value> {
    let mut out = object(
        [
            ("chains", array(&v.chains, write_chain, w)?),
            ("is_chimeric", v.is_chimeric.into()),
        ],
        w,
    )?;
    if let Some(name) = &v.name {
        field(&mut out, "name", string(name, w)?, w)?;
    }
    if let Some(charge) = &v.charge {
        field(&mut out, "charge", write_charge(charge, w)?, w)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_one_field_objects_charge_a_complete_map_node() {
        let mut work = Work::new();
        work.bytes = 512;
        assert!(object([("name", Value::Null)], &mut work).is_err());
        let mut work = Work::new();
        work.bytes = MAP_ENTRY_BYTES + 4;
        assert!(object([("name", Value::Null)], &mut work).is_ok());
        assert_eq!(work.bytes, 0);
    }

    #[test]
    fn raw_duplicates_and_ignored_fields_are_precharged_before_dom_creation() {
        let mut work = Work::new();
        work.items = 3;
        assert!(decode(r#"{"ignored":1,"ignored":2}"#, &mut work).is_err());
        let mut work = Work::new();
        work.bytes = 1;
        assert!(decode(r#"{"ignored":[]}"#, &mut work).is_err());
    }

    #[test]
    fn all_chains_and_output_share_remaining_work_and_capacity() {
        let ion = PeptidoformIon {
            chains: vec![Peptidoform::default(); 3],
            ..Default::default()
        };
        let mut complete = Work::new();
        let value = write_ion(&ion, &mut complete).unwrap();
        let before_output = complete.remaining;
        encode(&value, &mut complete).unwrap();
        let mut bounded = Work::new();
        bounded.remaining = MAX_PROFORMA_JSON_WORK - before_output - 1;
        assert!(write_ion(&ion, &mut bounded).is_err());
        let mut bounded = Work::new();
        bounded.bytes = 0;
        assert!(encode(&value, &mut bounded).is_err());
    }

    #[test]
    fn lexical_preflight_counts_containers_but_not_quoted_escapes() {
        let mut work = Work::new();
        let input = format!(
            "{}null{}",
            "[".repeat(MAX_PROFORMA_JSON_DEPTH),
            "]".repeat(MAX_PROFORMA_JSON_DEPTH)
        );
        assert!(decode(&input, &mut work).is_ok());
        let input = format!("[{input}]");
        assert!(decode(&input, &mut Work::new()).is_err());
        let quoted = serde_json::to_string(&"[\\\"-0]".repeat(100)).unwrap();
        assert_eq!(
            decode(&quoted, &mut Work::new()).unwrap().as_str(),
            Some("[\\\"-0]".repeat(100).as_str())
        );
    }
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Bounded header-only output buffer. External writers receive it only after
//! the complete registry/reference/representation preflight has succeeded.
use super::*;

pub(super) struct Xml<'a> {
    pub text: String,
    ids: BTreeSet<String>,
    pub work: &'a mut Work,
}
impl<'a> Xml<'a> {
    pub fn new(work: &'a mut Work) -> Self {
        Self {
            text: String::new(),
            ids: BTreeSet::new(),
            work,
        }
    }
    pub fn raw(&mut self, text: &str) -> Result<()> {
        self.work.charge(text.len(), 0)?;
        let required = self
            .text
            .len()
            .checked_add(text.len())
            .ok_or_else(resource)?;
        if required > self.text.capacity() {
            let capacity = self
                .text
                .capacity()
                .checked_mul(2)
                .ok_or_else(resource)?
                .max(required)
                .max(128);
            self.work.charge(0, capacity)?;
            self.text
                .try_reserve_exact(capacity - self.text.len())
                .map_err(|_| resource())?;
        }
        self.text.push_str(text);
        Ok(())
    }
    fn escaped(&mut self, text: &str) -> Result<()> {
        self.work.charge(text.len(), 0)?;
        xml_string(text)?;
        for c in text.chars() {
            match c {
                '&' => self.raw("&amp;")?,
                '<' => self.raw("&lt;")?,
                '>' => self.raw("&gt;")?,
                '"' => self.raw("&quot;")?,
                '\'' => self.raw("&apos;")?,
                '\n' => self.raw("&#10;")?,
                '\r' => self.raw("&#13;")?,
                '\t' => self.raw("&#9;")?,
                c => {
                    let mut bytes = [0; 4];
                    self.raw(c.encode_utf8(&mut bytes))?;
                }
            }
        }
        Ok(())
    }
    pub fn attribute(&mut self, name: &str, value: &str) -> Result<()> {
        self.raw(" ")?;
        self.raw(name)?;
        self.raw("=\"")?;
        self.escaped(value)?;
        self.raw("\"")
    }
    pub fn define_id(&mut self, id: &str) -> Result<()> {
        self.work.charge(id.len().saturating_mul(64), 0)?;
        if parameter_id(id)? != id {
            return Err(invalid(
                "native header ID must not contain surrounding whitespace",
            ));
        }
        self.work.meter().tree::<String>(1)?;
        let owned = self.work.copy(id)?;
        if !self.ids.insert(owned) {
            return Err(invalid("duplicate generated or native header ID"));
        }
        Ok(())
    }
    pub fn start(&mut self, tag: &str, attrs: &[(&str, &str)]) -> Result<()> {
        if let Some((_, id)) = attrs.iter().find(|(key, _)| *key == "id") {
            self.define_id(id)?;
        }
        self.raw("<")?;
        self.raw(tag)?;
        for &(name, value) in attrs {
            self.attribute(name, value)?;
        }
        self.raw(">\n")
    }
    pub fn end(&mut self, tag: &str) -> Result<()> {
        self.raw("</")?;
        self.raw(tag)?;
        self.raw(">\n")
    }
    pub fn integer(&mut self, value: usize) -> Result<String> {
        self.work.charge(64, 64)?;
        Ok(value.to_string())
    }
    pub fn signed(&mut self, value: i32) -> Result<String> {
        self.work.charge(64, 64)?;
        Ok(value.to_string())
    }
    pub fn scalar(&mut self, value: &MetaValue) -> Result<String> {
        match value.data() {
            MetaValueData::String(s) => self.work.copy(s),
            MetaValueData::Integer(n) => {
                self.work.charge(64, 64)?;
                Ok(n.to_string())
            }
            MetaValueData::Float(n) => {
                self.work.charge(1024, 1024)?;
                Ok(n.to_string())
            }
            _ => Err(Error::Unsupported(
                "Empty/list metadata has no lossless mzML scalar encoding".into(),
            )),
        }
    }
    pub fn cv(&mut self, id: &str, value: Option<&MetaValue>) -> Result<()> {
        let term = self.work.cv(id)?;
        let prefix = id.split_once(':').map_or("", |p| p.0);
        self.raw("<cvParam")?;
        self.attribute("cvRef", prefix)?;
        self.attribute("accession", id)?;
        self.attribute("name", &term.name)?;
        if let Some(value) = value {
            let scalar = self.scalar(value)?;
            self.attribute("value", &scalar)?;
            self.unit(value)?;
        }
        self.raw("/>\n")
    }
    fn unit(&mut self, value: &MetaValue) -> Result<()> {
        if let Some(unit) = value.unit() {
            product_unit(unit)?;
            // Retain the actual owned unit name. An omitted name remains
            // omitted instead of being silently replaced by dictionary text.
            self.work.cv(unit.accession())?;
            self.attribute("unitAccession", unit.accession())?;
            self.attribute("unitCvRef", unit.cv_ref())?;
            if !unit.name().is_empty() {
                self.attribute("unitName", unit.name())?;
            }
        }
        Ok(())
    }
    pub fn cv_text(&mut self, id: &str, value: &str) -> Result<()> {
        let owned = self.work.copy(value)?;
        self.cv(id, Some(&owned.into()))
    }
    pub fn cv_float(&mut self, id: &str, value: f64, unit: Option<&str>) -> Result<()> {
        let mut value = MetaValue::try_from(value)?;
        if let Some(id) = unit {
            let term = self.work.cv(id)?;
            self.work.charge(
                term.id.len().saturating_add(term.name.len()),
                term.id
                    .len()
                    .saturating_add(term.name.len())
                    .saturating_add(2),
            )?;
            let unit = Unit::new(&term.id, &term.name, id.split_once(':').unwrap().0)?;
            value = value.with_unit(unit)?;
        }
        self.cv(id, Some(&value))
    }
    pub fn user(&mut self, name: &str, value: &MetaValue) -> Result<()> {
        let scalar = self.scalar(value)?;
        let kind = match value.data() {
            MetaValueData::Integer(_) => "xsd:integer",
            MetaValueData::Float(_) => "xsd:double",
            _ => "xsd:string",
        };
        self.raw("<userParam")?;
        self.attribute("name", name)?;
        self.attribute("type", kind)?;
        self.attribute("value", &scalar)?;
        self.unit(value)?;
        self.raw("/>\n")
    }
    /// All CV entries precede all user entries. A name is promoted only if the
    /// exact read route returns that same metadata key/type; source promotion
    /// into an unrelated dedicated field would silently lose the metadata.
    pub fn metadata(&mut self, owner: &str, meta: &MetaInfo, exclude: &[&str]) -> Result<()> {
        let mut route = Vec::new();
        self.work.slots::<Option<&str>>(meta.len())?;
        route
            .try_reserve_exact(meta.len())
            .map_err(|_| resource())?;
        for (name, value) in meta {
            self.work.charge(
                name.len()
                    .saturating_mul(exclude.len().saturating_add(1))
                    .saturating_mul(3),
                0,
            )?;
            if exclude.contains(&name.as_str()) {
                route.push(None);
                continue;
            }
            if matches!(
                name.as_str(),
                "GO cellular component" | "brenda source tissue"
            ) && owner == "sample"
            {
                let text = value.as_str()?;
                let term = ControlledVocabulary::psi_ms()?
                    .find_term_by_name_with_budget(
                        text,
                        &mut self.work.remaining,
                        &mut self.work.bytes,
                    )?
                    .ok_or_else(|| invalid("sample ontology name is not in the pinned provider"))?;
                let expected = if name == "GO cellular component" {
                    "GO:"
                } else {
                    "BTO:"
                };
                if !term.id.starts_with(expected) || value.unit().is_some() {
                    return Err(Error::Unsupported(
                        "sample ontology key/value identity mismatch".into(),
                    ));
                }
                route.push(Some(term.id.as_str()));
                continue;
            }
            let term = ControlledVocabulary::psi_ms()?.find_term_by_name_with_budget(
                name,
                &mut self.work.remaining,
                &mut self.work.bytes,
            )?;
            let term = match term {
                Some(term)
                    if reversible_metadata(owner, name, &term.id)
                        && compatible_value(term.xref_type, name, value)
                        && mapping::permitted(owner, &term.id, self.work)? =>
                {
                    Some(term.id.as_str())
                }
                _ => None,
            };
            route.push(term);
        }
        for ((name, value), id) in meta.iter().zip(&route) {
            if exclude.contains(&name.as_str()) {
                continue;
            }
            if let Some(id) = id {
                if owner == "sample"
                    && matches!(
                        name.as_str(),
                        "GO cellular component" | "brenda source tissue"
                    )
                {
                    self.cv(id, None)?;
                } else {
                    self.cv(id, Some(value))?;
                }
            }
        }
        for ((name, value), id) in meta.iter().zip(route) {
            if !exclude.contains(&name.as_str()) && id.is_none() {
                self.user(name, value)?;
            }
        }
        Ok(())
    }
}
fn reversible_metadata(owner: &str, key: &str, id: &str) -> bool {
    match owner {
        "spectrum" => record_transport::SPECTRUM_CV
            .iter()
            .any(|row| row.0 == id && row.1 == key),
        "sample" => key == "sample batch" || id.starts_with("PATO:"),
        "instrumentConfiguration" => matches!(
            key,
            "instrument serial number"
                | "transmission"
                | "accelerating voltage"
                | "electric field strength"
                | "field-free region"
                | "space charge effect"
        ),
        "source" => matches!(
            key,
            "ionization efficiency"
                | "source potential"
                | "declustering potential"
                | "cone voltage"
                | "tube lens"
                | "wavelength"
                | "focus diameter x"
                | "focus diameter y"
                | "pulse energy"
                | "pulse duration"
                | "attenuation"
                | "impact angle"
                | "matrix solution"
                | "matrix concentration"
        ),
        _ => false,
    }
}
fn compatible_value(kind: XRefType, name: &str, value: &MetaValue) -> bool {
    if matches!(name, "field-free region" | "space charge effect") {
        return value.unit().is_none()
            && matches!(value.data(),MetaValueData::String(v) if v=="true");
    }
    match (kind, value.data()) {
        (XRefType::Decimal, MetaValueData::Float(_)) => true,
        (
            XRefType::Integer
            | XRefType::NegativeInteger
            | XRefType::PositiveInteger
            | XRefType::NonNegativeInteger
            | XRefType::NonPositiveInteger,
            MetaValueData::Integer(n),
        ) => i32::try_from(*n).is_ok(),
        (XRefType::String | XRefType::None | XRefType::AnyUri, MetaValueData::String(_)) => true,
        _ => false,
    }
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source record metadata, primary detector roles and independent noise grids.
use super::*;
use crate::metadata::MetaInfo;

pub(super) const COORDINATE: &str = "mzml coordinate array";
pub(super) const INTENSITY: &str = "mzml intensity array";
pub(super) const NOISE: [&str; 3] = [
    "sampled noise m/z array",
    "sampled noise intensity array",
    "sampled noise baseline array",
];
pub(super) const SPECTRUM_SKIP: [&str; 5] = [COORDINATE, INTENSITY, NOISE[0], NOISE[1], NOISE[2]];
pub(super) const CHROMATOGRAM_SKIP: [&str; 3] =
    [COORDINATE, INTENSITY, "chromatogram type accession"];

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Role {
    Wavelength,
    Absorption,
    Pressure,
    Flow,
    Detector,
    NoiseMz,
    NoiseIntensity,
    NoiseBaseline,
}
impl Role {
    pub(super) fn terms(self) -> (&'static str, &'static str, &'static str, &'static str) {
        match self {
            Self::Wavelength => ("MS:1000617", "wavelength array", "UO:0000018", "nanometer"),
            Self::Absorption => (
                "MS:1000515",
                "intensity array",
                "UO:0000269",
                "absorbance unit",
            ),
            Self::Pressure => (
                "MS:1000821",
                "pressure array",
                "UO:0000109",
                "pressure unit",
            ),
            Self::Flow => (
                "MS:1000820",
                "flow rate array",
                "UO:0000270",
                "volumetric flow rate unit",
            ),
            Self::Detector => (
                "MS:1000786",
                "non-standard data array",
                "UO:0000000",
                "unit",
            ),
            Self::NoiseMz => ("MS:1002743", NOISE[0], "MS:1000040", "m/z"),
            Self::NoiseIntensity => (
                "MS:1002744",
                NOISE[1],
                "MS:1000131",
                "number of detector counts",
            ),
            Self::NoiseBaseline => ("MS:1002745", NOISE[2], "", ""),
        }
    }
    pub(super) fn noise(self) -> bool {
        matches!(
            self,
            Self::NoiseMz | Self::NoiseIntensity | Self::NoiseBaseline
        )
    }
    fn selector(self) -> &'static str {
        match self {
            Self::Wavelength => "wavelength",
            Self::Absorption => "absorption",
            Self::Pressure => "pressure",
            Self::Flow => "flow",
            Self::Detector => "nonstandard",
            _ => unreachable!(),
        }
    }
}
pub(super) fn noise_role(index: usize) -> Role {
    [Role::NoiseMz, Role::NoiseIntensity, Role::NoiseBaseline][index]
}
pub(super) fn slot(budget: &mut ParameterBudget, key: &str) -> Result<()> {
    budget.spend(
        512usize
            .saturating_add(11 * (std::mem::size_of::<(String, MetaValue)>() + 32))
            .saturating_add(key.len()),
    )
}
pub(super) fn array_role(
    b: &mut Binary,
    attrs: &BTreeMap<String, String>,
    budget: &mut ParameterBudget,
) -> Result<Option<Kind>> {
    let id = required(attrs, "accession")?;
    let role = match id {
        "MS:1000617" if b.spectrum => Some(Role::Wavelength),
        "MS:1000515"
            if attrs
                .get("unitAccession")
                .is_some_and(|v| v == "UO:0000269") =>
        {
            Some(Role::Absorption)
        }
        "MS:1000821" if !b.spectrum => Some(Role::Pressure),
        "MS:1000820" if !b.spectrum => Some(Role::Flow),
        "MS:1000786"
            if !b.spectrum && attrs.get("value").is_some_and(|v| v == "detector signal") =>
        {
            Some(Role::Detector)
        }
        "MS:1002743" if b.spectrum => Some(Role::NoiseMz),
        "MS:1002744" if b.spectrum => Some(Role::NoiseIntensity),
        "MS:1002745" if b.spectrum => Some(Role::NoiseBaseline),
        _ => None,
    };
    let Some(role) = role else { return Ok(None) };
    let (_, _, unit, _) = role.terms();
    if attrs.get("unitAccession").is_some_and(|v| v != unit)
        || (unit.is_empty() && attrs.contains_key("unitAccession"))
    {
        return Err(Error::Unsupported(
            "conflicting primary/noise array unit".into(),
        ));
    }
    if (attrs.contains_key("unitName") || attrs.contains_key("unitCvRef"))
        && !attrs.contains_key("unitAccession")
    {
        return Err(invalid("array unit attributes require an accession"));
    }
    if let Some(prefix) = attrs.get("unitCvRef") {
        if unit.split_once(':').map(|p| p.0) != Some(prefix.as_str()) {
            return Err(invalid("array unit CV identity mismatch"));
        }
    }
    if !role.noise() {
        let key = if role == Role::Wavelength {
            COORDINATE
        } else {
            INTENSITY
        };
        slot(budget, key)?;
        if b.metadata
            .insert(key.into(), role.selector().into())
            .is_some()
        {
            return Err(invalid("duplicate array role metadata"));
        }
        if role != Role::Absorption {
            if let Some(unit) = attrs.get("unitAccession") {
                slot(budget, "unit_accession")?;
                if b.metadata
                    .insert("unit_accession".into(), unit.clone().into())
                    .is_some()
                {
                    return Err(invalid("duplicate array unit metadata"));
                }
            }
        }
    }
    Ok(Some(Kind::Role(role)))
}

// Exact source record routes. Types use XMLHandler's existing scalar decoder.
pub(super) const SPECTRUM_CV: &[(&str, &str, &str)] = &[
    ("MS:1000285", "total ion current", "xsd:double"),
    ("MS:1000504", "base peak m/z", "xsd:double"),
    ("MS:1000505", "base peak intensity", "xsd:double"),
    ("MS:1000527", "highest observed m/z", "xsd:double"),
    ("MS:1000528", "lowest observed m/z", "xsd:double"),
    ("MS:1000618", "highest observed wavelength", "xsd:double"),
    ("MS:1000619", "lowest observed wavelength", "xsd:double"),
    ("MS:1000796", "spectrum title", "xsd:string"),
    ("MS:1000797", "peak list scans", "xsd:string"),
    ("MS:1000798", "peak list raw scans", "xsd:string"),
];
const SCAN_CV: &[(&str, &str, &str)] = &[
    ("MS:1000502", "dwell time", "xsd:double"),
    ("MS:1000011", "mass resolution", "xsd:string"),
    ("MS:1000015", "scan rate", "xsd:double"),
    ("MS:1000512", "filter string", "xsd:string"),
    ("MS:1000803", "analyzer scan offset", "xsd:double"),
    ("MS:1000616", "preset scan configuration", "xsd:string"),
    ("MS:1000800", "mass resolving power", "xsd:string"),
    ("MS:1000880", "interchannel delay", "xsd:double"),
];
pub(super) fn read_cv(
    record: &mut Record,
    parent: &str,
    attrs: &BTreeMap<String, String>,
    budget: &mut ParameterBudget,
) -> Result<bool> {
    if precursor_metadata::read_metadata_cv(record, parent, attrs, budget)? {
        return Ok(true);
    }
    let id = required(attrs, "accession")?;
    let row = match parent {
        "spectrum" => SPECTRUM_CV,
        "scan" => SCAN_CV,
        _ => &[],
    }
    .iter()
    .find(|r| r.0 == id);
    let (key, value) = if let Some(&(_, key, kind)) = row {
        slot(budget, key)?;
        (key, scalar_user_value(attrs, kind)?)
    } else if parent == "scan" && id == "MS:1000826" {
        slot(budget, "elution time (seconds)")?;
        let value = finite(required(attrs, "value")?, "elution time")?
            * if attrs
                .get("unitAccession")
                .is_some_and(|v| v == "UO:0000031")
            {
                60.0
            } else {
                1.0
            };
        ("elution time (seconds)", MetaValue::try_from(value)?)
    } else if parent == "scan"
        && matches!(
            id,
            "MS:1000092" | "MS:1000093" | "MS:1000094" | "MS:1000095" | "MS:1000096"
        )
    {
        let (key, text) = match id {
            "MS:1000092" => ("scan direction", "decreasing"),
            "MS:1000093" => ("scan direction", "increasing"),
            "MS:1000094" => ("scan law", "exponential"),
            "MS:1000095" => ("scan law", "linear"),
            _ => ("scan law", "quadratic"),
        };
        slot(budget, key)?;
        (key, text.into())
    } else {
        return Ok(false);
    };
    if record.metadata().insert(key.into(), value).is_some() {
        return Err(invalid("duplicate routed record metadata"));
    }
    Ok(true)
}
pub(super) fn noise_values(meta: &MetaInfo, index: usize) -> Result<Option<&[f64]>> {
    meta.get(NOISE[index])
        .map(|v| {
            if v.unit().is_some() {
                return Err(Error::Unsupported(
                    "sampled noise list cannot have metadata units".into(),
                ));
            }
            v.as_float_list()
        })
        .transpose()
}
pub(super) fn primary_kind(meta: &MetaInfo, chrom: bool, coordinate: bool) -> Result<Kind> {
    let key = if coordinate { COORDINATE } else { INTENSITY };
    let fallback = if coordinate {
        if chrom { Kind::Time } else { Kind::Mz }
    } else {
        Kind::Intensity
    };
    let Some(value) = meta.get(key) else {
        return Ok(fallback);
    };
    if value.unit().is_some() {
        return Err(invalid("array selector cannot have a unit"));
    }
    let kind = match (coordinate, chrom, value.as_str()?) {
        (true, false, "wavelength") => Kind::Role(Role::Wavelength),
        (false, _, "absorption") => Kind::Role(Role::Absorption),
        (false, true, "pressure") => Kind::Role(Role::Pressure),
        (false, true, "flow") => Kind::Role(Role::Flow),
        (false, true, "nonstandard") => Kind::Role(Role::Detector),
        _ => {
            return Err(Error::Unsupported(
                "unsupported or redundant primary array selector".into(),
            ));
        }
    };
    Ok(kind)
}
pub(super) fn validate(meta: &MetaInfo, chrom: bool) -> Result<()> {
    let coordinate = primary_kind(meta, chrom, true)?;
    let intensity = primary_kind(meta, chrom, false)?;
    if let Some(unit) = meta.get("unit_accession") {
        if matches!(coordinate, Kind::Role(_)) || matches!(intensity, Kind::Role(_)) {
            if unit.unit().is_some() {
                return Err(invalid("array unit selector cannot have a unit"));
            }
            let role_unit = match (&coordinate, &intensity) {
                (Kind::Role(r), _) if *r != Role::Absorption => Some(r.terms().2),
                (_, Kind::Role(r)) if *r != Role::Absorption => Some(r.terms().2),
                _ => None,
            };
            if role_unit.is_some_and(|u| unit.as_str().ok() != Some(u)) {
                return Err(Error::Unsupported(
                    "retained primary unit conflicts with source output unit".into(),
                ));
            }
        }
    }
    for (key, value) in meta {
        if !chrom && NOISE.contains(&key.as_str()) {
            noise_values(meta, NOISE.iter().position(|n| n == key).unwrap())?;
            continue;
        }
        if key == NAME_KEY {
            return Err(invalid("reserved record name metadata"));
        }
        validate_scalar_value(key, value)?;
    }
    Ok(())
}

// A role promoted by this reader must have exactly one native owner.
pub(super) fn array_owners(
    meta: &MetaInfo,
    chrom: bool,
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
) -> Result<()> {
    for name in floats
        .iter()
        .map(|a| a.name.as_str())
        .chain(integers.iter().map(|a| a.name.as_str()))
        .chain(strings.iter().map(|a| a.name.as_str()))
    {
        if (!chrom && NOISE.contains(&name))
            || (chrom
                && matches!(
                    name,
                    "pressure array" | "flow rate array" | "detector signal"
                ))
            || (!chrom
                && name == "wavelength array"
                && matches!(
                    primary_kind(meta, false, true)?,
                    Kind::Role(Role::Wavelength)
                ))
        {
            return Err(Error::Unsupported(
                "promoted array role requires record metadata ownership".into(),
            ));
        }
    }
    Ok(())
}

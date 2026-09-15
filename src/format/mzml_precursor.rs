// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Typed precursor acquisition CV terms from the pinned MzMLHandler.
//!
//! Activation and intensity-unit transport: `docs/MZML_PRECURSOR_ACTIVATION_SUPPORT.md`.

use super::{
    BTreeMap, Error, Precursor, Result, Write, cv, escape, finite, invalid, number, xml_string,
};
use crate::metadata::{ActivationMethod as A, DriftTimeUnit as D};

const METHODS: [(A, &str); 19] = [
    (A::Cid, "MS:1000133"),
    (A::Psd, "MS:1000135"),
    (A::Pd, "MS:1000134"),
    (A::Sid, "MS:1000136"),
    (A::Bird, "MS:1000242"),
    (A::Ecd, "MS:1000250"),
    (A::Imd, "MS:1000262"),
    (A::Sori, "MS:1000282"),
    (A::Hcid, "MS:1002481"),
    (A::Lcid, "MS:1000433"),
    (A::Phd, "MS:1000435"),
    (A::Etd, "MS:1000598"),
    (A::Etcid, "MS:1003182"),
    (A::Ethcd, "MS:1002631"),
    (A::Pqd, "MS:1000599"),
    (A::Trap, "MS:1002472"),
    (A::Hcd, "MS:1000422"),
    (A::InSource, "MS:1001880"),
    (A::Lift, "MS:1002000"),
];
// Only explicit pinned activation routes are retained here; unknown CVs keep
// the established reader policy. The source calls MS:1000138 metadata by its
// historical key even though its current CV name is normalized collision energy.
const ACTIVATION_METADATA: [(&str, &str, &str); 9] = [
    ("MS:1000245", "charge stripping", "xsd:string"),
    ("MS:1000045", "collision energy", "xsd:double"),
    ("MS:1000412", "buffer gas", "xsd:string"),
    ("MS:1000419", "collision gas", "xsd:string"),
    ("MS:1000138", "percent collision energy", "xsd:double"),
    ("MS:1000869", "collision gas pressure", "xsd:double"),
    (
        "MS:1002679",
        "supplemental collision-induced dissociation",
        "xsd:string",
    ),
    (
        "MS:1002678",
        "supplemental beam-type collision-induced dissociation",
        "xsd:string",
    ),
    ("MS:1002680", "supplemental collision energy", "xsd:double"),
];
const INTENSITY_UNIT_KEY: &str = "peak intensity unit accession";

pub(super) fn read_metadata_cv(
    record: &mut super::Record,
    parent: &str,
    attrs: &BTreeMap<String, String>,
    budget: &mut super::ParameterBudget,
) -> Result<bool> {
    if record.precursor.is_none() {
        return Ok(false);
    }
    let id = super::required(attrs, "accession")?;
    if parent == "selectedIon" && id == "MS:1000042" {
        if !record.precursor_fields.insert("precursor_intensity") {
            return Err(invalid("duplicate precursor intensity"));
        }
        let p = record.precursor.as_mut().expect("checked precursor");
        p.intensity = super::intensity(finite(
            super::required(attrs, "value")?,
            "precursor intensity",
        )?)?;
        if let Some(unit) = attrs.get("unitAccession") {
            record.precursor_fields.insert("intensity_explicit_unit");
            super::record_transport::slot(budget, INTENSITY_UNIT_KEY)?;
            intensity_unit(unit)?;
            if attrs
                .get("unitCvRef")
                .is_some_and(|prefix| !unit.starts_with(&format!("{prefix}:")))
            {
                return Err(invalid(
                    "precursor intensity unit reference conflicts with accession",
                ));
            }
            if p.cv_terms.metadata.contains_key(INTENSITY_UNIT_KEY) {
                return Err(invalid("duplicate precursor intensity unit metadata"));
            }
            if unit != "MS:1000132" {
                p.cv_terms
                    .metadata
                    .insert(INTENSITY_UNIT_KEY.into(), unit.as_str().into());
            }
        } else if attrs.contains_key("unitCvRef") || attrs.contains_key("unitName") {
            return Err(invalid(
                "precursor intensity unit attributes require an accession",
            ));
        }
        return Ok(true);
    }
    if parent != "activation" {
        return Ok(false);
    }
    let Some(&(_, key, kind)) = ACTIVATION_METADATA.iter().find(|row| row.0 == id) else {
        return Ok(false);
    };
    super::record_transport::slot(budget, key)?;
    let value = if id == "MS:1000245" {
        if attrs.contains_key("unitAccession")
            || attrs.contains_key("unitCvRef")
            || attrs.contains_key("unitName")
        {
            return Err(invalid(
                "charge stripping flag cannot preserve unit metadata",
            ));
        }
        crate::metadata::MetaValue::from("true")
    } else {
        super::scalar_user_value(attrs, kind)?
    };
    let p = record.precursor.as_mut().expect("checked precursor");
    if p.cv_terms.metadata.insert(key.into(), value).is_some() {
        return Err(invalid("duplicate precursor activation metadata"));
    }
    match id {
        "MS:1002679" => {
            p.activation_methods.insert(A::Etcid);
        }
        "MS:1002678" => {
            p.activation_methods.insert(A::Ethcd);
        }
        _ => {}
    }
    Ok(true)
}

fn intensity_unit(
    accession: &str,
) -> Result<&'static crate::format::controlled_vocabulary::CVTermDefinition> {
    if !accession.starts_with("MS:") && !accession.starts_with("UO:") {
        return Err(Error::Unsupported(
            "precursor intensity unit must have MS or UO identity".into(),
        ));
    }
    crate::format::controlled_vocabulary::ControlledVocabulary::psi_ms()?.get_term(accession)
}

/// Write the selected ion's `MS:1000042` peak intensity term, or nothing.
///
/// Source `MzMLHandler.cpp:4596-4601` writes the term only for a positive
/// intensity (or a `peak intensity` meta value), always with unit attributes,
/// `MS:1000132` unless `peak intensity unit accession` names another unit. A
/// precursor without an intensity reads back as `0.0`, so the default `+0.0`
/// is omitted here as well; the benchmark slice showed this port adding
/// `value="0"` to 49 of 600 spectra that neither the input nor the C++ output
/// carries. Two native differences keep the term where the source drops it: a
/// negative or `-0.0` intensity, which would otherwise read back changed, and
/// an explicit non-default unit, whose identity would otherwise be lost.
pub(super) fn write_intensity(w: &mut impl Write, p: &Precursor) -> Result<()> {
    let explicit_unit = p.cv_terms.metadata.get(INTENSITY_UNIT_KEY);
    if explicit_unit.is_none() && p.intensity.to_bits() == 0 {
        return Ok(());
    }
    let accession = match explicit_unit {
        Some(value) => value.as_str()?,
        None => "MS:1000132",
    };
    let term = intensity_unit(accession)?;
    let prefix = term.id.split_once(':').expect("validated unit prefix").0;
    let unit = format!(
        " unitAccession=\"{}\" unitCvRef=\"{}\" unitName=\"{}\"",
        escape(&term.id),
        escape(prefix),
        escape(&term.name)
    );
    cv(
        w,
        "MS:1000042",
        "peak intensity",
        &p.intensity.to_string(),
        &unit,
    )
}

fn promoted(p: &Precursor, value: &crate::metadata::MetaValue, id: &str, kind: &str) -> bool {
    // These CVs derive a combined activation method on read. Keep a metadata-only
    // caller value as userParam instead of silently adding a method on reload.
    if (id == "MS:1002679" && !p.activation_methods.contains(&A::Etcid))
        || (id == "MS:1002678" && !p.activation_methods.contains(&A::Ethcd))
    {
        return false;
    }
    use crate::metadata::MetaValueData;
    match (kind, value.data()) {
        ("xsd:double", MetaValueData::Float(_)) => true,
        ("xsd:string", MetaValueData::String(text)) if id == "MS:1000245" => {
            text == "true" && value.unit().is_none()
        }
        ("xsd:string", MetaValueData::String(_)) => true,
        _ => false,
    }
}

const MZ_UNIT: &str = " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"";

pub(super) fn field(parent: &str, accession: &str) -> Option<&'static str> {
    match (parent, accession) {
        ("isolationWindow", "MS:1000828") => Some("isolation_lower"),
        ("isolationWindow", "MS:1000829") => Some("isolation_upper"),
        ("activation", "MS:1000509") => Some("activation_energy"),
        ("selectedIon", "MS:1002476" | "MS:1002815" | "MS:1001581" | "MS:1002954") => {
            Some("precursor_mobility")
        }
        _ => None,
    }
}

/// Drift-time unit and PSI-MS/UO unit accession of a source mobility term.
///
/// `MzMLHandler.cpp` maps the same four accessions on the selected-ion
/// (1839-1882), scan (2279-2311) and spectrum (1731-1738, FAIMS only) routes.
/// Every Rust route reads this one table, so the accession/unit pairs cannot
/// drift apart between them.
pub(super) fn mobility_term(accession: &str) -> Option<(D, &'static str)> {
    match accession {
        "MS:1002476" => Some((D::Millisecond, "UO:0000028")),
        "MS:1002815" => Some((D::InverseReducedMobility, "MS:1002814")),
        "MS:1001581" => Some((D::FaimsCompensationVoltage, "UO:0000218")),
        "MS:1002954" => Some((D::CollisionCrossSection, "UO:0000324")),
        _ => None,
    }
}

/// Accession, CV name and unit attributes the writers emit for a mobility
/// unit, or `None` for `DriftTimeUnit::None`, which has no mzML term.
///
/// Shared by the selected-ion writer and the scan writer
/// (`MzMLHandler.cpp:4607-4628` and 5412-5440).
pub(super) fn mobility_cv(unit: D) -> Option<(&'static str, &'static str, &'static str)> {
    match unit {
        D::Millisecond => Some((
            "MS:1002476",
            "ion mobility drift time",
            " unitCvRef=\"UO\" unitAccession=\"UO:0000028\" unitName=\"millisecond\"",
        )),
        D::InverseReducedMobility => Some((
            "MS:1002815",
            "inverse reduced ion mobility",
            " unitCvRef=\"MS\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\"",
        )),
        D::FaimsCompensationVoltage => Some((
            "MS:1001581",
            "FAIMS compensation voltage",
            " unitCvRef=\"UO\" unitAccession=\"UO:0000218\" unitName=\"volt\"",
        )),
        D::CollisionCrossSection => Some((
            "MS:1002954",
            "collisional cross sectional area",
            " unitCvRef=\"UO\" unitAccession=\"UO:0000324\" unitName=\"square angstrom\"",
        )),
        D::None => None,
    }
}
fn unit(attrs: &BTreeMap<String, String>, expected: &str) -> Result<()> {
    if attrs.get("unitAccession").is_some_and(|v| v != expected) {
        return Err(Error::Unsupported(
            "precursor CV unit does not match its typed quantity".into(),
        ));
    }
    Ok(())
}
pub(super) fn read_cv(
    p: &mut Precursor,
    parent: &str,
    accession: &str,
    value: &str,
    attrs: &BTreeMap<String, String>,
    selected_mz: bool,
) -> Result<bool> {
    match (parent, accession) {
        ("isolationWindow", "MS:1000827" | "MS:1000828" | "MS:1000829") => {
            unit(attrs, "MS:1000040")?;
            let v = finite(value, "isolation window quantity")?;
            if v < 0.0 {
                return Err(invalid("negative isolation window quantity"));
            }
            match accession {
                "MS:1000827" => {
                    p.isolation_target_mz = Some(v);
                    if !selected_mz {
                        p.mz = v;
                    }
                }
                "MS:1000828" => p.isolation_window_lower_offset = v,
                _ => p.isolation_window_upper_offset = v,
            }
        }
        ("selectedIon", "MS:1000633") => p
            .possible_charge_states
            .push(number(value, "possible charge state")?),
        ("selectedIon", "MS:1002476" | "MS:1002815" | "MS:1001581" | "MS:1002954") => {
            let Some((kind, cv_unit)) = mobility_term(accession) else {
                return Ok(false);
            };
            unit(attrs, cv_unit)?;
            p.drift_time = Some(finite(value, "precursor mobility")?);
            p.drift_time_unit = kind;
        }
        ("activation", "MS:1000509") => {
            unit(attrs, "UO:0000266")?;
            p.activation_energy = finite(value, "activation energy")?;
        }
        ("activation", _) => {
            let method = METHODS
                .iter()
                .find(|(_, a)| *a == accession)
                .map(|(m, _)| *m)
                .or(match accession {
                    "MS:1002679" => Some(A::Etcid),
                    "MS:1002678" => Some(A::Ethcd),
                    _ => None,
                });
            if let Some(method) = method {
                p.activation_methods.insert(method);
            } else {
                return Ok(false);
            }
        }
        _ => return Ok(false),
    }
    Ok(true)
}

pub(super) fn validate_write(p: &Precursor) -> Result<()> {
    p.validate()?;
    selected_mz(p)?;
    if (p.isolation_target_mz.is_some()
        || p.isolation_window_lower_offset != 0.0
        || p.isolation_window_upper_offset != 0.0)
        && p.isolation_target_mz.unwrap_or(p.mz) < 0.0
    {
        return Err(invalid("negative effective isolation target"));
    }
    if let Some(reference) = &p.spectrum_reference {
        xml_string(reference)?;
    }
    if p.drift_window_lower_offset != 0.0
        || p.drift_window_upper_offset != 0.0
        || !p.cv_terms.is_empty()
    {
        return Err(Error::Unsupported(
            "mzML writer cannot store precursor drift-window offsets or arbitrary CV metadata"
                .into(),
        ));
    }
    if (p.drift_time.is_some() && p.drift_time_unit == D::None)
        || (p.drift_time.is_none() && p.drift_time_unit != D::None)
    {
        return Err(Error::Unsupported(
            "precursor mobility requires both a value and an explicit unit".into(),
        ));
    }
    super::validate_scalar_metadata(&p.cv_terms.metadata)?;
    if let Some(value) = p.cv_terms.metadata.get(INTENSITY_UNIT_KEY) {
        if value.unit().is_some() {
            return Err(invalid(
                "intensity unit accession metadata cannot itself have a unit",
            ));
        }
        let accession = value.as_str()?;
        intensity_unit(accession)?;
        if accession == "MS:1000132" {
            return Err(invalid(
                "explicit default intensity unit metadata would normalize away",
            ));
        }
    }
    if let Some(value) = p.cv_terms.metadata.get("external_spectrum_id") {
        if value.unit().is_some() {
            return Err(invalid("external_spectrum_id cannot have a unit"));
        }
        value.as_str()?;
    }
    if p.cv_terms
        .metadata
        .contains_key("activation information unavailable")
    {
        return Err(invalid("reserved precursor fallback metadata name"));
    }
    Ok(())
}

/// Source-owned discriminator: emitted in selectedIon, never as activation metadata.
pub(super) fn selected_mz(p: &Precursor) -> Result<f64> {
    let Some(value) = p.cv_terms.metadata.get("selected ion m/z") else {
        return Ok(p.mz);
    };
    if value.unit().is_some() {
        return Err(invalid("selected ion m/z metadata cannot have a unit"));
    }
    let mz = value.as_f64()?;
    if !mz.is_finite() || mz < 0.0 {
        return Err(invalid("invalid selected ion m/z metadata"));
    }
    Ok(mz)
}

pub(super) fn write_start(w: &mut impl Write, p: &Precursor, tpp: bool) -> Result<()> {
    write!(w, "<precursor")?;
    if let Some(reference) = &p.spectrum_reference {
        write!(w, " spectrumRef=\"{}\"", escape(reference))?;
    }
    if let Some(value) = p.cv_terms.metadata.get("external_spectrum_id") {
        write!(w, " externalSpectrumID=\"{}\"", escape(value.as_str()?))?;
    }
    writeln!(w, ">")?;
    if !tpp
        && (p.isolation_target_mz.is_some()
            || p.cv_terms.metadata.contains_key("selected ion m/z")
            || p.isolation_window_lower_offset != 0.0
            || p.isolation_window_upper_offset != 0.0)
    {
        writeln!(w, "<isolationWindow>")?;
        for (accession, name, value) in [
            (
                "MS:1000827",
                "isolation window target m/z",
                p.isolation_target_mz.unwrap_or(p.mz),
            ),
            (
                "MS:1000828",
                "isolation window lower offset",
                p.isolation_window_lower_offset,
            ),
            (
                "MS:1000829",
                "isolation window upper offset",
                p.isolation_window_upper_offset,
            ),
        ] {
            cv(w, accession, name, &value.to_string(), MZ_UNIT)?;
        }
        writeln!(w, "</isolationWindow>")?;
    }
    writeln!(w, "<selectedIonList count=\"1\"><selectedIon>")?;
    Ok(())
}
pub(super) fn write_end(w: &mut impl Write, p: &Precursor) -> Result<()> {
    for charge in &p.possible_charge_states {
        cv(
            w,
            "MS:1000633",
            "possible charge state",
            &charge.to_string(),
            "",
        )?;
    }
    if let Some(value) = p.drift_time {
        let Some((accession, name, unit)) = mobility_cv(p.drift_time_unit) else {
            return Err(invalid("mobility unit missing after preflight"));
        };
        cv(w, accession, name, &value.to_string(), unit)?;
    }
    writeln!(w, "</selectedIon></selectedIonList><activation>")?;
    if p.activation_energy != 0.0 {
        cv(
            w,
            "MS:1000509",
            "activation energy",
            &p.activation_energy.to_string(),
            " unitCvRef=\"UO\" unitAccession=\"UO:0000266\" unitName=\"electronvolt\"",
        )?;
    }
    for method in &p.activation_methods {
        if (*method == A::Etcid
            && p.cv_terms
                .metadata
                .get("supplemental collision-induced dissociation")
                .is_some_and(|v| promoted(p, v, "MS:1002679", "xsd:string")))
            || (*method == A::Ethcd
                && p.cv_terms
                    .metadata
                    .get("supplemental beam-type collision-induced dissociation")
                    .is_some_and(|v| promoted(p, v, "MS:1002678", "xsd:string")))
        {
            continue;
        }
        let accession = METHODS.iter().find(|(m, _)| m == method).unwrap().1;
        let name = if *method == A::Lift {
            "LIFT".into()
        } else {
            method.name().to_lowercase()
        };
        cv(w, accession, &name, "", "")?;
    }
    let mut skip = [""; 12];
    skip[..3].copy_from_slice(&[
        "external_spectrum_id",
        "selected ion m/z",
        INTENSITY_UNIT_KEY,
    ]);
    let mut used = 3;
    for (id, key, kind) in ACTIVATION_METADATA {
        if let Some(value) = p
            .cv_terms
            .metadata
            .get(key)
            .filter(|v| promoted(p, v, id, kind))
        {
            let term = crate::format::controlled_vocabulary::ControlledVocabulary::psi_ms()?
                .get_term(id)?;
            let text = match value.data() {
                _ if id == "MS:1000245" => String::new(),
                // The source renders the promoted meta value with
                // `DataValue::toString`; see `mzml::float_text`.
                crate::metadata::MetaValueData::Float(number) => super::float_text(*number),
                _ => value.to_string(),
            };
            let unit = value
                .unit()
                .map(|u| {
                    format!(
                        " unitAccession=\"{}\" unitCvRef=\"{}\" unitName=\"{}\"",
                        escape(u.accession()),
                        escape(u.cv_ref()),
                        escape(u.name())
                    )
                })
                .unwrap_or_default();
            cv(w, id, &term.name, &text, &unit)?;
            skip[used] = key;
            used += 1;
        }
    }
    // After every cvParam: `activation` is a ParamGroupType, whose schema
    // sequence puts all cvParams before any userParam.
    if p.activation_methods.is_empty() && p.activation_energy == 0.0 {
        writeln!(
            w,
            "<userParam name=\"activation information unavailable\"/>"
        )?;
    }
    super::write_scalar_metadata_skipping(w, &p.cv_terms.metadata, &skip[..used])?;
    writeln!(w, "</activation></precursor>")?;
    Ok(())
}

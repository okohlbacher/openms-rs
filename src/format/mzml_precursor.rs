// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Typed precursor acquisition CV terms from the pinned MzMLHandler.

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
            let (kind, cv_unit) = match accession {
                "MS:1002476" => (D::Millisecond, "UO:0000028"),
                "MS:1002815" => (D::InverseReducedMobility, "MS:1002814"),
                "MS:1001581" => (D::FaimsCompensationVoltage, "UO:0000218"),
                _ => (D::CollisionCrossSection, "UO:0000324"),
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
        || !p.cv_terms.metadata.is_empty()
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
    Ok(())
}

pub(super) fn write_start(w: &mut impl Write, p: &Precursor) -> Result<()> {
    if let Some(reference) = &p.spectrum_reference {
        writeln!(w, "<precursor spectrumRef=\"{}\">", escape(reference))?;
    } else {
        writeln!(w, "<precursor>")?;
    }
    if p.isolation_target_mz.is_some()
        || p.isolation_window_lower_offset != 0.0
        || p.isolation_window_upper_offset != 0.0
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
        let (accession, name, unit) = match p.drift_time_unit {
            D::Millisecond => (
                "MS:1002476",
                "ion mobility drift time",
                " unitCvRef=\"UO\" unitAccession=\"UO:0000028\" unitName=\"millisecond\"",
            ),
            D::InverseReducedMobility => (
                "MS:1002815",
                "inverse reduced ion mobility",
                " unitCvRef=\"MS\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\"",
            ),
            D::FaimsCompensationVoltage => (
                "MS:1001581",
                "FAIMS compensation voltage",
                " unitCvRef=\"UO\" unitAccession=\"UO:0000218\" unitName=\"volt\"",
            ),
            D::CollisionCrossSection => (
                "MS:1002954",
                "collisional cross sectional area",
                " unitCvRef=\"UO\" unitAccession=\"UO:0000324\" unitName=\"square angstrom\"",
            ),
            D::None => return Err(invalid("mobility unit missing after preflight")),
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
        let accession = METHODS.iter().find(|(m, _)| m == method).unwrap().1;
        let name = if *method == A::Lift {
            "LIFT".into()
        } else {
            method.name().to_lowercase()
        };
        cv(w, accession, &name, "", "")?;
    }
    if p.activation_methods.is_empty() && p.activation_energy == 0.0 {
        writeln!(
            w,
            "<userParam name=\"activation information unavailable\"/>"
        )?;
    }
    writeln!(w, "</activation></precursor>")?;
    Ok(())
}

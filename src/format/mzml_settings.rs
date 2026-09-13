// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Record-local mzML acquisition settings, without an instrument registry.

use super::*;
use crate::metadata::{ChromatogramType, InstrumentSettings, Polarity, ScanMode};

const SCAN_MODES: &[(ScanMode, &str, &str)] = &[
    (ScanMode::MassSpectrum, "MS:1000294", "mass spectrum"),
    (ScanMode::Ms1Spectrum, "MS:1000579", "MS1 spectrum"),
    (ScanMode::MsnSpectrum, "MS:1000580", "MSn spectrum"),
    (
        ScanMode::ConsecutiveReactionMonitoring,
        "MS:1000581",
        "CRM spectrum",
    ),
    (
        ScanMode::SelectedIonMonitoring,
        "MS:1000582",
        "SIM spectrum",
    ),
    (
        ScanMode::SelectedReactionMonitoring,
        "MS:1000583",
        "SRM spectrum",
    ),
    (
        ScanMode::ElectromagneticRadiation,
        "MS:1000804",
        "electromagnetic radiation spectrum",
    ),
    (ScanMode::Emission, "MS:1000805", "emission spectrum"),
    (ScanMode::Absorption, "MS:1000806", "absorption spectrum"),
    (
        ScanMode::ConstantNeutralGain,
        "MS:1000325",
        "constant neutral gain spectrum",
    ),
    (
        ScanMode::ConstantNeutralLoss,
        "MS:1000326",
        "constant neutral loss spectrum",
    ),
    (ScanMode::Precursor, "MS:1000341", "precursor ion spectrum"),
    (
        ScanMode::EnhancedMultiplyCharged,
        "MS:1000789",
        "enhanced multiply charged spectrum",
    ),
    (
        ScanMode::TimeDelayedFragmentation,
        "MS:1000790",
        "time-delayed fragmentation spectrum",
    ),
];
const CHROMATOGRAM_TYPES: &[(ChromatogramType, &str, &str)] = &[
    (
        ChromatogramType::Mass,
        "MS:1000810",
        "ion current chromatogram",
    ),
    (
        ChromatogramType::TotalIonCurrent,
        "MS:1000235",
        "total ion current chromatogram",
    ),
    (
        ChromatogramType::SelectedIonCurrent,
        "MS:1000627",
        "selected ion current chromatogram",
    ),
    (
        ChromatogramType::BasePeak,
        "MS:1000628",
        "basepeak chromatogram",
    ),
    (
        ChromatogramType::SelectedIonMonitoring,
        "MS:1001472",
        "selected ion monitoring chromatogram",
    ),
    (
        ChromatogramType::SelectedReactionMonitoring,
        "MS:1001473",
        "selected reaction monitoring chromatogram",
    ),
    (
        ChromatogramType::ElectromagneticRadiation,
        "MS:1000811",
        "electromagnetic radiation chromatogram",
    ),
    (
        ChromatogramType::Absorption,
        "MS:1000812",
        "absorption chromatogram",
    ),
    (
        ChromatogramType::Emission,
        "MS:1000813",
        "emission chromatogram",
    ),
];

pub(super) fn read_cv(
    record: &mut Record,
    parent: &str,
    attrs: &BTreeMap<String, String>,
) -> Result<bool> {
    let accession = required(attrs, "accession")?;
    if parent == "spectrum" {
        if let Some(&(mode, _, _)) = SCAN_MODES.iter().find(|t| t.1 == accession) {
            if !record.seen_fields.insert("scan_mode") {
                return Err(invalid("duplicate scan mode"));
            }
            record
                .spectrum
                .as_mut()
                .ok_or_else(|| invalid("scan mode outside spectrum"))?
                .instrument_settings
                .scan_mode = mode;
            return Ok(true);
        }
        if matches!(accession, "MS:1000129" | "MS:1000130") {
            if !record.seen_fields.insert("polarity") {
                return Err(invalid("duplicate scan polarity"));
            }
            record
                .spectrum
                .as_mut()
                .ok_or_else(|| invalid("polarity outside spectrum"))?
                .instrument_settings
                .polarity = if accession == "MS:1000129" {
                Polarity::Negative
            } else {
                Polarity::Positive
            };
            return Ok(true);
        }
    }
    if matches!(parent, "spectrum" | "scan") && accession == "MS:1000497" {
        record
            .spectrum
            .as_mut()
            .ok_or_else(|| invalid("zoom outside spectrum"))?
            .instrument_settings
            .zoom_scan = true;
        return Ok(true);
    }
    if parent == "chromatogram" {
        let mode = if accession == "MS:1001474" {
            Some(ChromatogramType::SelectedReactionMonitoring)
        } else {
            CHROMATOGRAM_TYPES
                .iter()
                .find(|t| t.1 == accession)
                .map(|t| t.0)
        };
        if let Some(mode) = mode {
            if !record.seen_fields.insert("chromatogram_type") {
                return Err(invalid("duplicate chromatogram type"));
            }
            record
                .chromatogram
                .as_mut()
                .ok_or_else(|| invalid("chromatogram type outside chromatogram"))?
                .chromatogram_type = mode;
            return Ok(true);
        }
        if matches!(accession, "MS:1003019" | "MS:1003020" | "MS:1000626") {
            let c = record
                .chromatogram
                .as_mut()
                .ok_or_else(|| invalid("chromatogram type outside record"))?;
            if c.metadata
                .insert("chromatogram type accession".into(), accession.into())
                .is_some()
            {
                return Err(invalid("duplicate chromatogram type override"));
            }
            return Ok(true);
        }
    }
    if parent == "scanWindow" {
        let field = match accession {
            "MS:1000501" => "lower",
            "MS:1000500" => "upper",
            // Unknown CV terms retain the existing reader's ignored-metadata
            // policy; only the two represented quantities are interpreted.
            _ => return Ok(false),
        };
        if !record.scan_window_fields.insert(field) {
            return Err(invalid("duplicate scan window bound"));
        }
        let unit = attrs
            .get("unitAccession")
            .map(String::as_str)
            .unwrap_or("MS:1000040");
        unit_identity(unit)?;
        if let Some(prefix) = attrs.get("unitCvRef") {
            if unit.split_once(':').map(|p| p.0) != Some(prefix.as_str()) {
                return Err(invalid("scan window unit prefix mismatch"));
            }
        }
        if record
            .scan_window_unit
            .as_deref()
            .is_some_and(|prior| prior != unit)
        {
            return Err(Error::Unsupported(
                "mixed scan window endpoint units".into(),
            ));
        }
        record.scan_window_unit = Some(unit.into());
        let window = record
            .scan_window
            .as_mut()
            .ok_or_else(|| invalid("bound outside scan window"))?;
        if unit != "MS:1000040" {
            window.metadata.insert("unit_accession".into(), unit.into());
        }
        let value = finite(
            attrs.get("value").map(String::as_str).unwrap_or(""),
            "scan window bound",
        )?;
        if field == "lower" {
            window.begin = value;
        } else {
            window.end = value;
        }
        return Ok(true);
    }
    Ok(false)
}

// mzML permits omitting unitName. Keep opaque MS/UO identities rather than
// pretending to resolve an ontology that is not part of this adapter.
fn unit_identity(accession: &str) -> Result<(&str, Option<&'static str>)> {
    let Some((prefix, digits)) = accession.split_once(':') else {
        return Err(invalid("invalid scan window unit accession"));
    };
    if !matches!(prefix, "MS" | "UO")
        || digits.len() != 7
        || !digits.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(Error::Unsupported(
            "scan window units require an MS/UO numeric accession".into(),
        ));
    }
    Ok((
        prefix,
        match accession {
            "MS:1000040" => Some("m/z"),
            "UO:0000018" => Some("nanometer"),
            _ => None,
        },
    ))
}

fn instrument_unrepresented(settings: &InstrumentSettings) -> bool {
    settings.scan_mode != ScanMode::Unknown
        || settings.polarity != Polarity::Unknown
        || settings.zoom_scan
        || !settings.scan_windows.is_empty()
        || !settings.metadata.is_empty()
}
pub(super) fn spectrum_guard(s: &MSSpectrum) -> Result<()> {
    if !s.instrument_settings.metadata.is_empty() {
        return Err(Error::Unsupported("mzML spectrum SourceFile, DataProcessing or InstrumentSettings metadata are not represented".into()));
    }
    Ok(())
}
pub(super) fn chromatogram_guard(c: &MSChromatogram) -> Result<()> {
    if !c.acquisition_info.acquisitions.is_empty()
        || !c.acquisition_info.metadata.is_empty()
        || !c.acquisition_info.method_of_combination.is_empty()
        || instrument_unrepresented(&c.instrument_settings)
    {
        return Err(Error::Unsupported("mzML chromatogram instrument/acquisition/source/processing settings are not represented".into()));
    }
    if let Some(value) = c.metadata.get("chromatogram type accession") {
        if value.unit().is_some()
            || !matches!(value.as_str()?, "MS:1003019" | "MS:1003020" | "MS:1000626")
            || c.chromatogram_type != ChromatogramType::Mass
        {
            return Err(Error::Unsupported(
                "chromatogram type override conflicts with typed settings".into(),
            ));
        }
    }
    if c.chromatogram_type == ChromatogramType::Unknown {
        return Err(Error::Unsupported(
            "unknown chromatogram type cannot roundtrip through mzML".into(),
        ));
    }
    Ok(())
}
pub(super) fn validate_spectrum(s: &MSSpectrum) -> Result<()> {
    acquisition_metadata::validate(&s.acquisition_info)?;
    for window in &s.instrument_settings.scan_windows {
        validate_scalar_metadata(&window.metadata)?;
        if let Some(value) = window.metadata.get("unit_accession") {
            if value.unit().is_some() {
                return Err(Error::Unsupported(
                    "unit_accession metadata cannot itself have a unit".into(),
                ));
            }
            unit_identity(value.as_str()?)?;
            if value.as_str()? == "MS:1000040" {
                return Err(Error::Unsupported("explicit default unit_accession metadata cannot roundtrip through source mzML semantics".into()));
            }
        }
    }
    for product in &s.products {
        validate_product_write(product)?;
    }
    Ok(())
}
pub(super) fn write_spectrum(w: &mut impl Write, spectrum: &MSSpectrum) -> Result<()> {
    if let Some(&(_, accession, name)) = SCAN_MODES
        .iter()
        .find(|t| t.0 == spectrum.instrument_settings.scan_mode)
    {
        cv(w, accession, name, "", "")?;
    }
    match spectrum.instrument_settings.polarity {
        Polarity::Positive => cv(w, "MS:1000130", "positive scan", "", "")?,
        Polarity::Negative => cv(w, "MS:1000129", "negative scan", "", "")?,
        Polarity::Unknown => {}
    }
    Ok(())
}
pub(super) fn write_file_content(w: &mut impl Write, experiment: &MSExperiment) -> Result<()> {
    let mut present = [false; ScanMode::ALL.len()];
    for spectrum in &experiment.spectra {
        present[spectrum.instrument_settings.scan_mode as usize] = true;
    }
    // Fixed source emission order; no heap allocation or string-key map.
    for mode in [
        ScanMode::MassSpectrum,
        ScanMode::Ms1Spectrum,
        ScanMode::MsnSpectrum,
        ScanMode::SelectedIonMonitoring,
        ScanMode::SelectedReactionMonitoring,
        ScanMode::ConsecutiveReactionMonitoring,
        ScanMode::Precursor,
        ScanMode::ConstantNeutralGain,
        ScanMode::ConstantNeutralLoss,
        ScanMode::ElectromagneticRadiation,
        ScanMode::Emission,
        ScanMode::Absorption,
        ScanMode::EnhancedMultiplyCharged,
        ScanMode::TimeDelayedFragmentation,
    ] {
        if present[mode as usize] {
            let &(_, accession, name) =
                SCAN_MODES.iter().find(|t| t.0 == mode).expect("known mode");
            cv(w, accession, name, "", "")?;
        }
    }
    if present[ScanMode::Unknown as usize] || experiment.spectra.is_empty() {
        cv(w, "MS:1000294", "mass spectrum", "", "")?;
    }
    Ok(())
}
pub(super) fn write_scan(w: &mut impl Write, spectrum: &MSSpectrum) -> Result<()> {
    let settings = &spectrum.instrument_settings;
    let info = &spectrum.acquisition_info;
    // MzMLHandler.cpp:5412-5440 writes a FAIMS voltage even at the -1 sentinel,
    // and every other mobility unit only for a set drift time.
    let mobility = spectrum.has_drift_time()
        || spectrum.drift_time_unit == crate::metadata::DriftTimeUnit::FaimsCompensationVoltage;
    if spectrum.rt == -1.0
        && !mobility
        && !settings.zoom_scan
        && settings.scan_windows.is_empty()
        && info.acquisitions.is_empty()
        && info.metadata.is_empty()
        && info.method_of_combination.is_empty()
    {
        return Ok(());
    }
    let count = info.acquisitions.len().max(1);
    writeln!(w, "<scanList count=\"{count}\">")?;
    let (accession, name) = acquisition_metadata::COMBINATIONS
        .iter()
        .copied()
        .find(|t| t.1 == info.method_of_combination)
        .unwrap_or(("MS:1000795", "no combination"));
    cv(w, accession, name, "", "")?;
    write_scalar_metadata(w, &info.metadata, None)?;
    for index in 0..count {
        write!(w, "<scan")?;
        if let Some(scan) = info.acquisitions.get(index) {
            if !scan.identifier.is_empty() {
                write!(w, " externalSpectrumID=\"{}\"", escape(&scan.identifier))?;
            }
            if let Some(id) = scan.metadata.get("instrument_configuration_ref") {
                write!(
                    w,
                    " instrumentConfigurationRef=\"{}\"",
                    escape(id.as_str()?)
                )?;
            }
        }
        writeln!(w, ">")?;
        if index == 0 && spectrum.rt != -1.0 {
            cv(
                w,
                "MS:1000016",
                "scan start time",
                &spectrum.rt.to_string(),
                SECOND,
            )?;
        }
        if index == 0 && mobility {
            let Some((accession, name, unit)) =
                precursor_metadata::mobility_cv(spectrum.drift_time_unit)
            else {
                return Err(invalid("spectrum mobility unit missing after preflight"));
            };
            cv(w, accession, name, &spectrum.drift_time.to_string(), unit)?;
        }
        // ParamGroupType requires all CVs before userParams. The source emits
        // zoom later, producing invalid XML when acquisition metadata is present.
        if settings.zoom_scan {
            cv(w, "MS:1000497", "zoom scan", "", "")?;
        }
        if let Some(scan) = info.acquisitions.get(index) {
            write_scalar_metadata(w, &scan.metadata, Some("instrument_configuration_ref"))?;
        }
        if index == 0 && !settings.scan_windows.is_empty() {
            writeln!(
                w,
                "<scanWindowList count=\"{}\">",
                settings.scan_windows.len()
            )?;
            for window in &settings.scan_windows {
                writeln!(w, "<scanWindow>")?;
                let accession = window
                    .metadata
                    .get("unit_accession")
                    .map(|v| v.as_str().expect("validated unit"))
                    .unwrap_or("MS:1000040");
                let (prefix, name) = unit_identity(accession)?;
                for (value, term, label) in [
                    (window.begin, "MS:1000501", "scan window lower limit"),
                    (window.end, "MS:1000500", "scan window upper limit"),
                ] {
                    write!(
                        w,
                        "<cvParam cvRef=\"MS\" accession=\"{term}\" name=\"{label}\" value=\"{value}\" unitAccession=\"{accession}\" unitCvRef=\"{prefix}\""
                    )?;
                    if let Some(name) = name {
                        write!(w, " unitName=\"{name}\"")?;
                    }
                    writeln!(w, "/>")?;
                }
                write_scalar_metadata(w, &window.metadata, Some("unit_accession"))?;
                writeln!(w, "</scanWindow>")?;
            }
            writeln!(w, "</scanWindowList>")?;
        }
        writeln!(w, "</scan>")?;
    }
    writeln!(w, "</scanList>")?;
    Ok(())
}
pub(super) fn write_chromatogram(w: &mut impl Write, c: &MSChromatogram) -> Result<()> {
    if let Some(value) = c.metadata.get("chromatogram type accession") {
        let id = value.as_str()?;
        let name = match id {
            "MS:1003019" => "pressure chromatogram",
            "MS:1003020" => "flow rate chromatogram",
            "MS:1000626" => "chromatogram type",
            _ => return Err(invalid("unsupported chromatogram type override")),
        };
        return cv(w, id, name, "", "");
    }
    let &(_, accession, name) = CHROMATOGRAM_TYPES
        .iter()
        .find(|t| t.0 == c.chromatogram_type)
        .expect("validated chromatogram type");
    cv(w, accession, name, "", "")
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Instrument field routing. The owner precharges complete input/metadata before
//! these allocation-light operations and handles dictionary descendant fallback.

use super::instrument_terms as terms;
use crate::data_structures::list::ListParse;
use crate::{Error, Result, metadata::*};

pub(super) enum Target<'a> {
    Instrument(&'a mut Instrument),
    Source(&'a mut IonSource),
    Analyzer(&'a mut MassAnalyzer),
    Detector(&'a mut IonDetector),
}
fn invalid(text: &str) -> Error {
    Error::InvalidValue(text.into())
}
fn double(text: &str) -> Result<f64> {
    let value = f64::from_list_item(text)?;
    if !value.is_finite() {
        return Err(invalid("nonfinite instrument number"));
    }
    Ok(value)
}
fn integer(text: &str) -> Result<i32> {
    i32::from_list_item(text)
}
fn put(meta: &mut MetaInfo, name: &str, value: MetaValue) -> Result<()> {
    if meta.insert(name.into(), value).is_some() {
        return Err(invalid("duplicate instrument metadata key"));
    }
    Ok(())
}
/// Return false only when dictionary-dependent fallback or ignore handling is
/// still needed. A failed mutation belongs to the caller's unpublished draft.
pub(super) fn apply(target: Target<'_>, id: &str, raw: &str, value: MetaValue) -> Result<bool> {
    match target {
        Target::Instrument(i) => {
            if let Some(v) = terms::ion_optics(id) {
                i.ion_optics = v;
                return Ok(true);
            }
            let key = match id {
                "MS:1000031" => return Ok(true),
                "MS:1000032" => {
                    i.customizations = raw.into();
                    return Ok(true);
                }
                "MS:1000529" => "instrument serial number",
                "MS:1000236" => "transmission",
                "MS:1000304" => "accelerating voltage",
                "MS:1000308" => "electric field strength",
                "MS:1000216" => {
                    put(&mut i.metadata, "field-free region", "true".into())?;
                    return Ok(true);
                }
                "MS:1000319" => {
                    put(&mut i.metadata, "space charge effect", "true".into())?;
                    return Ok(true);
                }
                _ => return Ok(false),
            };
            put(&mut i.metadata, key, value)?;
        }
        Target::Source(i) => {
            if let Some(v) = terms::inlet(id) {
                i.inlet_type = v;
                return Ok(true);
            }
            if let Some(v) = terms::ionization(id) {
                i.ionization_method = v;
                return Ok(true);
            }
            let (key, fixed) = match id {
                "MS:1000392" => ("ionization efficiency", None),
                "MS:1000486" => ("source potential", None),
                "MS:1000875" => ("declustering potential", None),
                "MS:1000876" => ("cone voltage", None),
                "MS:1000877" => ("tube lens", None),
                "MS:1000843" => ("wavelength", None),
                "MS:1000844" => ("focus diameter x", None),
                "MS:1000845" => ("focus diameter y", None),
                "MS:1000846" => ("pulse energy", None),
                "MS:1000847" => ("pulse duration", None),
                "MS:1000848" => ("attenuation", None),
                "MS:1000849" => ("impact angle", None),
                "MS:1000850" => ("laser type", Some("gas laser")),
                "MS:1000851" => ("laser type", Some("solid-state laser")),
                "MS:1000852" => ("laser type", Some("dye-laser")),
                "MS:1000853" => ("laser type", Some("free electron laser")),
                "MS:1000834" => ("matrix solution", None),
                "MS:1000835" => ("matrix solution concentration", None),
                "MS:1000836" => ("matrix application type", Some("dried dropplet")),
                "MS:1000837" => ("matrix application type", Some("printed")),
                "MS:1000838" => ("matrix application type", Some("sprayed")),
                "MS:1000839" => ("matrix application type", Some(" precoated plate")),
                _ => return Ok(false),
            };
            put(&mut i.metadata, key, fixed.map_or(value, MetaValue::from))?;
        }
        Target::Analyzer(i) => {
            if let Some(v) = terms::analyzer(id) {
                i.analyzer_type = v;
                return Ok(true);
            }
            if let Some(v) = terms::reflectron(id) {
                i.reflectron_state = v;
                return Ok(true);
            }
            match id {
                "MS:1000014" => i.accuracy = double(raw)?,
                "MS:1000022" => i.tof_total_path_length = double(raw)?,
                "MS:1000024" => i.final_ms_exponent = integer(raw)?,
                "MS:1000025" => i.magnetic_field_strength = double(raw)?,
                _ => return Ok(false),
            }
        }
        Target::Detector(i) => {
            if let Some(v) = terms::detector(id) {
                i.detector_type = v;
                return Ok(true);
            }
            if let Some(v) = terms::acquisition(id) {
                i.acquisition_mode = v;
                return Ok(true);
            }
            match id {
                "MS:1000028" => i.resolution = double(raw)?,
                "MS:1000029" => i.adc_sampling_frequency = double(raw)?,
                _ => return Ok(false),
            }
        }
    }
    Ok(true)
}

/// Scalar/source-model representability, following complete payload preflight.
/// Source omitted fields remain explicit loss errors; scientific hardware is
/// never fabricated merely to satisfy an incomplete component list.
pub(super) fn validate(i: &Instrument) -> Result<()> {
    if !i.vendor.is_empty() || !i.model.is_empty() {
        return Err(Error::Unsupported(
            "mzML cannot preserve separate instrument vendor/model fields".into(),
        ));
    }
    let any =
        !i.ion_sources.is_empty() || !i.mass_analyzers.is_empty() || !i.ion_detectors.is_empty();
    if any
        && (i.ion_sources.is_empty() || i.mass_analyzers.is_empty() || i.ion_detectors.is_empty())
    {
        return Err(Error::Unsupported(
            "mzML instrument components require a source, analyzer and detector".into(),
        ));
    }
    terms::write_ion_optics(i.ion_optics)?;
    for v in &i.ion_sources {
        if v.polarity != Polarity::Unknown {
            return Err(Error::Unsupported(
                "mzML instrument source polarity has no source mapping".into(),
            ));
        }
        terms::write_inlet(v.inlet_type)?;
        if v.metadata.contains_key("ionization accession") {
            if v.ionization_method != IonizationMethod::Unknown {
                return Err(Error::Unsupported(
                    "ionization accession would override stored ionization method".into(),
                ));
            }
        } else {
            terms::write_ionization(v.ionization_method)?;
        }
    }
    for v in &i.mass_analyzers {
        if v.resolution_method != ResolutionMethod::Unknown
            || v.resolution_type != ResolutionType::Unknown
            || v.scan_direction != ScanDirection::Unknown
            || v.scan_law != ScanLaw::Unknown
            || [v.resolution, v.scan_rate, v.scan_time, v.isolation_width]
                .iter()
                .any(|v| v.to_bits() != 0)
        {
            return Err(Error::Unsupported(
                "mzML cannot preserve omitted mass-analyzer fields".into(),
            ));
        }
        if ![
            v.accuracy,
            v.tof_total_path_length,
            v.magnetic_field_strength,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(invalid("nonfinite mass-analyzer field"));
        }
        terms::write_reflectron(v.reflectron_state)?;
        if v.metadata.contains_key("mass analyzer accession") {
            if v.analyzer_type != AnalyzerType::Unknown {
                return Err(Error::Unsupported(
                    "mass analyzer accession would override stored analyzer type".into(),
                ));
            }
        } else {
            terms::write_analyzer(v.analyzer_type)?;
        }
    }
    for v in &i.ion_detectors {
        if ![v.resolution, v.adc_sampling_frequency]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(invalid("nonfinite ion-detector field"));
        }
        terms::write_acquisition(v.acquisition_mode)?;
        terms::write_detector(v.detector_type)?;
    }
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source scan-list combination and ordered acquisition metadata transport.
use super::*;
use crate::metadata::AcquisitionInfo;

/// Whether to retain the source reader's implicit dummy scan.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AcquisitionMode {
    /// Collapse only a single empty scan with no list metadata and no meaningful
    /// combination. Informative scans and every neighbor remain in source order.
    #[default]
    Canonical,
    /// Materialize every scan, including the dummy scan emitted by mzML writers.
    Source,
}

pub(super) const COMBINATIONS: [(&str, &str); 4] = [
    ("MS:1000571", "sum of spectra"),
    ("MS:1000573", "median of spectra"),
    ("MS:1000575", "mean of spectra"),
    ("MS:1000795", "no combination"),
];

// These nine terms are the acquisition-level fallback subset of the 28 scan
// attribute/direction/law descendants in the pinned PSI-MS semantic mapping.
// Other known scan terms feed RT/instrument or unrepresented spectrum fields.
const SCAN_VALUES: [(&str, &str); 9] = [
    ("MS:1000927", "xsd:double"), // ion injection time
    ("MS:1002082", "xsd:double"), // first column elution time
    ("MS:1002083", "xsd:double"), // second column elution time
    ("MS:1002527", "xsd:string"), // instrument specific scan attribute
    ("MS:1002528", "xsd:string"), // synchronous prefilter selection
    ("MS:1002892", "xsd:string"), // ion mobility attribute
    ("MS:1003057", "xsd:int"),    // scan number: source signed i32 cast
    ("MS:1003371", "xsd:double"), // SelexION compensation voltage
    ("MS:1003394", "xsd:double"), // SelexION separation voltage
];

pub(super) fn read_cv(
    record: &mut Record,
    parent: &str,
    attrs: &BTreeMap<String, String>,
) -> Result<bool> {
    let accession = required(attrs, "accession")?;
    if parent == "scanList" {
        if let Some(&(_, name)) = COMBINATIONS.iter().find(|t| t.0 == accession) {
            if !record.seen_fields.insert("acquisition_combination") {
                return Err(invalid("duplicate scan-list combination"));
            }
            record
                .spectrum
                .as_mut()
                .ok_or_else(|| invalid("scan-list CV outside spectrum"))?
                .acquisition_info
                .method_of_combination = name.into();
            return Ok(true);
        }
    } else if parent == "scan" {
        if let Some(&(_, kind)) = SCAN_VALUES.iter().find(|t| t.0 == accession) {
            let value = scalar_user_value(attrs, kind)?;
            insert(record, parent, accession, value)?;
            return Ok(true);
        }
    }
    Ok(false)
}
pub(super) fn insert(
    record: &mut Record,
    parent: &str,
    name: &str,
    value: MetaValue,
) -> Result<()> {
    let info = &mut record
        .spectrum
        .as_mut()
        .ok_or_else(|| invalid("acquisition metadata outside spectrum"))?
        .acquisition_info;
    let metadata = if parent == "scanList" {
        &mut info.metadata
    } else {
        if !record.scan_active {
            return Err(invalid("acquisition metadata outside active scan"));
        }
        &mut info
            .acquisitions
            .last_mut()
            .ok_or_else(|| invalid("missing acquisition"))?
            .metadata
    };
    if metadata.insert(name.into(), value).is_some() {
        return Err(invalid("duplicate acquisition metadata key"));
    }
    Ok(())
}
pub(super) fn normalize(info: &mut AcquisitionInfo, mode: AcquisitionMode) {
    if mode == AcquisitionMode::Canonical
        && info.acquisitions.len() == 1
        && info.acquisitions[0].identifier.is_empty()
        && info.acquisitions[0].metadata.is_empty()
        && info.metadata.is_empty()
        && matches!(info.method_of_combination.as_str(), "" | "no combination")
    {
        *info = AcquisitionInfo::default();
    }
}
pub(super) fn validate(info: &AcquisitionInfo) -> Result<()> {
    if !info.method_of_combination.is_empty()
        && !COMBINATIONS
            .iter()
            .any(|t| t.1 == info.method_of_combination)
    {
        return Err(Error::Unsupported(
            "unknown mzML acquisition combination".into(),
        ));
    }
    validate_scalar_metadata(&info.metadata)?;
    for scan in &info.acquisitions {
        xml_string(&scan.identifier)?;
        validate_scalar_metadata(&scan.metadata)?;
    }
    Ok(())
}

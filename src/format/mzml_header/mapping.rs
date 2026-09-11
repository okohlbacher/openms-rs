// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Pinned source mapping-path membership only; not full semantic validation.
use super::*;
pub(super) fn permitted(owner: &str, id: &str, work: &mut Work) -> Result<bool> {
    if owner == "spectrum" && record_transport::SPECTRUM_CV.iter().any(|row| row.0 == id) {
        return Ok(true);
    }
    let rules: &[(&str, bool, bool)] = match owner {
        "analyzer" => &[("MS:1000443", true, true), ("MS:1000480", false, true)],
        "binaryDataArray" => &[
            ("MS:1000513", false, true),
            ("MS:1000518", false, true),
            ("MS:1000572", false, true),
            ("MS:1000513", false, true),
            ("MS:1000518", false, true),
            ("MS:1000572", false, true),
        ],
        "contact" => &[
            ("MS:1000586", true, false),
            ("MS:1000590", true, false),
            ("MS:1000585", false, true),
        ],
        "detector" => &[
            ("MS:1000026", true, true),
            ("MS:1000027", false, true),
            ("MS:1000481", false, true),
        ],
        "instrumentConfiguration" => &[
            ("MS:1000031", true, true),
            ("MS:1000496", false, true),
            ("MS:1000597", false, true),
            ("MS:1000487", false, true),
        ],
        "processingMethod" => &[("MS:1000452", false, true), ("MS:1000630", false, true)],
        "run" => &[("MS:1000857", false, true)],
        "sample" => &[
            ("MS:1000548", false, true),
            ("PATO:0001241", false, true),
            ("GO:0005575", false, true),
            ("BTO:0000000", false, true),
        ],
        "software" => &[("MS:1000531", false, true)],
        "source" => &[
            ("MS:1000008", true, true),
            ("MS:1000007", false, true),
            ("MS:1000482", false, true),
            ("MS:1000841", false, true),
            ("MS:1000842", false, true),
            ("MS:1000832", false, true),
            ("MS:1000833", false, true),
        ],
        "sourceFile" => &[
            ("MS:1000560", false, true),
            ("MS:1000561", false, true),
            ("MS:1000767", false, true),
        ],
        _ => return Ok(false),
    };
    for &(term, direct, children) in rules {
        work.charge(1, 0)?;
        if (direct && id == term) || (children && work.child(id, term)?) {
            return Ok(true);
        }
    }
    Ok(false)
}

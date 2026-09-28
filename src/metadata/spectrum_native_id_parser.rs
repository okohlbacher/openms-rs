// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Parser for extracting scan numbers from spectrum native IDs.
//!
//! Ports `src/openms/include/OpenMS/METADATA/SpectrumNativeIDParser.h` and
//! `src/openms/source/METADATA/SpectrumNativeIDParser.cpp` (core SDK
//! `bc9cc12`). See `docs/SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md`.
//!
//! Native IDs are vendor-specific spectrum identifiers that encode a scan
//! number, a file index or an experiment number. The source recognises them by
//! the PSI-MS native-ID format accession of the run:
//!
//! | CV accession | Vendor/format | Native ID pattern | Example |
//! |---|---|---|---|
//! | MS:1000768 | Thermo | `scan=NUMBER` | `scan=42` |
//! | MS:1000769 | Waters | `function=X process=Y scan=NUMBER` | `function=2 process=1 scan=100` |
//! | MS:1000770 | WIFF (AB Sciex) | `sample=X period=Y cycle=Z experiment=W` | `sample=1 period=1 cycle=42 experiment=1` |
//! | MS:1000771 | Bruker/Agilent | `scan=NUMBER` | `scan=42` |
//! | MS:1000772 | Bruker BAF | `scan=NUMBER` | `scan=42` |
//! | MS:1000773 | Bruker FID | `file=NUMBER` | `file=42` |
//! | MS:1000774 | Index-based | `index=NUMBER` | `index=42` (returns 43) |
//! | MS:1000775 | Single peak list | `file=NUMBER` | `file=42` |
//! | MS:1000776 | Thermo/Bruker TDF | `scan=NUMBER` | `scan=42` |
//! | MS:1000777 | Generic spectrum | `spectrum=NUMBER` | `spectrum=42` |
//! | MS:1001508 | Agilent MassHunter | `scanId=NUMBER` | `scanId=42` |
//! | MS:1001530 | Numeric-only | `NUMBER` | `42` |
//! | MS:1002818 | Bruker TDF | `frame=F scan=NUMBER [...]` | `frame=2 scan=529` |
//!
//! The last row is not in the header's table but is in the implementation's
//! `scan=` list, and the class test pins it.
//!
//! The source's C++ class holds only static functions, so it becomes this
//! module of free functions. Its regular expressions are Boost's; they run here
//! through [`BoostRegex`], whose
//! [`tokens`](crate::concept::boost_regex::BoostRegex::tokens) is
//! `boost::sregex_token_iterator`, so "use the last matching sub-group" means
//! exactly what it means in the source. `\d` is ASCII there, as in Boost's C
//! locale.
//!
//! # Example
//!
//! ```
//! use openms::metadata::spectrum_native_id_parser::{
//!     extract_scan_number, is_native_id, regex_from_native_id,
//! };
//!
//! // Extract a scan number using the CV accession of the native-ID format.
//! let thermo = extract_scan_number("scan=42", "MS:1000768").unwrap();
//! assert_eq!(thermo.scan, Some(42));
//! // Determine the regular expression from the native ID itself.
//! assert_eq!(regex_from_native_id("scan=123"), r"scan=(?<GROUP>\d+)");
//! // Check whether a string is a native ID at all.
//! assert!(is_native_id("scan=123"));
//! ```

use std::sync::OnceLock;

use crate::concept::boost_regex::BoostRegex;
use crate::{Error, Result};

/// Native-ID prefixes [`is_native_id`] recognises, in the source's test order.
///
/// `frame=` (Bruker TDF, MS:1002818) is recognised by the implementation
/// although the header's documentation does not list it.
pub const NATIVE_ID_PREFIXES: [&str; 10] = [
    "scan=",
    "scanId=",
    "scanID=",
    "controllerType=",
    "function=",
    "sample=",
    "index=",
    "spectrum=",
    "file=",
    "frame=",
];

/// Accessions whose native IDs carry the scan number as `scan=NUMBER`.
///
/// Thermo (MS:1000768), Waters (MS:1000769), Bruker/Agilent YEP (MS:1000771),
/// Bruker BAF (MS:1000772), Thermo/Bruker TDF scan (MS:1000776) and Bruker TDF
/// (MS:1002818, `frame=<int> scan=<int>` plus optional trailing tokens such as
/// `precursor=`, which the source's regular expression ignores).
pub const SCAN_ACCESSIONS: [&str; 6] = [
    "MS:1000768",
    "MS:1000769",
    "MS:1000771",
    "MS:1000772",
    "MS:1000776",
    "MS:1002818",
];

/// Accessions whose native IDs carry the number as `file=NUMBER`: Bruker FID
/// (MS:1000773) and single peak list (MS:1000775).
pub const FILE_ACCESSIONS: [&str; 2] = ["MS:1000773", "MS:1000775"];

/// The WIFF native-ID format accession, whose scan number is
/// `cycle * 1000 + experiment`.
pub const WIFF_ACCESSION: &str = "MS:1000770";

/// The index-based native-ID format accession, whose scan number is `index + 1`.
pub const INDEX_ACCESSION: &str = "MS:1000774";

/// The Agilent MassHunter native-ID format accession (`scanId=NUMBER`).
pub const SCAN_ID_ACCESSION: &str = "MS:1001508";

/// The generic `spectrum=NUMBER` native-ID format accession.
pub const SPECTRUM_ACCESSION: &str = "MS:1000777";

/// The numeric-only native-ID format accession (`NUMBER`).
pub const NUMERIC_ACCESSION: &str = "MS:1001530";

const SCAN_PATTERN: &str = r"scan=(?<GROUP>\d+)";
const INDEX_PATTERN: &str = r"index=(?<GROUP>\d+)";
const SCAN_ID_LOWER_PATTERN: &str = r"scanId=(?<GROUP>\d+)";
const SCAN_ID_UPPER_PATTERN: &str = r"scanID=(?<GROUP>\d+)";
const SPECTRUM_PATTERN: &str = r"spectrum=(?<GROUP>\d+)";
const FILE_PATTERN: &str = r"file=(?<GROUP>\d+)";
const NUMBER_PATTERN: &str = r"(?<GROUP>\d+)";
/// Both groups carry the name `GROUP`, which Boost allows; the token iterator
/// addresses them by index, 1 and 2.
const WIFF_PATTERN: &str = r"cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+)";

/// Whether `id` is a native ID from a vendor file: it starts with one of
/// [`NATIVE_ID_PREFIXES`].
///
/// The test is a case-sensitive prefix test, so `SCAN=123`, ` scan=123` and
/// `scan =123` are not native IDs, and a bare number is not one either.
pub fn is_native_id(id: &str) -> bool {
    NATIVE_ID_PREFIXES
        .iter()
        .any(|prefix| id.starts_with(prefix))
}

/// The regular expression that extracts the scan or index number from native
/// IDs shaped like `native_id`, with the number in the named group `GROUP`.
///
/// Chosen by the prefix of `native_id`, first match wins:
///
/// - `scan=`, `controllerType=`, `function=` and `frame=` give
///   `scan=(?<GROUP>\d+)`: Thermo, Waters and Bruker TDF IDs carry their scan
///   number in a `scan=` token even when it is not the first;
/// - `index=` gives `index=(?<GROUP>\d+)`;
/// - `scanId=` and `scanID=` give `scanId=(?<GROUP>\d+)` and
///   `scanID=(?<GROUP>\d+)`, preserving the case of the ID;
/// - `spectrum=` gives `spectrum=(?<GROUP>\d+)`;
/// - `file=` gives `file=(?<GROUP>\d+)`;
/// - anything else, including an empty string, `sample=` (WIFF) and a plain
///   number, gives `(?<GROUP>\d+)`.
///
/// The returned text is Boost syntax, and [`BoostRegex::new`] accepts it.
pub fn regex_from_native_id(native_id: &str) -> &'static str {
    if ["scan=", "controllerType=", "function=", "frame="]
        .iter()
        .any(|prefix| native_id.starts_with(prefix))
    {
        SCAN_PATTERN
    } else if native_id.starts_with("index=") {
        INDEX_PATTERN
    } else if native_id.starts_with("scanId=") {
        SCAN_ID_LOWER_PATTERN
    } else if native_id.starts_with("scanID=") {
        SCAN_ID_UPPER_PATTERN
    } else if native_id.starts_with("spectrum=") {
        SPECTRUM_PATTERN
    } else if native_id.starts_with("file=") {
        FILE_PATTERN
    } else {
        NUMBER_PATTERN
    }
}

/// Extract the scan number from `native_id` using `scan_regexp`, failing when
/// none can be extracted; the source's `extractScanNumber(native_id,
/// scan_regexp, no_error = false)`.
///
/// Every non-overlapping match of `scan_regexp` contributes its first capture
/// group (group 1, whatever its name — the source's documentation asks for a
/// group named `SCAN` but the code reads the first group), and the **last**
/// such group is converted as `StringUtils::toInt32` does: surrounding space,
/// tab, newline and carriage return are skipped, one `+` may precede the
/// number, and the whole text must be a decimal `i32`. A group that did not
/// take part in the last match converts as the empty string, which fails.
///
/// # Errors
///
/// [`Error::InvalidValue`] where the source throws `Exception::ParseError`
/// ("Could not extract scan number"): no match, or a last group that does not
/// convert. Also [`Error::InvalidValue`] when a search exceeds the regular
/// expression's work budget (see [`BoostRegex::tokens`]).
pub fn extract_scan_number_with(native_id: &str, scan_regexp: &BoostRegex) -> Result<i32> {
    match try_extract_scan_number_with(native_id, scan_regexp)? {
        Some(scan) => Ok(scan),
        None => Err(Error::InvalidValue(format!(
            "{native_id}: Could not extract scan number"
        ))),
    }
}

/// Extract the scan number from `native_id` using `scan_regexp`, or `None`;
/// the source's `extractScanNumber(native_id, scan_regexp, no_error = true)`.
///
/// Matching and conversion are as for [`extract_scan_number_with`]. The source
/// returns `-1` where this returns `None`; a caller-supplied expression that
/// captures a literal `-1` therefore reads as a failure there and as
/// `Some(-1)` here.
///
/// # Errors
///
/// Only [`Error::InvalidValue`] when a search exceeds the regular expression's
/// work budget, which the source does not have; an unextractable scan number
/// is `Ok(None)`.
pub fn try_extract_scan_number_with(
    native_id: &str,
    scan_regexp: &BoostRegex,
) -> Result<Option<i32>> {
    let matches = sub_matches(native_id, scan_regexp, &[1])?;
    Ok(matches.last().and_then(|value| to_int32(value)))
}

/// The outcome of [`extract_scan_number`]: the scan number and the warning the
/// source would have logged.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScanNumberExtraction {
    /// The scan number, or `None` where the source returns `-1`.
    pub scan: Option<i32>,
    /// The text the source writes to `OPENMS_LOG_WARN` on this path, if any.
    /// The crate has no global log stream, so the caller decides whether to
    /// print it.
    pub warning: Option<String>,
}

impl ScanNumberExtraction {
    /// The value the source returns: the scan number, or `-1` for `None`.
    pub fn source_value(&self) -> i32 {
        self.scan.unwrap_or(-1)
    }

    fn found(scan: i32) -> Self {
        Self {
            scan: Some(scan),
            warning: None,
        }
    }

    fn failed(warning: String) -> Self {
        Self {
            scan: None,
            warning: Some(warning),
        }
    }
}

/// Extract the scan number from `native_id` using the native-ID format named
/// by the PSI-MS accession `native_id_type_accession`; the source's
/// `extractScanNumber(native_id, native_id_type_accession)`.
///
/// The accession selects the expression, as the module table lists; the last
/// match wins, so a merged spectrum whose native ID repeats the format (for
/// example `... scan=100 merged ... scan=200`) yields the merged scan's number.
/// Two formats transform the number, as the source's notes state:
///
/// - WIFF (MS:1000770) returns `cycle * 1000 + experiment`, from the `cycle=`
///   and `experiment=` numbers of the last match, separated by any whitespace;
/// - index-based IDs (MS:1000774) return `index + 1`, the scan number pepXML
///   and similar consumers expect.
///
/// An unknown accession, no match, or a number that does not convert returns
/// `scan: None` with the source's warning text.
///
/// # Notes
///
/// The source computes `cycle * 1000 + experiment` and `index + 1` in `int`, so
/// a large enough value is signed overflow, which C++ leaves undefined. Here the
/// arithmetic is checked, and an overflow returns `scan: None` with a warning
/// naming the overflow instead of a wrapped number.
///
/// # Errors
///
/// [`Error::InvalidValue`] for a WIFF ID whose experiment number is 1000 or
/// more, where the source throws `Exception::InvalidValue` because the
/// `cycle * 1000 + experiment` encoding cannot represent it; and when a search
/// exceeds the regular expression's work budget.
pub fn extract_scan_number(
    native_id: &str,
    native_id_type_accession: &str,
) -> Result<ScanNumberExtraction> {
    let accession = native_id_type_accession;
    let (pattern, wiff) = if SCAN_ACCESSIONS.contains(&accession) {
        (SCAN_PATTERN, false)
    } else if accession == WIFF_ACCESSION {
        (WIFF_PATTERN, true)
    } else if FILE_ACCESSIONS.contains(&accession) {
        (FILE_PATTERN, false)
    } else if accession == INDEX_ACCESSION {
        (INDEX_PATTERN, false)
    } else if accession == SCAN_ID_ACCESSION {
        (SCAN_ID_LOWER_PATTERN, false)
    } else if accession == SPECTRUM_ACCESSION {
        (SPECTRUM_PATTERN, false)
    } else if accession == NUMERIC_ACCESSION {
        (NUMBER_PATTERN, false)
    } else {
        return Ok(ScanNumberExtraction::failed(format!(
            "native_id: {native_id} accession: {accession} Could not extract scan number - no valid native_id_type_accession was provided"
        )));
    };
    let regex = fixed_regex(pattern)?;
    let subgroups: &[i32] = if wiff { &[1, 2] } else { &[1] };
    let matches = sub_matches(native_id, regex, subgroups)?;
    if matches.len() < subgroups.len() {
        return Ok(ScanNumberExtraction::failed(format!(
            "native_id '{native_id}' is invalid. Could not extract scan number."
        )));
    }
    if wiff {
        return wiff_scan_number(native_id, accession, &matches);
    }
    let value = &matches[matches.len() - 1];
    let Some(number) = to_int32(value) else {
        return Ok(ScanNumberExtraction::failed(format!(
            "Value: '{value}' could not be converted to int in string. Native ID='{native_id}'"
        )));
    };
    if accession != INDEX_ACCESSION {
        return Ok(ScanNumberExtraction::found(number));
    }
    Ok(match number.checked_add(1) {
        Some(scan) => ScanNumberExtraction::found(scan),
        None => ScanNumberExtraction::failed(format!(
            "native_id '{native_id}' scan number {number} + 1 overflows a 32-bit integer."
        )),
    })
}

/// The WIFF branch: `cycle * 1000 + experiment` from the last match's two
/// groups, with the source's check order (experiment converted and range
/// checked first, then cycle).
fn wiff_scan_number(
    native_id: &str,
    accession: &str,
    matches: &[String],
) -> Result<ScanNumberExtraction> {
    let cycle_text = &matches[matches.len() - 2];
    let experiment_text = &matches[matches.len() - 1];
    let conversion_warning = || {
        ScanNumberExtraction::failed(format!(
            "Values: '{cycle_text}', '{experiment_text}' could not be converted to int in string. Native ID='{native_id}' accession='{accession}'"
        ))
    };
    let Some(experiment) = to_int32(experiment_text) else {
        return Ok(conversion_warning());
    };
    if experiment >= 1000 {
        return Err(Error::InvalidValue(format!(
            "The value of experiment is too large and can not be handled properly.: '{experiment_text}'"
        )));
    }
    let Some(cycle) = to_int32(cycle_text) else {
        return Ok(conversion_warning());
    };
    Ok(
        match cycle
            .checked_mul(1000)
            .and_then(|value| value.checked_add(experiment))
        {
            Some(scan) => ScanNumberExtraction::found(scan),
            None => ScanNumberExtraction::failed(format!(
                "native_id '{native_id}' encodes a scan number that overflows a 32-bit integer."
            )),
        },
    )
}

/// The texts `boost::sregex_token_iterator(native_id, regex, subgroups)`
/// yields, in order; an unmatched group is the empty string, as Boost's
/// `sub_match` converts.
fn sub_matches(native_id: &str, regex: &BoostRegex, subgroups: &[i32]) -> Result<Vec<String>> {
    let haystack = native_id.as_bytes();
    let mut out = Vec::new();
    for token in regex.tokens(haystack, subgroups)? {
        let token = token?;
        out.push(String::from_utf8_lossy(token.as_bytes(haystack)).into_owned());
    }
    Ok(out)
}

/// The compiled form of one of this module's fixed expressions, built once.
fn fixed_regex(pattern: &'static str) -> Result<&'static BoostRegex> {
    static CACHE: [OnceLock<Option<BoostRegex>>; 7] = [
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
        OnceLock::new(),
    ];
    const PATTERNS: [&str; 7] = [
        SCAN_PATTERN,
        INDEX_PATTERN,
        SCAN_ID_LOWER_PATTERN,
        SPECTRUM_PATTERN,
        FILE_PATTERN,
        NUMBER_PATTERN,
        WIFF_PATTERN,
    ];
    PATTERNS
        .iter()
        .position(|candidate| *candidate == pattern)
        .and_then(|slot| {
            CACHE[slot]
                .get_or_init(|| BoostRegex::new(pattern).ok())
                .as_ref()
        })
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "native-ID expression {pattern} could not be compiled"
            ))
        })
}

/// `StringUtils::toInt32`, `None` where the source throws `ConversionError`.
fn to_int32(text: &str) -> Option<i32> {
    crate::concept::string_conversion::to_int32(text).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_fixed_expression_compiles() {
        for pattern in [
            SCAN_PATTERN,
            INDEX_PATTERN,
            SCAN_ID_LOWER_PATTERN,
            SPECTRUM_PATTERN,
            FILE_PATTERN,
            NUMBER_PATTERN,
            WIFF_PATTERN,
        ] {
            assert!(fixed_regex(pattern).is_ok(), "{pattern}");
        }
        // `scanID=` is only ever returned as text, never compiled here, but a
        // caller compiles it.
        assert!(BoostRegex::new(SCAN_ID_UPPER_PATTERN).is_ok());
    }
}

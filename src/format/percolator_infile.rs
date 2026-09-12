// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Percolator tab-separated input (`.pin`): writing, reading and the PIN
//! feature set.
//!
//! Ports `FORMAT/PercolatorInfile.h`. A `.pin` file is one header line naming
//! the columns, then one line per peptide-spectrum match. The column contract
//! Percolator parses is three mandatory leading columns — `SpecId`, `Label`,
//! `ScanNr` — then the per-PSM feature columns, then `Peptide` and `Proteins`
//! last. [`standard_feature_set`](crate::format::percolator_infile::standard_feature_set)
//! returns the leading columns plus the standard
//! features; a caller appends its search-engine-specific columns and finally
//! `Peptide` and `Proteins` before writing.
//!
//! [`stamp_pin_features`](crate::format::percolator_infile::stamp_pin_features)
//! computes those features onto the hits without writing a file, so in-process
//! Percolator training sees the same feature vectors the file round trip would
//! have produced.
//!
//! See `docs/PERCOLATOR_INFILE_SUPPORT.md` for the API mapping, the preserved
//! source conventions and the native differences.

use crate::chemistry::{AASequence, SequenceModification};
use crate::concept::constants::user_param::{ID_MERGE_INDEX, IM, ISOTOPE_ERROR};
use crate::concept::constants::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use crate::format::csv::{self, CsvFile};
use crate::format::text::TextFile;
use crate::identification::{PeptideEvidence, PeptideHit, PeptideIdentification, TargetDecoyType};
use crate::metadata::{MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::path::Path;

/// The three header columns Percolator requires first, in order.
pub const MANDATORY_COLUMNS: [&str; 3] = ["SpecId", "Label", "ScanNr"];
/// Mass and length features written between the mandatory and charge columns.
pub const MASS_COLUMNS: [&str; 4] = ["ExpMass", "CalcMass", "mass", "peplen"];
/// Enzyme and mass-delta features written after the charge columns.
pub const ENZYME_COLUMNS: [&str; 5] = ["enzN", "enzC", "enzInt", "dm", "absdm"];
/// The two columns Percolator requires last, in order.
pub const TRAILING_COLUMNS: [&str; 2] = ["Peptide", "Proteins"];
/// The last column, whose value is itself Percolator's tab-separated list of
/// the protein accessions a PSM maps to.
const PROTEINS_COLUMN: &str = TRAILING_COLUMNS[1];
/// `ScanNr` value stamped when no scan number can be extracted from the scan
/// identifier: what `SpectrumNativeIDParser::extractScanNumber` returns with
/// its `no_error` flag set, which is how the source calls it.
const NO_SCAN_NUMBER: i32 = -1;
/// Largest number of one-hot `charge<c>` columns a feature set may declare.
pub const MAX_CHARGE_COLUMNS: usize = 1024;
/// Largest number of peptide hits one stamping or writing call may process.
pub const MAX_HITS: usize = 5_000_000;
/// Largest number of rows one `.pin` file may contribute.
pub const MAX_ROWS: usize = 5_000_000;
/// Largest number of columns one `.pin` header may declare.
pub const MAX_COLUMNS: usize = 10_000;
/// Largest cumulative owned payload one call may allocate.
pub const MAX_BYTES: usize = 256 * 1024 * 1024;

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn missing(message: impl Into<String>) -> Error {
    Error::MissingInformation(message.into())
}
fn parse(line: usize, message: impl Into<String>) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}
/// `StringUtils::toStr(double)`, which is what `DataValue::toString()` calls
/// for a floating meta value and therefore what the `.pin` row contains.
fn text(value: f64) -> String {
    crate::param::value::format_float(value, true)
}
fn finite(value: f64, label: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("percolator {label} must be finite")))
    }
}

/// The standard Percolator feature columns every `.pin` file should declare.
///
/// The list is the three mandatory header columns ([`MANDATORY_COLUMNS`]),
/// then [`MASS_COLUMNS`], then one `charge<c>` column for every charge state in
/// `min_charge..=max_charge`, then [`ENZYME_COLUMNS`]. Callers append their
/// search-engine-specific extra features and finally [`TRAILING_COLUMNS`]
/// before calling [`store`]. This is the single source of truth used by the
/// Percolator adapter and any other tool that emits `.pin` for external
/// Percolator consumption.
///
/// # Arguments
///
/// * `min_charge` — lower bound of the one-hot charge columns, inclusive.
/// * `max_charge` — upper bound, inclusive. A `max_charge` below `min_charge`
///   yields no charge column at all, as the source loop does.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when the inclusive span would exceed
/// [`MAX_CHARGE_COLUMNS`] columns. The source builds the list unbounded and a
/// caller-supplied range is external input, so the span is checked here.
pub fn standard_feature_set(min_charge: i32, max_charge: i32) -> Result<Vec<String>> {
    let span = if max_charge < min_charge {
        0
    } else {
        usize::try_from(i64::from(max_charge) - i64::from(min_charge) + 1).map_err(|_| {
            Error::InvalidRange("percolator charge range is not representable".into())
        })?
    };
    if span > MAX_CHARGE_COLUMNS {
        return Err(Error::InvalidRange(format!(
            "percolator charge range spans {span} columns, more than {MAX_CHARGE_COLUMNS}"
        )));
    }
    let mut columns = Vec::with_capacity(
        MANDATORY_COLUMNS.len() + MASS_COLUMNS.len() + span + ENZYME_COLUMNS.len(),
    );
    columns.extend(MANDATORY_COLUMNS.iter().map(|c| (*c).to_owned()));
    columns.extend(MASS_COLUMNS.iter().map(|c| (*c).to_owned()));
    for charge in 0..span {
        let charge = i64::from(min_charge) + charge as i64;
        columns.push(format!("charge{charge}"));
    }
    columns.extend(ENZYME_COLUMNS.iter().map(|c| (*c).to_owned()));
    Ok(columns)
}

/// The scan identifier used for the `SpecId` column.
///
/// Prefers the spectrum reference, which is what MS-GF+ writes. When that is
/// empty the `spectrum_id` meta value is used as `scan=<value>`, which is what
/// X!Tandem writes — those identifiers are one-based, unlike the zero-based
/// `index`, which the source warns makes them risky for merging. When neither
/// is present the fallback is `index=<index>`; the source logs "no known
/// spectrum identifiers, using index \[1,n\] - use at own risk" there.
///
/// All space, tab, carriage-return and newline characters are removed from the
/// result, as the source `StringUtils::removeWhitespaces` does.
pub fn scan_identifier(identification: &PeptideIdentification, index: usize) -> String {
    let mut value = identification.spectrum_reference();
    if value.is_empty() {
        value = match identification.metadata.get("spectrum_id") {
            Some(id) if !id.to_string().is_empty() => format!("scan={id}"),
            _ => format!("index={index}"),
        };
    }
    value.retain(|c| !matches!(c, ' ' | '\t' | '\n' | '\r'));
    value
}

/// Whether the bond between the residues `n` and `c` is a cleavage site of
/// `enzyme`, as Percolator's own `Enzyme.h` defines it.
///
/// `n` is the residue before the bond and `c` the residue after it. A `-`
/// marker on either side denotes a protein terminus and always counts as
/// enzymatic. An enzyme name the source does not know returns `true`, which
/// makes every bond enzymatic; that is the source's `else` branch and is
/// preserved because the adapters pass `no_enzyme` deliberately.
///
/// Accepted names are `trypsin`, `trypsinp`, `chymotrypsin`, `thermolysin`,
/// `proteinasek`, `pepsin`, `elastase`, `lys-n`, `lys-c`, `arg-c`, `asp-n` and
/// `glu-c`.
pub fn is_enzymatic(n: char, c: char, enzyme: &str) -> bool {
    let terminus = n == '-' || c == '-';
    match enzyme {
        "trypsin" => ((n == 'K' || n == 'R') && c != 'P') || terminus,
        "trypsinp" => (n == 'K' || n == 'R') || terminus,
        "chymotrypsin" => ((n == 'F' || n == 'W' || n == 'Y' || n == 'L') && c != 'P') || terminus,
        "thermolysin" => {
            ((c == 'A'
                || c == 'F'
                || c == 'I'
                || c == 'L'
                || c == 'M'
                || c == 'V'
                || (n == 'R' && c == 'G'))
                && n != 'D'
                && n != 'E')
                || terminus
        }
        "proteinasek" => {
            (n == 'A'
                || n == 'E'
                || n == 'F'
                || n == 'I'
                || n == 'L'
                || n == 'T'
                || n == 'V'
                || n == 'W'
                || n == 'Y')
                || terminus
        }
        "pepsin" => {
            ((c == 'F'
                || c == 'L'
                || c == 'W'
                || c == 'Y'
                || n == 'F'
                || n == 'L'
                || n == 'W'
                || n == 'Y')
                && n != 'R')
                || terminus
        }
        "elastase" => ((n == 'L' || n == 'V' || n == 'A' || n == 'G') && c != 'P') || terminus,
        "lys-n" => (c == 'K') || terminus,
        "lys-c" => (n == 'K' && c != 'P') || terminus,
        "arg-c" => (n == 'R' && c != 'P') || terminus,
        "asp-n" => (c == 'D') || terminus,
        "glu-c" => (n == 'E' && c != 'P') || terminus,
        _ => true,
    }
}

/// Number of internal bonds of `peptide` that [`is_enzymatic`] accepts, the
/// `enzInt` feature.
///
/// Bonds are counted between consecutive characters, so an empty or
/// single-residue peptide has none. The source indexes the string by byte;
/// this iterates characters, so a non-ASCII input is counted rather than split
/// mid-codepoint.
pub fn count_enzymatic(peptide: &str, enzyme: &str) -> usize {
    let mut count = 0;
    let mut previous = None;
    for c in peptide.chars() {
        if let Some(n) = previous {
            if is_enzymatic(n, c, enzyme) {
                count += 1;
            }
        }
        previous = Some(c);
    }
    count
}

/// The pattern `SpectrumNativeIDParser::getRegExFromNativeID` selects for a
/// native identifier, reduced to the literal prefix its regex looks for.
///
/// The source builds a `boost::regex`; this port has no regular-expression
/// dependency, and every pattern the source can return has the shape
/// `<prefix>=(\d+)` or a bare `(\d+)`, so the prefix alone reproduces it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScanNumberPattern(Option<&'static str>);

impl ScanNumberPattern {
    /// The pattern for `identifier`, chosen by its leading token exactly as
    /// the source does. Thermo, Waters and Bruker TDF identifiers all resolve
    /// to the `scan=` token, which remains the most meaningful scan-number
    /// proxy for their trailing tokens.
    pub fn for_native_id(identifier: &str) -> Self {
        for prefix in ["scan=", "controllerType=", "function=", "frame="] {
            if identifier.starts_with(prefix) {
                return Self(Some("scan="));
            }
        }
        for prefix in ["index=", "scanId=", "scanID=", "spectrum=", "file="] {
            if identifier.starts_with(prefix) {
                return Self(Some(prefix));
            }
        }
        Self(None)
    }

    /// The scan number of `identifier`, or `None` when none can be extracted.
    ///
    /// The last match is used, as in the source, and a digit run that does not
    /// fit `i32` yields `None`; the source catches its conversion error and
    /// falls through to the same outcome. The source's `no_error` argument
    /// selects between returning `-1` and throwing `Exception::ParseError`;
    /// this returns `None` and leaves the choice to the caller.
    pub fn extract(self, identifier: &str) -> Option<i32> {
        let mut last = None;
        match self.0 {
            Some(prefix) => {
                let mut rest = identifier;
                while let Some(at) = rest.find(prefix) {
                    let tail = rest.get(at + prefix.len()..).unwrap_or("");
                    let digits: String = tail.chars().take_while(char::is_ascii_digit).collect();
                    if !digits.is_empty() {
                        last = Some(digits);
                    }
                    rest = rest.get(at + prefix.len()..).unwrap_or("");
                    if rest.is_empty() {
                        break;
                    }
                }
            }
            None => {
                let mut digits = String::new();
                for c in identifier.chars() {
                    if c.is_ascii_digit() {
                        digits.push(c);
                    } else if !digits.is_empty() {
                        last = Some(std::mem::take(&mut digits));
                    }
                }
                if !digits.is_empty() {
                    last = Some(digits);
                }
            }
        }
        last.and_then(|value| value.parse::<i32>().ok())
    }
}

/// Percolator-formatted sequence in bracket notation with relative mass deltas,
/// the `Peptide` column's middle part.
///
/// Reproduces `AASequence::toBracketString(false, true)`: accurate masses,
/// written as signed deltas, with `n[...]` and `c[...]` for terminal
/// annotations. Values are rendered with the source's full-precision
/// `StringUtils::toStr`, so `Oxidation` becomes `M[+15.9949]`.
///
/// # Errors
///
/// Returns [`Error::Unsupported`] when an annotation has no known delta mass,
/// and for a modified `X` residue, whose absolute internal mass the source
/// substitutes but which this port cannot resolve without a known residue mass.
pub fn bracket_sequence(sequence: &AASequence) -> Result<String> {
    fn delta(modification: &SequenceModification) -> Result<f64> {
        modification.diff_mono_mass()
    }
    fn signed(value: f64) -> String {
        let rendered = text(value);
        if value > 0.0 {
            format!("+{rendered}")
        } else {
            rendered
        }
    }
    let mut result = String::new();
    if sequence.is_empty() {
        return Ok(result);
    }
    if let Some(modification) = sequence.n_terminal_modification() {
        result.push_str(&format!("n[{}]", signed(delta(modification)?)));
    }
    for (index, residue) in sequence.as_str().chars().enumerate() {
        result.push(residue);
        let Some(modification) = sequence.residue_modification(index)? else {
            continue;
        };
        if residue == 'X' {
            // The source cannot express a delta for X and writes the absolute
            // internal residue mass without a sign.
            let internal = sequence
                .subsequence(index..index + 1)?
                .mono_mass_for(crate::chemistry::PeptideFragmentType::Internal, 0)?;
            result.push_str(&format!("[{}]", text(internal)));
            continue;
        }
        result.push_str(&format!("[{}]", signed(delta(modification)?)));
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        result.push_str(&format!("c[{}]", signed(delta(modification)?)));
    }
    Ok(result)
}

/// Enzyme and charge-range settings shared by feature stamping and writing.
#[derive(Clone, Debug)]
pub struct PinOptions {
    /// Enzyme name passed to [`is_enzymatic`]; an unknown name makes every
    /// bond enzymatic, as in the source.
    pub enzyme: String,
    /// Lower bound of the one-hot `charge<c>` features, inclusive.
    pub min_charge: i32,
    /// Upper bound of the one-hot `charge<c>` features, inclusive.
    pub max_charge: i32,
    /// Ceiling on the total peptide hits one call may process.
    pub max_hits: usize,
    /// Ceiling on the cumulative owned payload one call may allocate.
    pub max_bytes: usize,
}
impl Default for PinOptions {
    fn default() -> Self {
        Self {
            enzyme: "trypsin".into(),
            min_charge: 2,
            max_charge: 5,
            max_hits: MAX_HITS,
            max_bytes: MAX_BYTES,
        }
    }
}
impl PinOptions {
    fn validate(&self) -> Result<()> {
        if self.max_hits == 0 || self.max_bytes == 0 {
            return Err(invalid("percolator limits must be positive"));
        }
        // Rejects an unrepresentable span before any hit is touched.
        standard_feature_set(self.min_charge, self.max_charge)?;
        Ok(())
    }
}

/// What [`stamp_pin_features`] did, per hit.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StampReport {
    /// `(peptide identification index, hit index)` pairs left untouched
    /// because the hit had no protein evidence or no target/decoy status.
    pub skipped: BTreeSet<(usize, usize)>,
    /// Meta-value keys that did not exist on a hit before stamping, so a
    /// caller that stamps temporarily can remove only what it added.
    pub added_meta_values: BTreeMap<(usize, usize), BTreeSet<String>>,
    /// Diagnostics the source writes to its warning log.
    pub warnings: Vec<String>,
}

struct Budget {
    bytes: usize,
}
impl Budget {
    fn spend(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("percolator payload limit exceeded"))?;
        Ok(())
    }
}

fn count_hits(identifications: &[PeptideIdentification], options: &PinOptions) -> Result<usize> {
    let mut total = 0usize;
    for identification in identifications {
        total = total
            .checked_add(identification.hits.len())
            .filter(|&n| n <= options.max_hits)
            .ok_or_else(|| invalid("percolator peptide hit limit exceeded"))?;
    }
    Ok(total)
}

/// Compute and stamp the PIN meta values on every peptide hit.
///
/// Runs the per-hit computation [`prepare_pin`] applies when writing a `.pin`
/// file, but mutates the identifications in place instead. After a successful
/// call each kept hit carries `SpecId`, `ScanNr`, `Label`, `CalcMass`,
/// `ExpMass`, `deltamass`, `retentiontime`, `mass`, `score`, `peplen`,
/// `charge<min>`..`charge<max>`, `enzN`, `enzC`, `enzInt`, `dm`, `absdm`,
/// `Peptide` and `Proteins`. This is what in-process Percolator training needs
/// so it sees the feature vectors the subprocess path would have seen.
///
/// Hits with no protein evidence, and hits whose target/decoy status is
/// unknown, are left untouched and reported in [`StampReport::skipped`]; the
/// source logs a warning naming incomplete peptide indexing for both.
///
/// An existing `CalcMass` is reused rather than recomputed, as in the source.
/// An `IsotopeError` meta value — the legacy MS-GF+ adapter's spelling before
/// OpenMS 2.6 — or the current `isotope_error` shifts `ExpMass` down by
/// `isotope_error * C13C12_MASSDIFF_U / charge`.
///
/// The scan-number pattern is derived from the **first** identification's scan
/// identifier and then applied to all of them, which is what the source does;
/// an identification whose identifier has a different shape therefore gets the
/// `ScanNr` `-1`, the sentinel `extractScanNumber` returns when its `no_error`
/// flag is set — and the source sets it. See the support document.
///
/// The `Proteins` value joins the hit's accessions with tabs, because that is
/// Percolator's trailing protein list; [`prepare_pin`] writes it into the last
/// column unescaped, as the source does.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the hit count or payload ceiling in
/// `options` is reached, when a coordinate is not finite, when a charge is
/// zero and an isotope-error correction would divide by it, or when a kept
/// hit's sequence is empty, whose first and last residue the source reads past
/// the end of; [`Error::MissingInformation`] when an identification carries no
/// m/z or no retention time — the source's defaults for both are NaN, which it
/// would write into the `ExpMass`, `mass`, `dm` and `retentiontime` columns as
/// `nan`; [`Error::Unsupported`] when a sequence has no computable mass or an
/// annotation no delta mass. On any error the identifications are unchanged,
/// because the computation runs on a temporary that is committed last.
pub fn stamp_pin_features(
    identifications: &mut Vec<PeptideIdentification>,
    options: &PinOptions,
) -> Result<StampReport> {
    options.validate()?;
    count_hits(identifications, options)?;
    let (stamped, report) = stamped_copy(identifications, options)?;
    *identifications = stamped;
    Ok(report)
}

fn stamped_copy(
    identifications: &[PeptideIdentification],
    options: &PinOptions,
) -> Result<(Vec<PeptideIdentification>, StampReport)> {
    let mut report = StampReport::default();
    let mut budget = Budget {
        bytes: options.max_bytes,
    };
    let mut result = identifications.to_vec();
    if result.is_empty() {
        return Ok((result, report));
    }
    let pattern = ScanNumberPattern::for_native_id(&scan_identifier(&result[0], 0));
    for (pid_index, identification) in result.iter_mut().enumerate() {
        // The source numbers identifications from one when it builds the scan
        // identifier, but from zero when it keys the skipped/added maps.
        let one_based = pid_index
            .checked_add(1)
            .ok_or_else(|| invalid("percolator identification index overflows"))?;
        let identifier = scan_identifier(identification, one_based);
        let mut file_identifier = match identification.metadata.get("file_origin") {
            Some(value) => value.to_string(),
            None => String::new(),
        };
        if let Some(value) = identification.metadata.get(ID_MERGE_INDEX) {
            file_identifier.push_str(&value.to_string());
        }
        let scan_number = pattern.extract(&identifier);
        let experimental_mass = identification.mz;
        let retention_time = identification.rt;

        for (hit_index, hit) in identification.hits.iter_mut().enumerate() {
            let location = (pid_index, hit_index);
            if hit.evidences.is_empty() {
                report.warnings.push(
                    "PSM (PeptideHit) without protein reference found. This may indicate \
                     incomplete mapping during PeptideIndexing (e.g., wrong enzyme settings). \
                     Will skip this PSM."
                        .into(),
                );
                report.skipped.insert(location);
                continue;
            }
            if hit.target_decoy_type()? == TargetDecoyType::Unknown {
                report.warnings.push(
                    "PSM without target/decoy information found. This may indicate incomplete \
                     mapping during PeptideIndexing (e.g., wrong decoy prefix settings). Will \
                     skip this PSM."
                        .into(),
                );
                report.skipped.insert(location);
                continue;
            }

            let mut added = BTreeSet::new();
            let stamp = |hit: &mut PeptideHit,
                         budget: &mut Budget,
                         added: &mut BTreeSet<String>,
                         key: &str,
                         value: MetaValue|
             -> Result<()> {
                if !hit.metadata.contains_key(key) {
                    budget.spend(key.len().saturating_add(64))?;
                    added.insert(key.to_owned());
                }
                budget.spend(key.len().saturating_add(64))?;
                hit.metadata.insert(key.to_owned(), value);
                Ok(())
            };

            let charge = hit.charge;
            let unmodified = hit.sequence.as_str().to_owned();
            // The source calls extractScanNumber with no_error = true, which
            // returns -1 when the pattern matches nothing, and stamps that -1.
            // An identification whose identifier has a different shape from
            // the first one's therefore gets ScanNr -1 rather than aborting
            // the whole call.
            let scan = scan_number.unwrap_or(NO_SCAN_NUMBER);

            stamp(
                hit,
                &mut budget,
                &mut added,
                "SpecId",
                format!("{file_identifier}{identifier}").into(),
            )?;
            stamp(hit, &mut budget, &mut added, "ScanNr", scan.into())?;
            let label: i32 = if hit.is_decoy()? { -1 } else { 1 };
            stamp(hit, &mut budget, &mut added, "Label", label.into())?;

            let calculated = match hit.metadata.get("CalcMass") {
                Some(value) => finite(value.as_f64()?, "CalcMass")?,
                None => {
                    let value = finite(hit.sequence.mz(charge)?, "CalcMass")?;
                    stamp(
                        hit,
                        &mut budget,
                        &mut added,
                        "CalcMass",
                        MetaValue::try_from(value)?,
                    )?;
                    value
                }
            };

            let mut row_mass = finite(
                experimental_mass.ok_or_else(|| missing("percolator ExpMass requires an m/z"))?,
                "ExpMass",
            )?;
            let isotope_error = match (
                hit.metadata.get("IsotopeError"),
                hit.metadata.get(ISOTOPE_ERROR),
            ) {
                (Some(value), _) | (None, Some(value)) => Some(numeric_meta(value)?),
                (None, None) => None,
            };
            if let Some(error) = isotope_error {
                if charge == 0 {
                    return Err(invalid(
                        "percolator isotope-error correction requires a non-zero charge",
                    ));
                }
                row_mass = finite(
                    row_mass - (error * C13C12_MASSDIFF_U) / f64::from(charge),
                    "ExpMass",
                )?;
            }
            stamp(
                hit,
                &mut budget,
                &mut added,
                "ExpMass",
                MetaValue::try_from(row_mass)?,
            )?;

            let delta_mass = finite(row_mass - calculated, "deltamass")?;
            for (key, value) in [
                ("deltamass", delta_mass),
                (
                    "retentiontime",
                    finite(
                        retention_time
                            .ok_or_else(|| missing("percolator retentiontime requires an RT"))?,
                        "retentiontime",
                    )?,
                ),
                ("mass", row_mass),
                ("score", finite(hit.score, "score")?),
            ] {
                stamp(
                    hit,
                    &mut budget,
                    &mut added,
                    key,
                    MetaValue::try_from(value)?,
                )?;
            }
            let length = i64::try_from(unmodified.chars().count())
                .map_err(|_| invalid("percolator peplen exceeds i64"))?;
            stamp(hit, &mut budget, &mut added, "peplen", length.into())?;

            let mut one_hot = options.min_charge;
            while one_hot <= options.max_charge {
                // The source stores a bool, which has no DataValue alternative
                // and is therefore promoted to an int: the column reads 1 or 0.
                let value: i32 = i32::from(charge == one_hot);
                stamp(
                    hit,
                    &mut budget,
                    &mut added,
                    &format!("charge{one_hot}"),
                    value.into(),
                )?;
                one_hot = match one_hot.checked_add(1) {
                    Some(next) => next,
                    None => break,
                };
            }

            // The source takes the flanking residues from the first evidence only.
            let first = hit
                .evidences
                .first()
                .ok_or_else(|| missing("percolator flanks require protein evidence"))?;
            let mut before = first.aa_before.code();
            let mut after = first.aa_after.code();
            let n_terminal = unmodified.chars().next();
            let c_terminal = unmodified.chars().next_back();
            let enz_n = match n_terminal {
                Some(residue) => is_enzymatic(before, residue, &options.enzyme),
                // The source indexes prefix(1)[0] of the unmodified sequence,
                // which reads past the end of an empty peptide.
                None => return Err(invalid("percolator enzN requires a non-empty sequence")),
            };
            let enz_c = match c_terminal {
                Some(residue) => is_enzymatic(residue, after, &options.enzyme),
                None => return Err(invalid("percolator enzC requires a non-empty sequence")),
            };
            stamp(
                hit,
                &mut budget,
                &mut added,
                "enzN",
                i32::from(enz_n).into(),
            )?;
            stamp(
                hit,
                &mut budget,
                &mut added,
                "enzC",
                i32::from(enz_c).into(),
            )?;
            let internal = i64::try_from(count_enzymatic(&unmodified, &options.enzyme))
                .map_err(|_| invalid("percolator enzInt exceeds i64"))?;
            stamp(hit, &mut budget, &mut added, "enzInt", internal.into())?;
            stamp(
                hit,
                &mut budget,
                &mut added,
                "dm",
                MetaValue::try_from(delta_mass)?,
            )?;
            stamp(
                hit,
                &mut budget,
                &mut added,
                "absdm",
                MetaValue::try_from(finite(delta_mass.abs(), "absdm")?)?,
            )?;

            // The source maps its N/C-terminus markers onto Percolator's '-'.
            if before == '[' {
                before = '-';
            }
            if after == ']' {
                after = '-';
            }
            let bracket = bracket_sequence(&hit.sequence)?;
            budget.spend(bracket.len().saturating_add(4))?;
            stamp(
                hit,
                &mut budget,
                &mut added,
                "Peptide",
                format!("{before}.{bracket}.{after}").into(),
            )?;

            let mut proteins = String::new();
            for (index, evidence) in hit.evidences.iter().enumerate() {
                if index != 0 {
                    proteins.push('\t');
                }
                budget.spend(evidence.protein_accession.len().saturating_add(1))?;
                proteins.push_str(&evidence.protein_accession);
            }
            stamp(hit, &mut budget, &mut added, "Proteins", proteins.into())?;

            if !added.is_empty() {
                report.added_meta_values.insert(location, added);
            }
        }
    }
    Ok((result, report))
}

fn numeric_meta(value: &MetaValue) -> Result<f64> {
    match value.data() {
        MetaValueData::Integer(number) => Ok(*number as f64),
        MetaValueData::Float(number) => Ok(*number),
        // The source converts the meta value to text and then to a float, so a
        // string-typed isotope error is accepted.
        MetaValueData::String(number) => number
            .trim()
            .parse::<f64>()
            .map_err(|_| invalid("percolator isotope error is not numeric")),
        _ => Err(invalid("percolator isotope error is not numeric")),
    }
}

/// The `.pin` lines a feature set and a set of identifications produce.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PreparedPin {
    /// The header line followed by one line per written hit.
    pub lines: Vec<String>,
    /// Hits dropped because a declared feature had no meta value.
    pub hits_missing_features: usize,
    /// Names of the features that were missing on at least one dropped hit.
    pub missing_features: BTreeSet<String>,
    /// Diagnostics the source writes to its warning log, including those from
    /// feature stamping.
    pub warnings: Vec<String>,
}

/// Whether `position` is the `Proteins` column at the end of `feature_set`,
/// the one column whose value is itself a tab-separated list.
fn is_trailing_protein_column(feature_set: &[String], position: usize) -> bool {
    position.saturating_add(1) == feature_set.len()
        && feature_set.get(position).map(String::as_str) == Some(PROTEINS_COLUMN)
}

/// Build the `.pin` text: the tab-joined header, then the stamped features of
/// each kept hit in the declared column order.
///
/// The identifications are not modified; the source copies them too. A hit
/// whose stamped meta values do not cover every declared feature is dropped
/// and counted in [`PreparedPin::hits_missing_features`] rather than written
/// with blank fields, which is what the source does.
///
/// A float feature is rendered by the source's full-precision
/// `StringUtils::toStr`; an integer feature by its decimal digits; a `bool`
/// feature by `1` or `0`, because the source's `DataValue` has no boolean
/// alternative and promotes it to an int.
///
/// The final `Proteins` column is Percolator's trailing protein list, which is
/// tab-separated *inside* that one column: a hit with several protein
/// evidences therefore renders a row with one field per accession beyond the
/// declared column count. That is what
/// [`stamp_pin_features`] writes and what the source writes
/// (`PercolatorInfile.cpp:549`), so a tab is accepted there and nowhere else.
///
/// # Errors
///
/// As [`stamp_pin_features`], plus [`Error::InvalidValue`] when a feature name
/// or a rendered value contains a CR or LF, which would end the row early, and
/// when a rendered value contains a tab in any column but a trailing
/// `Proteins`, which would silently shift every column after it. The source
/// writes all of those through unescaped.
pub fn prepare_pin(
    identifications: &[PeptideIdentification],
    feature_set: &[String],
    options: &PinOptions,
) -> Result<PreparedPin> {
    options.validate()?;
    if feature_set.len() > MAX_COLUMNS {
        return Err(invalid("percolator column limit exceeded"));
    }
    for name in feature_set {
        if name.contains(['\t', '\n', '\r']) {
            return Err(invalid("percolator feature name contains a separator"));
        }
    }
    let mut result = PreparedPin {
        lines: vec![feature_set.join("\t")],
        ..Default::default()
    };
    if identifications.is_empty() {
        result
            .warnings
            .push("No identifications provided. Creating empty percolator input.".into());
        return Ok(result);
    }
    count_hits(identifications, options)?;
    let (stamped, report) = stamped_copy(identifications, options)?;
    result.warnings = report.warnings;
    for (pid_index, identification) in stamped.iter().enumerate() {
        for (hit_index, hit) in identification.hits.iter().enumerate() {
            if report.skipped.contains(&(pid_index, hit_index)) {
                continue;
            }
            let mut fields = Vec::with_capacity(feature_set.len());
            for (position, name) in feature_set.iter().enumerate() {
                let Some(value) = hit.metadata.get(name) else {
                    continue;
                };
                let rendered = match value.data() {
                    MetaValueData::Float(number) => text(*number),
                    other => other_meta_text(other),
                };
                // A CR or LF ends the row early whatever the column is.
                if rendered.contains(['\n', '\r']) {
                    return Err(invalid("percolator feature value contains a line break"));
                }
                // Percolator's trailing protein list is tab-separated inside
                // the last column, which is why stamp_pin_features joins the
                // accessions with tabs. Everywhere else a tab would shift the
                // columns that follow it.
                if rendered.contains('\t') && !is_trailing_protein_column(feature_set, position) {
                    return Err(invalid("percolator feature value contains a separator"));
                }
                fields.push(rendered);
            }
            if fields.len() == feature_set.len() {
                result.lines.push(fields.join("\t"));
            } else {
                result.hits_missing_features += 1;
                for name in feature_set {
                    if !hit.metadata.contains_key(name) {
                        result.missing_features.insert(name.clone());
                    }
                }
            }
        }
    }
    if result.hits_missing_features != 0 {
        result.warnings.push(format!(
            "There were peptide hits with missing features/meta values. Skipped peptide hits: {}",
            result.hits_missing_features
        ));
    }
    Ok(result)
}

fn other_meta_text(data: &MetaValueData) -> String {
    match data {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(value) => value.clone(),
        MetaValueData::Integer(value) => value.to_string(),
        MetaValueData::Float(value) => text(*value),
        MetaValueData::StringList(values) => format!("[{}]", values.join(", ")),
        MetaValueData::IntegerList(values) => format!(
            "[{}]",
            values
                .iter()
                .map(i64::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
        MetaValueData::FloatList(values) => format!(
            "[{}]",
            values
                .iter()
                .map(|value| text(*value))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// Write the `.pin` text produced by [`prepare_pin`] to `writer`.
///
/// # Errors
///
/// As [`prepare_pin`], plus [`Error::Io`] from the writer.
pub fn write(
    writer: impl Write,
    identifications: &[PeptideIdentification],
    feature_set: &[String],
    options: &PinOptions,
) -> Result<PreparedPin> {
    let prepared = prepare_pin(identifications, feature_set, options)?;
    let mut file = TextFile::new();
    for line in &prepared.lines {
        file.add_line(line)?;
    }
    file.write(writer)?;
    Ok(prepared)
}

/// Write a `.pin` file, publishing it only once complete.
///
/// # Errors
///
/// As [`write()`], plus [`Error::Io`] when the output cannot be created; the
/// source's `TextFile::store` throws `Exception::UnableToCreateFile` there.
pub fn store(
    path: impl AsRef<Path>,
    identifications: &[PeptideIdentification],
    feature_set: &[String],
    options: &PinOptions,
) -> Result<PreparedPin> {
    let prepared = prepare_pin(identifications, feature_set, options)?;
    let mut file = TextFile::new();
    for line in &prepared.lines {
        file.add_line(line)?;
    }
    file.store(path)?;
    Ok(prepared)
}

/// Reading settings, one per `load` argument of the source.
#[derive(Clone, Debug)]
pub struct ReadOptions {
    /// Whether a higher primary score ranks better.
    pub higher_score_better: bool,
    /// Column holding the primary score, and the score type recorded on each
    /// identification.
    pub score_name: String,
    /// Additional score columns stored verbatim on each hit as string meta
    /// values. A name that the header does not declare produces a warning.
    pub extra_scores: Vec<String>,
    /// Accession prefix that marks a decoy protein. When set, the target/decoy
    /// status is recomputed from the accessions and the file's `Label` column
    /// is overridden; when empty the `Label` column is trusted.
    pub decoy_prefix: String,
    /// Upper bound on the Sage `spectrum_q` of a kept row.
    pub spectrum_q_threshold: f64,
    /// Read the sibling Sage `.tsv` and `matched_fragments.sage.tsv` files.
    /// Only [`load`] can honour this, because the siblings are derived from the
    /// `.pin` path.
    pub sage_annotation: bool,
    /// Ceiling on the rows one file may contribute.
    pub max_rows: usize,
    /// Ceiling on the columns the header may declare.
    pub max_columns: usize,
    /// Ceiling on parsed payload and, separately, each staged input file.
    /// Applied while reading bytes, before the CSV table is materialized.
    pub max_bytes: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            higher_score_better: true,
            score_name: String::new(),
            extra_scores: Vec::new(),
            decoy_prefix: String::new(),
            spectrum_q_threshold: 0.01,
            sage_annotation: false,
            max_rows: MAX_ROWS,
            max_columns: MAX_COLUMNS,
            max_bytes: MAX_BYTES,
        }
    }
}

/// What reading a `.pin` file produced.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PinDocument {
    /// One identification per run of rows that share a `SpecId`.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// Raw file names in order of first appearance, from the `FileName` column.
    pub filenames: Vec<String>,
    /// Diagnostics the source writes to its warning log.
    pub warnings: Vec<String>,
}

struct Header {
    names: Vec<String>,
    index: BTreeMap<String, usize>,
}
impl Header {
    fn new(names: Vec<String>) -> Result<Self> {
        let mut index = BTreeMap::new();
        for (position, name) in names.iter().enumerate() {
            // The source builds an unordered_map, so a duplicate column name
            // silently resolves to the last occurrence. Reject it instead.
            if index.insert(name.clone(), position).is_some() {
                return Err(parse(1, format!("duplicate .pin column {name:?}")));
            }
        }
        Ok(Self { names, index })
    }
    fn position(&self, name: &str) -> Result<usize> {
        self.index
            .get(name)
            .copied()
            .ok_or_else(|| missing(format!(".pin file has no {name} column")))
    }
    fn optional(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }
}

fn field<'a>(row: &'a [String], position: usize, line: usize, name: &str) -> Result<&'a str> {
    row.get(position)
        .map(String::as_str)
        .ok_or_else(|| parse(line, format!("row is missing the {name} column")))
}
fn row_f64(row: &[String], position: usize, line: usize, name: &str) -> Result<f64> {
    let raw = field(row, position, line, name)?;
    raw.trim()
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| parse(line, format!("invalid finite {name}: {raw:?}")))
}
fn row_i32(row: &[String], position: usize, line: usize, name: &str) -> Result<i32> {
    let raw = field(row, position, line, name)?;
    raw.trim()
        .parse::<i32>()
        .map_err(|_| parse(line, format!("invalid integer {name}: {raw:?}")))
}

/// Charge columns a `.pin` header declares, in ascending charge order.
///
/// Percolator does not standardise the spelling: OpenMS and most engines write
/// `charge<c>`, Sage writes `z=<c>` plus a `z=other` column for matches whose
/// charge falls outside the searched range.
fn charge_columns(header: &Header) -> Vec<(usize, i32)> {
    let mut columns = Vec::new();
    for (position, name) in header.names.iter().enumerate() {
        for prefix in ["charge", "z="] {
            let Some(digits) = name.strip_prefix(prefix) else {
                continue;
            };
            if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
                continue;
            }
            if let Ok(charge) = digits.parse::<i32>() {
                columns.push((position, charge));
            }
        }
    }
    // The source iterates an unordered_map and breaks at the first column set
    // to 1, so a row with two one-hot columns set picks an unspecified charge.
    // Ascending charge order makes the same row deterministic here.
    columns.sort_by_key(|(_, charge)| *charge);
    columns
}

/// Read a `.pin` file from `reader`.
///
/// Each run of consecutive rows sharing a `SpecId` becomes one
/// [`PeptideIdentification`]; its `RT` is the `retentiontime` column times 60,
/// because search engines such as Sage write minutes, and its spectrum
/// reference is the `ScanNr` column, which may be an integer or a full vendor
/// identifier. The original `SpecId` is kept as the `PinSpecId` meta value.
///
/// Every hit gets a `target_decoy` meta value, the extra score columns as
/// string meta values, and a `DeltaMass` meta value holding
/// `ExpMass - CalcMass`. A `ln(-poisson)` value of `inf` is replaced by `3.5`,
/// a workaround the source carries for Sage.
///
/// The `Proteins` column is split on `;`, which is how Sage writes a protein
/// list, and a row must hold exactly as many fields as the header — both
/// straight from the source. Neither accepts the tab-separated trailing
/// protein list that [`prepare_pin`] and the source's own `store` write for a
/// hit with several protein evidences: such a row is
/// [`Error::Parse`] here and `Exception::ParseError` there. The asymmetry is
/// the source's and is recorded in the support document.
///
/// # Errors
///
/// Returns [`Error::Parse`] when a row does not declare the same number of
/// columns as the header — the source throws `Exception::ParseError` there —
/// or when a numeric field does not parse, including a non-positive `rank`;
/// [`Error::MissingInformation`] when a
/// required column (`SpecId`, `ScanNr`, `Label`, `Peptide`, `Proteins`,
/// `retentiontime`, `ExpMass`, `CalcMass`, `FileName` and the configured score)
/// is absent; [`Error::InvalidValue`] when a limit is reached or
/// `sage_annotation` is requested without a path. The source turns rank zero
/// into the rank `-1`; this reader reports a parse error.
///
/// The source reads `FileName` optionally but then looks the resulting name up
/// with `std::map::at`, which throws `std::out_of_range` when no `FileName`
/// column exists; this requires the column instead. The caller may catch the
/// source exception. See the support document.
pub fn read(reader: impl BufRead, options: &ReadOptions) -> Result<PinDocument> {
    if options.sage_annotation {
        return Err(invalid(
            "Sage annotation needs the .pin path to find its sibling files",
        ));
    }
    let csv = CsvFile::from_reader(reader, &tab_options(options)?)?;
    read_csv(&csv, options, None)
}

/// Read a `.pin` file from `path`.
///
/// With [`ReadOptions::sage_annotation`] the sibling Sage files are read too:
/// `<stem>.tsv` for the `spectrum_q` column that filters rows, and
/// `<prefix>matched_fragments.sage.tsv` for the peak annotations attached to
/// each hit, where `<prefix>` is `path` with the trailing `results.sage.pin`
/// removed.
///
/// # Errors
///
/// As [`read`], plus [`Error::Io`] when a file cannot be opened and
/// [`Error::InvalidValue`] when `sage_annotation` is set and the path is too
/// short to carry the expected Sage suffixes. The source computes both sibling
/// paths with unchecked `size() - N` subtraction on the path string, so a short
/// path can wrap the requested substring length; `substr` clamps that length
/// and yields an incorrect sibling name. A positive byte cutoff can also split
/// a UTF-8 character when the expected ASCII suffix is absent.
pub fn load(path: impl AsRef<Path>, options: &ReadOptions) -> Result<PinDocument> {
    let path = path.as_ref();
    let csv = CsvFile::from_path(path, &tab_options(options)?)?;
    let sage = if options.sage_annotation {
        Some(SageSiblings::read(path, options)?)
    } else {
        None
    };
    read_csv(&csv, options, sage.as_ref())
}

fn tab_options(options: &ReadOptions) -> Result<csv::ReadOptions> {
    if options.max_rows == 0 || options.max_columns == 0 || options.max_bytes == 0 {
        return Err(invalid("percolator limits must be positive"));
    }
    let cap = csv::Limits::default();
    let limits = csv::Limits {
        max_input_bytes: cap.max_input_bytes.min(options.max_bytes),
        max_storage_bytes: cap.max_storage_bytes.min(options.max_bytes),
        max_line_bytes: cap.max_line_bytes.min(options.max_bytes),
        ..cap
    };
    Ok(csv::ReadOptions {
        separator: b'\t',
        item_enclosed: false,
        first_n: -1,
        limits,
    })
}

struct SageSiblings {
    spectrum_q: Vec<f64>,
    annotations: BTreeMap<i32, Vec<crate::identification::PeakAnnotation>>,
}
impl SageSiblings {
    fn read(path: &Path, options: &ReadOptions) -> Result<Self> {
        let name = path
            .to_str()
            .ok_or_else(|| invalid(".pin path must be UTF-8 to derive the Sage siblings"))?;
        // The source slices the path by byte offset; these strip the expected
        // suffixes instead, so neither a short nor a non-ASCII path can split a
        // character or wrap the subtraction.
        let stem = name
            .strip_suffix("pin")
            .ok_or_else(|| invalid("Sage annotation expects a .pin path"))?;
        let results = name
            .strip_suffix("results.sage.pin")
            .ok_or_else(|| invalid("Sage annotation expects a results.sage.pin path"))?;
        let tsv = CsvFile::from_path(format!("{stem}tsv"), &tab_options(options)?)?;
        let annotations = CsvFile::from_path(
            format!("{results}matched_fragments.sage.tsv"),
            &tab_options(options)?,
        )?;

        let header = Header::new(tsv.row(0)?.1)?;
        let position = header.position("spectrum_q")?;
        let mut spectrum_q = Vec::new();
        for row in 1..tsv.row_count() {
            if row > options.max_rows {
                return Err(invalid("Sage tsv row limit exceeded"));
            }
            let fields = tsv.row(row)?.1;
            spectrum_q.push(row_f64(&fields, position, row + 1, "spectrum_q")?);
        }

        let header = Header::new(annotations.row(0)?.1)?;
        let psm = header.position("psm_id")?;
        let kind = header.position("fragment_type")?;
        let ordinal = header.position("fragment_ordinals")?;
        let charge = header.position("fragment_charge")?;
        let intensity = header.position("fragment_intensity")?;
        let mz = header.position("fragment_mz_experimental")?;
        let mut mapping: BTreeMap<i32, Vec<crate::identification::PeakAnnotation>> =
            BTreeMap::new();
        let mut budget = Budget {
            bytes: options.max_bytes,
        };
        for row in 1..annotations.row_count() {
            if row > options.max_rows {
                return Err(invalid("Sage annotation row limit exceeded"));
            }
            let fields = annotations.row(row)?.1;
            let line = row + 1;
            let annotation = format!(
                "{}{}",
                field(&fields, kind, line, "fragment_type")?,
                field(&fields, ordinal, line, "fragment_ordinals")?
            );
            budget.spend(annotation.len().saturating_add(64))?;
            mapping
                .entry(row_i32(&fields, psm, line, "psm_id")?)
                .or_default()
                .push(crate::identification::PeakAnnotation {
                    mz: row_f64(&fields, mz, line, "fragment_mz_experimental")?,
                    intensity: row_f64(&fields, intensity, line, "fragment_intensity")?,
                    charge: row_i32(&fields, charge, line, "fragment_charge")?,
                    annotation,
                });
        }
        Ok(Self {
            spectrum_q,
            annotations: mapping,
        })
    }
}

fn read_csv(
    csv: &CsvFile,
    options: &ReadOptions,
    sage: Option<&SageSiblings>,
) -> Result<PinDocument> {
    if options.max_rows == 0 || options.max_columns == 0 || options.max_bytes == 0 {
        return Err(invalid("percolator limits must be positive"));
    }
    let mut document = PinDocument::default();
    if csv.row_count() == 0 {
        return Err(parse(0, ".pin file has no header line"));
    }
    if csv.row_count().saturating_sub(1) > options.max_rows {
        return Err(invalid("percolator row limit exceeded"));
    }
    let header = Header::new(csv.row(0)?.1)?;
    if header.names.len() > options.max_columns {
        return Err(invalid("percolator column limit exceeded"));
    }

    let spec_id_at = header.position("SpecId")?;
    let scan_at = header.position("ScanNr")?;
    let label_at = header.position("Label")?;
    let peptide_at = header.position("Peptide")?;
    let proteins_at = header.position("Proteins")?;
    let rt_at = header.position("retentiontime")?;
    let experimental_at = header.position("ExpMass")?;
    let calculated_at = header.position("CalcMass")?;
    let file_at = header.position("FileName")?;
    let score_at = header.position(&options.score_name)?;
    let rank_at = header.optional("rank");
    let mobility_at = header.optional("ion_mobility");
    let other_charge_at = header.optional("z=other");
    let charges = charge_columns(&header);

    let mut extra = BTreeSet::new();
    for name in &options.extra_scores {
        if header.optional(name).is_some() {
            extra.insert(name.clone());
        } else {
            document.warnings.push(format!(
                "Extra score: {name} not found in Percolator input file."
            ));
        }
    }

    let mut budget = Budget {
        bytes: options.max_bytes,
    };
    let mut filename_index: BTreeMap<String, usize> = BTreeMap::new();
    let mut current_spec_id: Option<String> = None;
    for row in 1..csv.row_count() {
        let line = row + 1;
        if let Some(sage) = sage {
            let q = sage
                .spectrum_q
                .get(row - 1)
                .copied()
                .ok_or_else(|| parse(line, "Sage tsv has fewer rows than the .pin file"))?;
            if q > options.spectrum_q_threshold {
                continue;
            }
        }
        let fields = csv.row(row)?.1;
        if fields.len() != header.names.len() {
            // A surplus field is what a tab-separated trailing protein list
            // looks like to a reader that expects a rectangular table, which
            // is what the source expects; the hint names that case rather than
            // leaving the caller to guess.
            let hint = if fields.len() > header.names.len()
                && header.names.last().map(String::as_str) == Some(PROTEINS_COLUMN)
            {
                "; a multi-accession Proteins column produces exactly this"
            } else {
                ""
            };
            return Err(parse(
                line,
                format!(
                    "line {row} does not have the same number of columns as the pin_header ({} vs {}){hint}",
                    fields.len(),
                    header.names.len()
                ),
            ));
        }

        let raw_file_name = field(&fields, file_at, line, "FileName")?.to_owned();
        let merge_index = match filename_index.get(&raw_file_name) {
            Some(index) => *index,
            None => {
                budget.spend(raw_file_name.len().saturating_add(64))?;
                document.filenames.push(raw_file_name.clone());
                let index = document.filenames.len() - 1;
                filename_index.insert(raw_file_name, index);
                index
            }
        };

        let spec_id = field(&fields, spec_id_at, line, "SpecId")?.to_owned();
        let mobility = match mobility_at {
            Some(at) => row_f64(&fields, at, line, "ion_mobility")?,
            None => 0.0,
        };
        let scan_number = field(&fields, scan_at, line, "ScanNr")?.to_owned();
        // The source compares against a string that starts empty, so a first
        // row with an empty SpecId reads past the end of the empty vector. A
        // new identification is opened whenever none is current.
        if current_spec_id.as_deref() != Some(spec_id.as_str()) {
            let mut identification = PeptideIdentification {
                higher_score_better: options.higher_score_better,
                score_type: options.score_name.clone(),
                ..Default::default()
            };
            budget.spend(spec_id.len().saturating_add(256))?;
            identification.metadata.insert(
                ID_MERGE_INDEX.into(),
                i64::try_from(merge_index)
                    .map_err(|_| invalid("percolator id_merge_index exceeds i64"))?
                    .into(),
            );
            // Search engines typically write minutes.
            identification.rt = Some(finite(
                row_f64(&fields, rt_at, line, "retentiontime")? * 60.0,
                "retentiontime",
            )?);
            identification
                .metadata
                .insert("PinSpecId".into(), spec_id.clone().into());
            if mobility > 0.0 {
                identification
                    .metadata
                    .insert(IM.into(), MetaValue::try_from(mobility)?);
            }
            identification.set_spectrum_reference(scan_number.clone());
            document.peptide_identifications.push(identification);
            current_spec_id = Some(spec_id.clone());
        }
        let identification = document
            .peptide_identifications
            .last_mut()
            .ok_or_else(|| parse(line, "no current peptide identification"))?;

        let experimental = row_f64(&fields, experimental_at, line, "ExpMass")?;
        let calculated = row_f64(&fields, calculated_at, line, "CalcMass")?;
        let score = row_f64(&fields, score_at, line, &options.score_name)?;
        let label = row_i32(&fields, label_at, line, "Label")?;
        let rank = match rank_at {
            Some(at) => row_i32(&fields, at, line, "rank")?,
            None => 1,
        };
        // The column carries the whole i32 range, so the decrement is checked
        // before it is narrowed: i32::MIN would otherwise overflow the
        // subtraction (a debug panic, and i32::MAX after a release wrap).
        let rank = rank
            .checked_sub(1)
            .and_then(|rank| u32::try_from(rank).ok())
            .ok_or_else(|| parse(line, "percolator rank must be one or greater"))?;

        let mut charge = 0;
        for (at, value) in &charges {
            if field(&fields, *at, line, "charge")? == "1" {
                charge = *value;
                break;
            }
        }
        if charge == 0 {
            if let Some(at) = other_charge_at {
                charge = row_i32(&fields, at, line, "z=other")?;
            }
        }
        if charge != 0 {
            identification.mz = Some(finite(
                experimental / f64::from(charge).abs() + PROTON_MASS_U,
                "m/z",
            )?);
        }

        let accessions: Vec<&str> = field(&fields, proteins_at, line, "Proteins")?
            .split(';')
            .collect();
        if accessions.len() > options.max_columns {
            return Err(invalid("percolator accession list limit exceeded"));
        }
        let mut status = if label == 1 {
            TargetDecoyType::Target
        } else {
            TargetDecoyType::Decoy
        };
        if !options.decoy_prefix.is_empty() {
            let mut decoy = false;
            let mut target = false;
            for accession in &accessions {
                if accession.starts_with(options.decoy_prefix.as_str()) {
                    decoy = true;
                } else {
                    target = true;
                }
                if decoy && target {
                    break;
                }
            }
            status = match (target, decoy) {
                (true, true) => TargetDecoyType::TargetAndDecoy,
                (true, false) => TargetDecoyType::Target,
                _ => TargetDecoyType::Decoy,
            };
        }

        // Percolator writes terminal annotations as "[+42]-SEQ" and "SEQ-[+1]";
        // the crate's sequence parser reads the '.' spelling.
        let peptide = field(&fields, peptide_at, line, "Peptide")?
            .replace("]-", "].")
            .replace("-[", ".[");
        budget.spend(peptide.len().saturating_add(256))?;
        let sequence = AASequence::parse(&peptide)?;
        let mut hit = PeptideHit::new(score, rank, charge, sequence)?;
        hit.set_target_decoy_type(status);

        for name in &extra {
            let at = header.position(name)?;
            let mut value = field(&fields, at, line, name)?.to_owned();
            if name == "ln(-poisson)" && value == "inf" {
                value = "3.5".into();
            }
            budget.spend(name.len().saturating_add(value.len()).saturating_add(64))?;
            hit.metadata.insert(name.clone(), value.into());
        }
        if let Some(sage) = sage {
            let q = sage
                .spectrum_q
                .get(row - 1)
                .copied()
                .ok_or_else(|| parse(line, "Sage tsv has fewer rows than the .pin file"))?;
            hit.metadata
                .insert("spectrum_q".into(), MetaValue::try_from(q)?);
        }
        hit.metadata.insert(
            "DeltaMass".into(),
            MetaValue::try_from(finite(experimental - calculated, "DeltaMass")?)?,
        );
        if let Some(sage) = sage {
            if let Ok(id) = spec_id.trim().parse::<i32>() {
                if let Some(annotations) = sage.annotations.get(&id) {
                    hit.peak_annotations = annotations.clone();
                }
            }
        }
        for accession in accessions {
            budget.spend(accession.len().saturating_add(64))?;
            hit.evidences.push(PeptideEvidence {
                protein_accession: accession.to_owned(),
                ..Default::default()
            });
        }
        identification.hits.push(hit);
    }
    Ok(document)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_feature_set_rejects_an_unbounded_charge_span() {
        assert!(standard_feature_set(i32::MIN, i32::MAX).is_err());
        assert!(
            standard_feature_set(4, 2)
                .unwrap()
                .iter()
                .all(|c| !c.starts_with("charge"))
        );
    }

    #[test]
    fn scan_number_uses_the_last_match_and_rejects_overflow() {
        let pattern = ScanNumberPattern::for_native_id("scan=1");
        assert_eq!(pattern.extract("scan=1 scan=42"), Some(42));
        assert_eq!(pattern.extract("index=7"), None);
        assert_eq!(
            ScanNumberPattern::for_native_id("9").extract("99999999999999"),
            None
        );
    }

    #[test]
    fn scan_identifier_strips_whitespace_from_non_ascii_references() {
        let mut identification = PeptideIdentification::new();
        identification.set_spectrum_reference(" 日本 語 ");
        assert_eq!(scan_identifier(&identification, 3), "日本語");
    }

    #[test]
    fn count_enzymatic_counts_characters_not_bytes() {
        assert_eq!(count_enzymatic("KKK", "trypsin"), 2);
        // A non-ASCII residue letter is one character, not several bonds.
        assert_eq!(count_enzymatic("Kä", "trypsin"), 1);
        assert_eq!(count_enzymatic("", "trypsin"), 0);
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Mascot XML search-result reader, ported from `FORMAT/MascotXMLFile.h` and
//! the `FORMAT/HANDLERS/MascotXMLHandler.h` SAX handler it drives.
//!
//! A Mascot XML document carries a `<header>`, the search parameters, the
//! protein/peptide `<hits>`, the `<unassigned>` peptides and optionally a
//! `<queries>` section. Retention times are **not** part of the format in
//! general: they are recovered from each query's `<pep_scan_title>` or
//! `<StringTitle>` through a caller-supplied
//! [`SpectrumTitleLookup`](crate::format::mascot_xml::SpectrumTitleLookup),
//! which is why the source signature takes that helper separately.
//!
//! Entry points are
//! [`MascotXmlFile`](crate::format::mascot_xml::MascotXmlFile) and the free
//! functions [`read`](crate::format::mascot_xml::read) and
//! [`load`](crate::format::mascot_xml::load). Everything this module does and
//! does not reproduce is listed in `docs/MASCOT_XML_SUPPORT.md`.

use crate::chemistry::{AASequence, ModificationsDB, ProteaseDB};
use crate::identification::{
    FlankingResidue, PeakMassType, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification, SearchParameters,
};
use crate::metadata::{CompletionTime, MetaValue};
use crate::{Error, MSExperiment, Result};
use quick_xml::Reader;
use quick_xml::events::Event;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Meta key under which the handler stores each hit's Mascot E-value.
pub const EVALUE_KEY: &str = "EValue";
/// Meta key for the Mascot homology threshold recorded on a peptide hit.
pub const HOMOLOGY_THRESHOLD_KEY: &str = "homology_threshold";
/// Meta key for the Mascot identity threshold recorded on a peptide hit.
pub const IDENTITY_THRESHOLD_KEY: &str = "identity_threshold";
/// Score type the handler stamps on every protein and peptide identification.
pub const SCORE_TYPE: &str = "Mascot";

/// Resource ceilings for one document. Counts include elements and text that
/// are read and then discarded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum decoded document size.
    pub max_bytes: usize,
    /// Maximum XML events.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum `<NumQueries>`, which sizes the identification vector.
    pub max_queries: usize,
    /// Maximum peptide hits per query and protein hits per document.
    pub max_hits: usize,
    /// Maximum characters accumulated for one element's text.
    pub max_text_bytes: usize,
    /// Maximum entries in one modification list.
    pub max_modifications: usize,
}
impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 256 * 1024 * 1024,
            max_events: 20_000_000,
            max_depth: 64,
            max_queries: 5_000_000,
            max_hits: 100_000,
            max_text_bytes: 1024 * 1024,
            max_modifications: 10_000,
        }
    }
}
impl ReadLimits {
    fn validate(&self) -> Result<()> {
        let ceiling = Self::default();
        if self.max_bytes > ceiling.max_bytes
            || self.max_events > ceiling.max_events
            || self.max_depth > ceiling.max_depth
            || self.max_queries > ceiling.max_queries
            || self.max_hits > ceiling.max_hits
            || self.max_text_bytes > ceiling.max_text_bytes
            || self.max_modifications > ceiling.max_modifications
        {
            return Err(bad("Mascot XML limits exceed their hard ceilings"));
        }
        Ok(())
    }
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn parse(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// `StringUtils::trim`: only space, tab, CR and LF are stripped.
fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
}
/// `StringUtils::toDouble` restricted to finite values.
///
/// One leading `+` is consumed and the remainder goes to `std::from_chars`,
/// which refuses a *second* `+` even though Rust's own parser accepts it; an
/// overflowing decimal literal is `result_out_of_range` there and a conversion
/// error here, not an infinity.
fn number(text: &str, label: &str) -> Result<f64> {
    let text = trim(text);
    let body = text.strip_prefix('+').unwrap_or(text);
    let value = if body.starts_with('+') {
        None
    } else {
        body.parse::<f64>().ok().filter(|v| v.is_finite())
    };
    value.ok_or_else(|| parse(format!("Could not convert {text:?} to a finite {label}")))
}
/// `StringUtils::toInt32`, which likewise refuses a second `+`.
fn integer(text: &str, label: &str) -> Result<i32> {
    let text = trim(text);
    let body = text.strip_prefix('+').unwrap_or(text);
    let value = if body.starts_with('+') {
        None
    } else {
        body.parse::<i32>().ok()
    };
    value.ok_or_else(|| parse(format!("Could not convert {text:?} to an integer {label}")))
}
/// `StringUtils::split(s, c, out)` semantics.
fn split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// Per-spectrum meta data the title lookup can supply, the ported subset of
/// `SpectrumMetaDataLookup::SpectrumMetaData`.
///
/// Every floating-point field defaults to `None`, which stands for the source's
/// quiet NaN: assigning it to a [`PeptideIdentification`] leaves the value
/// unset rather than storing zero.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumMetaData {
    /// Retention time in seconds.
    pub rt: Option<f64>,
    /// Precursor mass-to-charge ratio.
    pub precursor_mz: Option<f64>,
    /// Precursor charge; zero when unknown, as in the source.
    pub precursor_charge: i32,
    /// MS level; zero when unknown, as in the source.
    pub ms_level: u32,
    /// Scan number, or `None` for the source's `-1` sentinel. The source reads
    /// it with `toInt32`, so the width is 32 bits and a longer digit run is a
    /// conversion failure rather than a large scan number.
    pub scan_number: Option<i32>,
    /// Spectrum native identifier.
    pub native_id: String,
}

/// One of the spectrum-reference formats `MascotXMLFile::initializeLookup`
/// registers.
///
/// The source registers Boost regular expressions. This crate has no regular
/// expression engine and may not add a dependency, so the three default
/// formats are hand-coded and a caller-supplied `scan_regex` is refused; see
/// `docs/MASCOT_XML_SUPPORT.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TitleReferenceFormat {
    /// `[Ss]can( [Nn]umber)?s?[=:]? *(?<SCAN>\d+)`: `scan=818`,
    /// `Spectrum136 scans:712,`, `6860: Scan 10668 (rt=5380.57)` and
    /// `Scan Number: 1460` all yield their scan number.
    ScanNumber,
    /// `\.(?<SCAN>\d+)\.\d+\.(?<CHARGE>\d+)(\.dta)?`: a DTA file name such as
    /// `/path/to/FTAC05_13.673.673.2.dta` yields scan 673 and charge 2.
    DtaFileName,
    /// `^(?<MZ>\d+(\.\d+)?)_(?<RT>\d+(\.\d+)?)`: a title that starts with
    /// precursor m/z and retention time joined by an underscore, as
    /// `MascotGenericFile` writes them.
    MzThenRt,
}

/// What one reference format extracted from a title.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct TitleMatch {
    rt: Option<f64>,
    mz: Option<f64>,
    charge: Option<i32>,
    scan: Option<i32>,
}

/// The ported subset of `SpectrumMetaDataLookup` that `MascotXMLFile` needs.
///
/// A default-constructed lookup holds no spectra and no reference formats, so
/// [`Self::spectrum_meta_data`] extracts nothing and the retention times of
/// every identification stay unset — which is exactly what the upstream class
/// test does and why its expected output carries no `RT` attribute. Call
/// [`MascotXmlFile::initialize_lookup`] to register the default formats.
#[derive(Clone, Debug, Default)]
pub struct SpectrumTitleLookup {
    spectra: Vec<SpectrumMetaData>,
    by_scan: BTreeMap<i32, usize>,
    by_native_id: BTreeMap<String, usize>,
    rts: Vec<(f64, usize)>,
    formats: Vec<TitleReferenceFormat>,
    /// Half-width of the retention-time window a title-derived RT may match,
    /// the source `SpectrumLookup::rt_tolerance` default.
    pub rt_tolerance: f64,
}

impl SpectrumTitleLookup {
    /// An empty lookup: no spectra, no formats, `rt_tolerance` 0.01 s.
    pub fn new() -> Self {
        Self {
            rt_tolerance: 0.01,
            ..Default::default()
        }
    }

    /// Record every spectrum's meta data, replacing anything read before.
    ///
    /// The scan number comes from the source default expression
    /// `=(?<SCAN>\d+)$`: the digits following the last `=` at the end of the
    /// native ID. A native ID with no such suffix contributes no scan-number
    /// entry, which the source reports as a warning and this port returns.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when the experiment
    /// holds more spectra than [`ReadLimits::max_queries`] allows.
    pub fn read_spectra(&mut self, experiment: &MSExperiment) -> Result<Vec<String>> {
        let limits = ReadLimits::default();
        if experiment.spectra.len() > limits.max_queries {
            return Err(bad("Mascot XML spectrum lookup limit exceeded"));
        }
        let mut spectra = Vec::new();
        spectra
            .try_reserve_exact(experiment.spectra.len())
            .map_err(|_| bad("spectrum lookup allocation failed"))?;
        let mut warnings = Vec::new();
        for spectrum in &experiment.spectra {
            let scan_number = trailing_scan_number(&spectrum.native_id);
            if scan_number.is_none() {
                warnings.push(format!(
                    "Warning: Could not extract scan number from spectrum native ID '{}' using regular expression '=(?<SCAN>\\d+)$'. Look-up by scan number may not work properly.",
                    spectrum.native_id
                ));
            }
            spectra.push(SpectrumMetaData {
                rt: Some(spectrum.rt),
                precursor_mz: spectrum.precursors.first().map(|p| p.mz),
                precursor_charge: spectrum.precursors.first().map_or(0, |p| p.charge),
                ms_level: spectrum.ms_level,
                scan_number,
                native_id: spectrum.native_id.clone(),
            });
        }
        self.by_scan.clear();
        self.by_native_id.clear();
        self.rts.clear();
        for (index, meta) in spectra.iter().enumerate() {
            if let Some(scan) = meta.scan_number {
                self.by_scan.insert(scan, index);
            }
            self.by_native_id.insert(meta.native_id.clone(), index);
            if let Some(rt) = meta.rt {
                self.rts.push((rt, index));
            }
        }
        self.rts
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        self.spectra = spectra;
        Ok(warnings)
    }

    /// Append a reference format, as `addReferenceFormat` does.
    pub fn add_reference_format(&mut self, format: TitleReferenceFormat) {
        self.formats.push(format);
    }

    /// Whether no spectra have been read. Mirrors `SpectrumLookup::empty`,
    /// which tests the spectrum count and not the format list.
    pub fn is_empty(&self) -> bool {
        self.spectra.is_empty()
    }

    /// Number of spectra recorded.
    pub fn len(&self) -> usize {
        self.spectra.len()
    }

    /// The registered reference formats, in the order they will be tried.
    pub fn reference_formats(&self) -> &[TitleReferenceFormat] {
        &self.formats
    }

    /// Recorded meta data of one spectrum by index.
    pub fn spectrum(&self, index: usize) -> Option<&SpectrumMetaData> {
        self.spectra.get(index)
    }

    /// Look up a spectrum reference, as `getSpectrumMetaData(ref, meta, flags)`.
    ///
    /// The first registered format that matches `title` wins. Fields the match
    /// itself carries are used directly; if either requested field is still
    /// missing, the matched scan number, native ID or retention time selects a
    /// recorded spectrum and its whole meta data is returned. A title matching
    /// no format yields `Ok(SpectrumMetaData::default())` — all fields unset,
    /// no error — which is what the source does, since its loop simply ends.
    ///
    /// `want_mz` mirrors the handler's conditional `MDF_PRECURSORMZ` flag: the
    /// m/z is only requested when the identification has none yet.
    ///
    /// # Errors
    ///
    /// Returns [`MissingInformation`](crate::Error::MissingInformation) when a
    /// format matched but carried nothing usable, and
    /// [`InvalidRange`](crate::Error::InvalidRange) when the referenced
    /// spectrum is not among those read — the two conditions for which the
    /// source throws `MissingInformation` and `ElementNotFound`. The handler
    /// catches both and reports them as warnings.
    pub fn spectrum_meta_data(&self, title: &str, want_mz: bool) -> Result<SpectrumMetaData> {
        for &format in &self.formats {
            // A format that matched but whose captured value does not convert is
            // an error, not a miss: the source returns after the first matching
            // expression, so the later formats are never tried.
            let Some(matched) = match_title(title, format)? else {
                continue;
            };
            let mut meta = SpectrumMetaData::default();
            let mut need_rt = true;
            let mut need_mz = want_mz;
            if let Some(rt) = matched.rt {
                meta.rt = Some(rt);
                need_rt = false;
            }
            if let Some(mz) = matched.mz {
                meta.precursor_mz = Some(mz);
                need_mz = false;
            }
            if let Some(charge) = matched.charge {
                meta.precursor_charge = charge;
            }
            if let Some(scan) = matched.scan {
                meta.scan_number = Some(scan);
            }
            if !need_rt && !need_mz {
                return Ok(meta);
            }
            let index = self.find_by_match(title, &matched)?;
            return Ok(self.spectra[index].clone());
        }
        Ok(SpectrumMetaData::default())
    }

    fn find_by_match(&self, title: &str, matched: &TitleMatch) -> Result<usize> {
        if let Some(scan) = matched.scan {
            return self.by_scan.get(&scan).copied().ok_or_else(|| {
                Error::InvalidRange(format!("spectrum with scan number {scan} not found"))
            });
        }
        if let Some(rt) = matched.rt {
            return self.find_by_rt(rt);
        }
        Err(Error::MissingInformation(format!(
            "Unexpected format of spectrum reference '{title}'. The reference format matched, but no usable information could be extracted."
        )))
    }

    fn find_by_rt(&self, rt: f64) -> Result<usize> {
        let mut best: Option<(f64, usize)> = None;
        for &(value, index) in &self.rts {
            let difference = (value - rt).abs();
            if difference <= self.rt_tolerance
                && best.is_none_or(|(previous, _)| difference < previous)
            {
                best = Some((difference, index));
            }
        }
        best.map(|(_, index)| index)
            .ok_or_else(|| Error::InvalidRange(format!("spectrum with RT {rt} not found")))
    }
}

/// The byte offsets at which Boost's `^` matches.
///
/// `initializeLookup` compiles its expressions with the default perl flags, so
/// `^` becomes `syntax_element_start_line` rather than `..._buffer_start`
/// (`basic_regex_parser.hpp`, `syntax_caret`): it matches at the buffer start
/// *and* after every line separator. For `char` Boost's separators are `\n`,
/// `\r` and `\f`, and the position between a `\r` and a `\n` is not a line
/// start (`perl_matcher_common.hpp`, `match_start_line`).
fn line_starts(text: &str) -> impl Iterator<Item = usize> + '_ {
    let bytes = text.as_bytes();
    std::iter::once(0).chain((1..=bytes.len()).filter(move |&index| {
        let previous = bytes[index - 1];
        matches!(previous, b'\n' | b'\r' | 0x0c)
            && !(previous == b'\r' && bytes.get(index) == Some(&b'\n'))
    }))
}

/// Whether Boost's `$` matches at byte offset `at`: the buffer end, or a
/// position holding a line separator.
fn is_line_end(bytes: &[u8], at: usize) -> bool {
    match bytes.get(at) {
        None => true,
        Some(&byte) => matches!(byte, b'\n' | b'\r' | 0x0c),
    }
}

/// The digits following the last `=` at a line end of `native_id`, the source
/// `SpectrumLookup::default_scan_regexp` `=(?<SCAN>\d+)$`.
///
/// `$` is a line anchor, so a native ID with a trailing annotation line still
/// yields the scan number of its first line, and the token iterator takes the
/// last match. The conversion is `toInt32`, whose failure the source answers
/// with its `-1` sentinel — reported here as `None`, which records no
/// scan-number entry and produces the source's warning.
fn trailing_scan_number(native_id: &str) -> Option<i32> {
    let bytes = native_id.as_bytes();
    let mut found = None;
    for (index, byte) in bytes.iter().enumerate() {
        if *byte != b'=' {
            continue;
        }
        let digits_at = index + 1;
        let mut end = digits_at;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
        if end > digits_at && is_line_end(bytes, end) {
            found = native_id.get(digits_at..end);
        }
    }
    found?.parse::<i32>().ok()
}

/// A hand-coded equivalent of one registered Boost expression.
///
/// `Ok(None)` is "the expression did not match"; an error is "it matched, but a
/// captured value did not convert", which the source reports through the
/// handler's `catch (...)` without trying another format.
fn match_title(title: &str, format: TitleReferenceFormat) -> Result<Option<TitleMatch>> {
    match format {
        TitleReferenceFormat::ScanNumber => match_scan_number(title),
        TitleReferenceFormat::DtaFileName => match_dta_name(title),
        TitleReferenceFormat::MzThenRt => match_mz_then_rt(title),
    }
}

/// `[Ss]can( [Nn]umber)?s?[=:]? *(?<SCAN>\d+)`, leftmost match.
fn match_scan_number(title: &str) -> Result<Option<TitleMatch>> {
    let bytes = title.as_bytes();
    for start in 0..bytes.len() {
        if !matches!(bytes[start], b'S' | b's') {
            continue;
        }
        let mut at = start + 1;
        if !title.get(at..).is_some_and(|rest| rest.starts_with("can")) {
            continue;
        }
        at += 3;
        // optional " Number" / " number"
        if title
            .get(at..)
            .is_some_and(|rest| rest.starts_with(" Number") || rest.starts_with(" number"))
        {
            at += 7;
        }
        // optional trailing 's'
        if bytes.get(at) == Some(&b's') {
            at += 1;
        }
        // optional '=' or ':'
        if matches!(bytes.get(at), Some(b'=') | Some(b':')) {
            at += 1;
        }
        // any number of spaces
        while bytes.get(at) == Some(&b' ') {
            at += 1;
        }
        let digits_at = at;
        while bytes.get(at).is_some_and(u8::is_ascii_digit) {
            at += 1;
        }
        if at == digits_at {
            continue;
        }
        let digits = title.get(digits_at..at).unwrap_or("");
        return Ok(Some(TitleMatch {
            scan: Some(scan_number(digits)?),
            ..Default::default()
        }));
    }
    Ok(None)
}

/// A matched `?<SCAN>` group, converted as `getSpectrumMetaData` does.
///
/// # Errors
///
/// Returns [`Parse`](crate::Error::Parse) for a digit run too long for the
/// source's `toInt32`, which throws rather than trying another format.
fn scan_number(digits: &str) -> Result<i32> {
    digits.parse::<i32>().map_err(|_| {
        parse(format!(
            "Could not convert {digits:?} to an integer scan number"
        ))
    })
}

/// `\.(?<SCAN>\d+)\.\d+\.(?<CHARGE>\d+)(\.dta)?`, leftmost match.
fn match_dta_name(title: &str) -> Result<Option<TitleMatch>> {
    let bytes = title.as_bytes();
    let digits = |from: usize| -> (usize, usize) {
        let mut at = from;
        while bytes.get(at).is_some_and(u8::is_ascii_digit) {
            at += 1;
        }
        (from, at)
    };
    for start in 0..bytes.len() {
        if bytes[start] != b'.' {
            continue;
        }
        let (scan_from, scan_to) = digits(start + 1);
        if scan_to == scan_from || bytes.get(scan_to) != Some(&b'.') {
            continue;
        }
        let (middle_from, middle_to) = digits(scan_to + 1);
        if middle_to == middle_from || bytes.get(middle_to) != Some(&b'.') {
            continue;
        }
        let (charge_from, charge_to) = digits(middle_to + 1);
        if charge_to == charge_from {
            continue;
        }
        let scan = scan_number(title.get(scan_from..scan_to).unwrap_or(""))?;
        // `getSpectrumMetaData` converts the charge group only when
        // `MDF_PRECURSORCHARGE` is requested, which `MascotXMLHandler` never
        // does, so the source never converts it at all and a charge that does
        // not fit is not an error. It is recorded here when it converts and
        // left unset otherwise; either way the spectrum look-up that a
        // scan-only match falls through to replaces the whole record.
        let charge = title
            .get(charge_from..charge_to)
            .and_then(|digits| digits.parse::<i32>().ok());
        return Ok(Some(TitleMatch {
            scan: Some(scan),
            charge,
            ..Default::default()
        }));
    }
    Ok(None)
}

/// `^(?<MZ>\d+(\.\d+)?)_(?<RT>\d+(\.\d+)?)`, anchored at a line start.
///
/// Boost's `^` is a line anchor by default, so a wrapped title whose *second*
/// line starts with the m/z-underscore-RT pair matches too.
fn match_mz_then_rt(title: &str) -> Result<Option<TitleMatch>> {
    let bytes = title.as_bytes();
    let unsigned = |from: usize| -> Option<usize> {
        let mut at = from;
        while bytes.get(at).is_some_and(u8::is_ascii_digit) {
            at += 1;
        }
        if at == from {
            return None;
        }
        if bytes.get(at) == Some(&b'.') {
            let fraction = at + 1;
            let mut end = fraction;
            while bytes.get(end).is_some_and(u8::is_ascii_digit) {
                end += 1;
            }
            if end > fraction {
                at = end;
            }
        }
        Some(at)
    };
    for start in line_starts(title) {
        let Some(mz_end) = unsigned(start) else {
            continue;
        };
        if bytes.get(mz_end) != Some(&b'_') {
            continue;
        }
        let Some(rt_end) = unsigned(mz_end + 1) else {
            continue;
        };
        let mz = coordinate(title.get(start..mz_end).unwrap_or(""), "precursor m/z")?;
        let rt = coordinate(
            title.get(mz_end + 1..rt_end).unwrap_or(""),
            "retention time",
        )?;
        return Ok(Some(TitleMatch {
            mz: Some(mz),
            rt: Some(rt),
            ..Default::default()
        }));
    }
    Ok(None)
}

/// A matched `?<MZ>` or `?<RT>` group, converted as `getSpectrumMetaData` does.
///
/// # Errors
///
/// Returns [`Parse`](crate::Error::Parse) when the digits overflow a `double`,
/// which is the source's `toDouble` conversion error; Rust's own parser would
/// return an infinity and defeat the finite-coordinate invariant.
fn coordinate(digits: &str, label: &str) -> Result<f64> {
    digits
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| parse(format!("Could not convert {digits:?} to a finite {label}")))
}

/// What one Mascot XML document yielded.
///
/// The source `load` signature writes through three out-parameters; this is
/// the returned equivalent.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MascotXmlResult {
    /// Protein identifications belonging to the whole experiment.
    pub protein_identification: ProteinIdentification,
    /// Peptide identifications with m/z and, where a lookup supplied one, RT.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// Messages the source sends to `OPENMS_LOG_WARN` or the XML handler's
    /// non-fatal `error`/`warning` reporters.
    pub warnings: Vec<String>,
}

/// Used to load Mascot XML files.
///
/// This class loads documents that implement the schema of Mascot XML files.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MascotXmlFile;

impl MascotXmlFile {
    /// A reader. The source class carries only its `XMLFile` base state, whose
    /// schema location and version are empty for Mascot XML.
    pub fn new() -> Self {
        Self
    }

    /// Load data from a Mascot XML file.
    ///
    /// `lookup` is the helper object used to look up spectrum meta data; an
    /// empty one leaves every retention time unset.
    ///
    /// # Errors
    ///
    /// Returns [`Io`](crate::Error::Io) if the file does not exist, where the
    /// source throws `Exception::FileNotFound`;
    /// [`Parse`](crate::Error::Parse) if the file does not suit the standard,
    /// where the source throws `Exception::ParseError`; and
    /// [`MissingInformation`](crate::Error::MissingInformation) when a
    /// non-empty `lookup` was supplied and no identification received a
    /// retention time.
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        lookup: &SpectrumTitleLookup,
    ) -> Result<MascotXmlResult> {
        load(path, lookup)
    }

    /// Load data from a Mascot XML file, replacing hit sequences from a map of
    /// modified peptides identified by their `<StringTitle>`.
    ///
    /// The map is the pepXML-derived one `MascotAdapter` passes: for every
    /// query title it holds the modified sequences, and each Mascot hit whose
    /// unmodified sequence matches one of them adopts that modified sequence.
    /// A non-empty map also suppresses the fixed-modification pass on
    /// `<pep_seq>`, because the supplied sequences already carry them.
    ///
    /// # Errors
    ///
    /// See [`Self::load`].
    pub fn load_with_peptides(
        &self,
        path: impl AsRef<Path>,
        peptides: &BTreeMap<String, Vec<AASequence>>,
        lookup: &SpectrumTitleLookup,
    ) -> Result<MascotXmlResult> {
        load_with_peptides(path, peptides, lookup)
    }

    /// Initialise a helper object for looking up spectrum meta data (RT, m/z).
    ///
    /// Reads the experiment's spectra and registers the default reference
    /// formats. With raw data available, scan-number and DTA-file-name titles
    /// are recognised as well as the m/z-underscore-RT form; without it only
    /// the last is, because it needs no spectrum to resolve.
    ///
    /// `scan_regex` is the source's optional override. This port cannot compile
    /// an arbitrary expression and refuses a non-empty one; pass `None` to use
    /// the default formats, or build the lookup directly with
    /// [`SpectrumTitleLookup::add_reference_format`].
    ///
    /// # Errors
    ///
    /// Returns [`Unsupported`](crate::Error::Unsupported) for a non-empty
    /// `scan_regex`, and [`InvalidValue`](crate::Error::InvalidValue) when the
    /// experiment exceeds the spectrum ceiling.
    pub fn initialize_lookup(
        experiment: &MSExperiment,
        scan_regex: Option<&str>,
    ) -> Result<(SpectrumTitleLookup, Vec<String>)> {
        let mut lookup = SpectrumTitleLookup::new();
        let warnings = lookup.read_spectra(experiment)?;
        match scan_regex.map(str::trim).filter(|value| !value.is_empty()) {
            Some(_) => {
                return Err(Error::Unsupported(
                    "a caller-supplied Mascot scan_regex needs a regular-expression engine, which this crate does not have; register a TitleReferenceFormat instead".into(),
                ));
            }
            None => {
                // Raw data given, so a spectrum look-up is possible. Possible
                // formats and resulting scan numbers:
                //   <pep_scan_title>scan=818</pep_scan_title>            -> 818
                //   <pep_scan_title>Spectrum136 scans:712,</...>         -> 712
                //   <pep_scan_title>Spectrum3411 scans: 2975,</...>      -> 2975
                //   <...>File773 Spectrum198145 scans: 6094</...>        -> 6094
                //   <...>6860: Scan 10668 (rt=5380.57)</...>             -> 10668
                //   <pep_scan_title>Scan Number: 1460</pep_scan_title>   -> 1460
                // and, with .dta input to Mascot,
                //   <...>/path/to/FTAC05_13.673.673.2.dta</...>          -> 673
                if !lookup.is_empty() {
                    lookup.add_reference_format(TitleReferenceFormat::ScanNumber);
                    lookup.add_reference_format(TitleReferenceFormat::DtaFileName);
                }
                // A title containing RT and m/z instead of a scan number:
                //   575.848571777344_5018.0811_controllerType=0 ... scan=11515_EcoliMS2small
                lookup.add_reference_format(TitleReferenceFormat::MzThenRt);
            }
        }
        Ok((lookup, warnings))
    }
}

/// Read a Mascot XML document from a stream.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn read(reader: impl BufRead, lookup: &SpectrumTitleLookup) -> Result<MascotXmlResult> {
    read_with_options(reader, &BTreeMap::new(), lookup, &ReadLimits::default())
}

/// Read a Mascot XML document from a stream, with a modified-peptide map.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn read_with_peptides(
    reader: impl BufRead,
    peptides: &BTreeMap<String, Vec<AASequence>>,
    lookup: &SpectrumTitleLookup,
) -> Result<MascotXmlResult> {
    read_with_options(reader, peptides, lookup, &ReadLimits::default())
}

/// Read a Mascot XML document under explicit resource ceilings.
///
/// Modification names are resolved through the shared [`ModificationsDB`], as
/// the source's `ModificationsDB::getInstance()` does.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn read_with_options(
    reader: impl BufRead,
    peptides: &BTreeMap<String, Vec<AASequence>>,
    lookup: &SpectrumTitleLookup,
    limits: &ReadLimits,
) -> Result<MascotXmlResult> {
    read_with_registry(
        reader,
        peptides,
        lookup,
        limits,
        ModificationsDB::global(),
        ProteaseDB::global(),
    )
}

/// Read a Mascot XML document with caller-owned registries.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn read_with_registry(
    reader: impl BufRead,
    peptides: &BTreeMap<String, Vec<AASequence>>,
    lookup: &SpectrumTitleLookup,
    limits: &ReadLimits,
    modifications: &ModificationsDB,
    proteases: &ProteaseDB,
) -> Result<MascotXmlResult> {
    limits.validate()?;
    let text = document(reader, limits.max_bytes)?;
    let mut handler = Handler::new(peptides, lookup, limits, modifications, proteases);
    handler.run(&text)?;
    handler.finish()
}

/// Load a Mascot XML file.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn load(path: impl AsRef<Path>, lookup: &SpectrumTitleLookup) -> Result<MascotXmlResult> {
    read(BufReader::new(File::open(path)?), lookup)
}

/// Load a Mascot XML file with a modified-peptide map.
///
/// # Errors
///
/// See [`MascotXmlFile::load`].
pub fn load_with_peptides(
    path: impl AsRef<Path>,
    peptides: &BTreeMap<String, Vec<AASequence>>,
    lookup: &SpectrumTitleLookup,
) -> Result<MascotXmlResult> {
    read_with_peptides(BufReader::new(File::open(path)?), peptides, lookup)
}

/// Read at most `limit` bytes and require valid UTF-8.
///
/// The source hands the bytes to Xerces, which honours the XML declaration's
/// encoding. This port accepts UTF-8 only, which covers every Mascot export
/// seen upstream; a different encoding is an explicit error rather than
/// silently mangled text.
fn document(input: impl Read, limit: usize) -> Result<String> {
    let count = u64::try_from(limit)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| bad("Mascot XML byte limit overflows"))?;
    let mut bytes = Vec::new();
    input.take(count).read_to_end(&mut bytes)?;
    if bytes.len() > limit {
        return Err(bad("Mascot XML byte limit exceeded"));
    }
    String::from_utf8(bytes).map_err(|_| parse("Mascot XML input is not valid UTF-8"))
}

/// Element name without any namespace prefix.
///
/// The source handler matches on the Xerces qname, which for a document using
/// a default namespace — as every Mascot export does — is the local name. A
/// prefixed document would not be recognised by the source at all; stripping
/// the prefix here is strictly more permissive and is recorded as such.
fn local_name(raw: &[u8]) -> Result<&str> {
    let name =
        std::str::from_utf8(raw).map_err(|_| parse("Mascot XML element name is not UTF-8"))?;
    Ok(match name.rfind(':') {
        Some(position) => name.get(position + 1..).unwrap_or(name),
        None => name,
    })
}

struct Handler<'a> {
    protein: ProteinIdentification,
    ids: Vec<PeptideIdentification>,
    declared_queries: usize,
    protein_hit: ProteinHit,
    peptide_hit: PeptideHit,
    evidence: PeptideEvidence,
    index: usize,
    query: u32,
    search: SearchParameters,
    identifier: String,
    tags_open: Vec<String>,
    tag: String,
    buffer: String,
    major_version: String,
    remove_fixed_mods: Vec<String>,
    peptides: &'a BTreeMap<String, Vec<AASequence>>,
    lookup: &'a SpectrumTitleLookup,
    limits: ReadLimits,
    modifications: &'a ModificationsDB,
    proteases: &'a ProteaseDB,
    no_rt_error: bool,
    warnings: Vec<String>,
}

impl<'a> Handler<'a> {
    fn new(
        peptides: &'a BTreeMap<String, Vec<AASequence>>,
        lookup: &'a SpectrumTitleLookup,
        limits: &ReadLimits,
        modifications: &'a ModificationsDB,
        proteases: &'a ProteaseDB,
    ) -> Self {
        Self {
            protein: ProteinIdentification::default(),
            ids: Vec::new(),
            declared_queries: 0,
            protein_hit: ProteinHit::default(),
            peptide_hit: PeptideHit::default(),
            evidence: PeptideEvidence::default(),
            index: 0,
            query: 0,
            search: SearchParameters::default(),
            identifier: String::new(),
            tags_open: Vec::new(),
            tag: String::new(),
            buffer: String::new(),
            major_version: String::new(),
            remove_fixed_mods: Vec::new(),
            peptides,
            lookup,
            limits: *limits,
            modifications,
            proteases,
            no_rt_error: false,
            warnings: Vec::new(),
        }
    }

    fn run(&mut self, text: &str) -> Result<()> {
        let mut reader = Reader::from_str(text);
        reader.config_mut().expand_empty_elements = true;
        reader.config_mut().check_end_names = true;
        let mut events = 0usize;
        // `check_end_names` only pairs the tags it sees, so the one-complete-
        // document rule Xerces enforces is checked here: exactly one root
        // element, closed, and no character data outside it.
        let mut root_opened = false;
        let mut root_closed = false;
        loop {
            let event = reader
                .read_event()
                .map_err(|error| parse(format!("Mascot XML is not well formed: {error}")))?;
            events += 1;
            if events > self.limits.max_events {
                return Err(bad("Mascot XML event limit exceeded"));
            }
            match event {
                Event::Start(start) => {
                    let name = local_name(start.name().as_ref())?.to_owned();
                    if self.tags_open.is_empty() {
                        if root_opened {
                            return Err(parse("Mascot XML holds more than one root element"));
                        }
                        root_opened = true;
                    }
                    if self.tags_open.len() >= self.limits.max_depth {
                        return Err(bad("Mascot XML nesting depth limit exceeded"));
                    }
                    let mut attributes = BTreeMap::new();
                    for attribute in start.attributes() {
                        let attribute = attribute
                            .map_err(|error| parse(format!("Mascot XML attribute: {error}")))?;
                        let key = local_name(attribute.key.as_ref())?.to_owned();
                        let value = attribute
                            .unescape_value()
                            .map_err(|error| parse(format!("Mascot XML attribute: {error}")))?
                            .into_owned();
                        attributes.insert(key, value);
                    }
                    self.tag = name.clone();
                    self.tags_open.push(name);
                    self.start_element(&attributes)?;
                }
                Event::End(end) => {
                    let raw = end.name();
                    let name = local_name(raw.as_ref())?;
                    self.tag = trim(name).to_owned();
                    if self.tags_open.pop().is_none() {
                        return Err(parse(format!(
                            "Closing tag {} not matched by opening tag",
                            self.tag
                        )));
                    }
                    self.end_element()?;
                    self.tag.clear();
                    self.buffer.clear();
                    root_closed = self.tags_open.is_empty();
                }
                Event::Text(body) => {
                    let decoded = body
                        .xml_content()
                        .map_err(|error| parse(format!("Mascot XML text: {error}")))?;
                    if self.tags_open.is_empty() {
                        // Character data before or after the root element is not
                        // a well-formed document; only white space is allowed
                        // there.
                        if trim(&decoded).is_empty() {
                            continue;
                        }
                        return Err(parse(
                            "Mascot XML holds character data outside the root element",
                        ));
                    }
                    // Source `onCharacters` ignores text that follows a child
                    // element's end tag, because `tag_` is cleared there.
                    if self.tag.is_empty() {
                        continue;
                    }
                    self.append_text(&decoded)?;
                }
                // quick-xml does not expand references inside text: it splits
                // the text at every `&...;` and emits the reference as its own
                // event, and `BytesText::xml_content()` only decodes and
                // normalises line endings. Swallowing these in the catch-all
                // arm DELETED the reference and concatenated the surrounding
                // fragments, so `<pep_score>1&#46;5</pep_score>` read as the
                // score 15. `src/format/mzml.rs:2465` already refuses both for
                // the same reason, and a DTD is refused there too.
                // Predefined and numeric references are ORDINARY XML and must
                // be expanded, not refused: a protein description holding
                // `&amp;` is valid Mascot output. Only a reference this reader
                // cannot resolve without a DTD is refused. `src/format/mzml.rs`
                // and `src/format/imzml_handler.rs` refuse every reference
                // because the only text they consume is Base64, which never
                // contains `&`; that reasoning does not carry to this reader,
                // whose text nodes hold sequences and free-text descriptions.
                Event::GeneralRef(reference) => {
                    if root_opened && !root_closed && !self.tag.is_empty() {
                        let expanded = expand_reference(reference.as_ref())?;
                        self.append_text(&expanded)?;
                    }
                }
                Event::CData(_) => {
                    return Err(Error::Unsupported(
                        "CDATA sections in Mascot XML are not supported".into(),
                    ));
                }
                Event::DocType(_) => {
                    return Err(Error::Unsupported("XML DTDs are not supported".into()));
                }
                Event::Eof => {
                    if !root_opened || !root_closed || !self.tags_open.is_empty() {
                        return Err(parse("incomplete Mascot XML document"));
                    }
                    return Ok(());
                }
                _ => {}
            }
        }
    }

    fn append_text(&mut self, decoded: &str) -> Result<()> {
        if self.buffer.len() + decoded.len() > self.limits.max_text_bytes {
            return Err(bad("Mascot XML element text limit exceeded"));
        }
        self.buffer.push_str(decoded);
        Ok(())
    }

    fn start_element(&mut self, attributes: &BTreeMap<String, String>) -> Result<()> {
        match self.tag.as_str() {
            "mascot_search_results" => {
                self.major_version = attributes
                    .get("majorVersion")
                    .cloned()
                    .ok_or_else(|| parse("mascot_search_results requires majorVersion"))?;
                self.no_rt_error = false;
            }
            "protein" => {
                self.protein_hit.accession = attributes
                    .get("accession")
                    .cloned()
                    .ok_or_else(|| parse("protein requires accession"))?;
            }
            "query" => {
                let number = attributes
                    .get("number")
                    .ok_or_else(|| parse("query requires number"))?;
                self.query = u32::try_from(integer(number, "query number")?)
                    .map_err(|_| parse("query number must be positive"))?;
                if self.query == 0 {
                    return Err(parse("query numbers count from one"));
                }
            }
            "peptide" | "u_peptide" | "q_peptide" => {
                let number = attributes
                    .get("query")
                    .ok_or_else(|| parse("peptide requires query"))?;
                let value = integer(number, "peptide query")?;
                // The source computes `query - 1` into an unsigned member and
                // then compares it with `>` against the size, so query 0 wraps
                // to a huge index — which that guard does catch — and
                // query == size + 1 indexes one past the end, which it does
                // not. Both are refused here.
                let index = value
                    .checked_sub(1)
                    .and_then(|v| usize::try_from(v).ok())
                    .ok_or_else(|| parse("peptide query numbers count from one"))?;
                if index >= self.declared_queries {
                    return Err(parse(
                        "No or conflicting header information present (make sure to use the 'show_header=1' option in the ./export_dat.pl script)",
                    ));
                }
                self.index = index;
            }
            _ => {}
        }
        Ok(())
    }

    fn end_element(&mut self) -> Result<()> {
        let body = std::mem::take(&mut self.buffer);
        let text = trim(&body);
        match self.tag.as_str() {
            "NumQueries" => {
                let count = integer(text, "NumQueries")?;
                let count =
                    usize::try_from(count).map_err(|_| parse("NumQueries must not be negative"))?;
                if count > self.limits.max_queries {
                    return Err(bad("Mascot XML query limit exceeded"));
                }
                // The source resizes the identification vector here, which
                // commits the whole declared count — hundreds of megabytes for
                // a five-million-query header, and again for every repetition
                // of the element. The count is recorded as the index space
                // instead and entries are materialised when a query actually
                // references them; a *smaller* repeat still truncates, which is
                // the data loss `resize` causes.
                self.declared_queries = count;
                self.ids.truncate(count);
            }
            "prot_score" => {
                // The source converts a protein score with toInt32, so a
                // fractional score is a conversion error, not a rounded value.
                self.protein_hit.score = f64::from(integer(text, "prot_score")?);
            }
            "pep_exp_mz" => {
                let mz = number(text, "pep_exp_mz")?;
                self.identification_mut()?.mz = Some(mz);
            }
            "pep_scan_title" => self.scan_title(text)?,
            "pep_exp_z" => self.peptide_hit.charge = integer(text, "pep_exp_z")?,
            "pep_score" => self.peptide_hit.score = number(text, "pep_score")?,
            "pep_expect" => {
                let value = number(text, "pep_expect")?;
                self.peptide_hit
                    .metadata
                    .insert(EVALUE_KEY.to_owned(), MetaValue::try_from(value)?);
            }
            "pep_homol" => {
                let value = number(text, "pep_homol")?;
                self.identification_mut()?.significance_threshold = value;
            }
            "pep_ident" => {
                // Matrix Science: the homology threshold is used only when it
                // exists and is smaller than the identity threshold.
                let homology = self.identification_mut()?.significance_threshold;
                let identity = number(text, "pep_ident")?;
                self.peptide_hit.metadata.insert(
                    HOMOLOGY_THRESHOLD_KEY.to_owned(),
                    MetaValue::try_from(homology)?,
                );
                self.peptide_hit.metadata.insert(
                    IDENTITY_THRESHOLD_KEY.to_owned(),
                    MetaValue::try_from(identity)?,
                );
                if homology > identity || homology == 0.0 {
                    self.identification_mut()?.significance_threshold = identity;
                }
            }
            "pep_seq" => self.peptide_sequence(text)?,
            "pep_res_before" => {
                if let Some(residue) = text.chars().next() {
                    self.evidence.aa_before = flanking_residue(residue, true)?;
                }
            }
            "pep_res_after" => {
                if let Some(residue) = text.chars().next() {
                    self.evidence.aa_after = flanking_residue(residue, false)?;
                }
            }
            "pep_var_mod_pos" => self.variable_modification_positions(text)?,
            "Date" => {
                let parts = split(text, 'T');
                if parts.len() == 2 {
                    let time = match parts[1].find('Z') {
                        Some(position) => parts[1].get(..position).unwrap_or(parts[1]),
                        None => parts[1],
                    };
                    let value = format!("{} {time}", parts[0]);
                    value.parse::<CompletionTime>()?;
                    self.identifier = format!("Mascot_{value}");
                    self.protein.date_time = Some(value);
                }
            }
            "StringTitle" => self.string_title(text)?,
            "RTINSECONDS" => {
                let rt = number(text, "RTINSECONDS")?;
                self.query_identification_mut()?.rt = Some(rt);
            }
            "MascotVer" => self.protein.search_engine_version = text.to_owned(),
            "DB" => self.search.database = text.to_owned(),
            "FastaVer" => self.search.database_version = text.to_owned(),
            "TAXONOMY" => self.search.taxonomy = text.to_owned(),
            "CHARGE" => self.search.charges = text.to_owned(),
            "PFA" => {
                self.search.missed_cleavages = u32::try_from(integer(text, "PFA")?)
                    .map_err(|_| parse("PFA must not be negative"))?;
            }
            "MASS" => match text {
                "Monoisotopic" => self.search.mass_type = PeakMassType::Monoisotopic,
                "Average" => self.search.mass_type = PeakMassType::Average,
                // Any other spelling leaves the default untouched.
                _ => {}
            },
            "MODS" => {
                // Read here only when no <fixed_mods> section was present.
                if self.search.fixed_modifications.is_empty() {
                    let mut collected = Vec::new();
                    for entry in split(text, ',') {
                        if self.remove_fixed_mods.iter().any(|m| m == entry) {
                            continue;
                        }
                        collected.extend(self.split_modification(entry)?);
                    }
                    self.push_modifications(collected, true)?;
                }
            }
            "IT_MODS" => {
                // Read here only when no <variable_mods> section was present,
                // because Mascot sometimes forces a user-set fixed
                // modification to be variable.
                if self.search.variable_modifications.is_empty() {
                    let mut collected = Vec::new();
                    for entry in split(text, ',') {
                        collected.extend(self.split_modification(entry)?);
                    }
                    self.push_modifications(collected, false)?;
                }
            }
            "CLE" => {
                if self.proteases.has_enzyme(text) {
                    self.search.digestion_enzyme =
                        self.proteases.get_enzyme(text)?.name().to_owned();
                }
            }
            "TOL" => {
                self.search.precursor_tolerance = tolerance(
                    number(text, "TOL")?,
                    matches!(
                        self.search.precursor_tolerance,
                        crate::comparison::Tolerance::Ppm(_)
                    ),
                )
            }
            "ITOL" => {
                self.search.fragment_tolerance = tolerance(
                    number(text, "ITOL")?,
                    matches!(
                        self.search.fragment_tolerance,
                        crate::comparison::Tolerance::Ppm(_)
                    ),
                )
            }
            "TOLU" => {
                let value = current(self.search.precursor_tolerance);
                self.search.precursor_tolerance = tolerance(value, text == "ppm");
            }
            "ITOLU" => {
                let value = current(self.search.fragment_tolerance);
                self.search.fragment_tolerance = tolerance(value, text == "ppm");
            }
            "name" => self.modification_name(text)?,
            "warning" => self.mascot_warning(&body),
            "protein" => {
                self.protein.score_type = SCORE_TYPE.to_owned();
                if self.protein.hits.len() >= self.limits.max_hits {
                    return Err(bad("Mascot XML protein hit limit exceeded"));
                }
                self.protein
                    .hits
                    .push(std::mem::take(&mut self.protein_hit));
            }
            "peptide" => self.close_peptide()?,
            "u_peptide" | "q_peptide" => {
                let identifier = self.identifier.clone();
                let hit = std::mem::take(&mut self.peptide_hit);
                let limit = self.limits.max_hits;
                let identification = self.identification_mut()?;
                identification.identifier = identifier;
                identification.score_type = SCORE_TYPE.to_owned();
                if identification.hits.len() >= limit {
                    return Err(bad("Mascot XML peptide hit limit exceeded"));
                }
                identification.hits.push(hit);
                self.evidence = PeptideEvidence::default();
            }
            "mascot_search_results" => {
                self.protein.search_engine = SCORE_TYPE.to_owned();
                self.protein.identifier = self.identifier.clone();
                // Split variable modifications now that the per-index order is
                // no longer needed, e.g. Phospho (ST) -> Phospho (S), Phospho (T).
                let mut split_mods = Vec::new();
                for modification in std::mem::take(&mut self.search.variable_modifications) {
                    split_mods.extend(self.split_modification(&modification)?);
                }
                self.push_modifications(split_mods, false)?;
                self.protein.search_parameters = self.search.clone();
            }
            _ => {}
        }
        Ok(())
    }

    /// The identification at `index`, materialised if `<NumQueries>` declared
    /// it but no earlier query reached it.
    fn identification_at(
        &mut self,
        index: usize,
        message: &'static str,
    ) -> Result<&mut PeptideIdentification> {
        if index >= self.declared_queries {
            return Err(parse(message));
        }
        if index >= self.ids.len() {
            let additional = index + 1 - self.ids.len();
            self.ids
                .try_reserve(additional)
                .map_err(|_| bad("identification allocation failed"))?;
            self.ids.resize(index + 1, PeptideIdentification::default());
        }
        self.ids
            .get_mut(index)
            .ok_or_else(|| bad("identification allocation failed"))
    }

    fn identification_mut(&mut self) -> Result<&mut PeptideIdentification> {
        let index = self.index;
        self.identification_at(
            index,
            "No or conflicting header information present (make sure to use the 'show_header=1' option in the ./export_dat.pl script)",
        )
    }

    /// The identification the current `<query number=...>` refers to.
    ///
    /// The source indexes `id_data_[actual_query_ - 1]` with no check at all,
    /// and `actual_query_` is unsigned, so a `<query number="0">` reads at
    /// index 4294967295 and a number beyond `<NumQueries>` past the end. Both
    /// are refused here.
    fn query_identification_mut(&mut self) -> Result<&mut PeptideIdentification> {
        let index = usize::try_from(self.query)
            .ok()
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| parse("query numbers count from one"))?;
        self.identification_at(index, "query number exceeds NumQueries")
    }

    fn scan_title(&mut self, title: &str) -> Result<()> {
        let want_mz = self.identification_mut()?.mz.is_none();
        match self.lookup.spectrum_meta_data(title, want_mz) {
            Ok(meta) => {
                let identification = self.identification_mut()?;
                identification.rt = meta.rt;
                if want_mz {
                    identification.mz = meta.precursor_mz;
                }
            }
            Err(error) => {
                self.warnings.push(format!(
                    "<pep_scan_title> element has unexpected format '{title}'. Could not extract spectrum meta data. ({error})"
                ));
            }
        }
        // The source guard is `if (!id_data_[i].getRT())`, which is false for
        // the NaN it just assigned, so an unresolved title reports nothing and
        // a title that legitimately encodes RT 0 reports an error. This port
        // reports the unresolved case, which is the one a caller can act on.
        let resolved = self.identification_mut()?.rt.is_some();
        if !resolved && !self.no_rt_error {
            let reference = if self.lookup.is_empty() {
                String::new()
            } else {
                "or a matching spectrum reference ".to_owned()
            };
            self.warnings.push(format!(
                "Could not extract RT value {reference}from <pep_scan_title> element with format '{title}'. Try adjusting the 'scan_regex' parameter."
            ));
            self.no_rt_error = true;
        }
        Ok(())
    }

    fn peptide_sequence(&mut self, text: &str) -> Result<()> {
        let mut sequence = AASequence::parse(text)?;
        if self.peptides.is_empty() {
            for modification in self.search.fixed_modifications.clone() {
                self.apply_fixed_modification(&mut sequence, &modification)?;
            }
        }
        self.peptide_hit.sequence = sequence;
        Ok(())
    }

    /// One entry of `search_parameters.fixed_modifications` applied to a hit.
    ///
    /// The name and its specificity are separated by spaces, so
    /// `Carboxymethyl (C)` modifies every C, `Acetyl (Protein N-term)` and
    /// `Acetyl (N-term)` the N terminus, and `Acetyl (N-term C)` the N terminus
    /// only when the first residue is C.
    fn apply_fixed_modification(
        &mut self,
        sequence: &mut AASequence,
        modification: &str,
    ) -> Result<()> {
        let parts = split(modification, ' ');
        if parts.len() < 2 || parts.len() > 3 {
            self.warnings
                .push(format!("Cannot parse fixed modification '{modification}'"));
            return Ok(());
        }
        let name = parts[0];
        let terminal_protein = |keyword: &str| {
            parts[1] == "(Protein" && parts.len() == 3 && parts[2] == format!("{keyword})")
        };
        if parts[1] == "(C-term)" || terminal_protein("C-term") {
            sequence.set_c_terminal_modification_with_registry(name, self.modifications)?;
        } else if parts[1] == "(N-term)" || terminal_protein("N-term") {
            sequence.set_n_terminal_modification_with_registry(name, self.modifications)?;
        } else if parts[1] == "(C-term" && parts.len() == 3 {
            // The source dereferences end() - 1 without checking, which is
            // undefined behaviour on an empty sequence.
            let residue = parts[2].replace(')', "");
            if last_residue(sequence).is_some_and(|code| code == residue) {
                sequence.set_c_terminal_modification_with_registry(name, self.modifications)?;
            }
        } else if parts[1] == "(N-term" && parts.len() == 3 {
            let residue = parts[2].replace(')', "");
            if first_residue(sequence).is_some_and(|code| code == residue) {
                sequence.set_n_terminal_modification_with_registry(name, self.modifications)?;
            }
        } else {
            let residue = parts[1].replace([')', '('], "");
            let codes: Vec<String> = sequence
                .as_str()
                .chars()
                .map(|c| c.to_string())
                .collect::<Vec<_>>();
            for (position, code) in codes.iter().enumerate() {
                if *code == residue {
                    sequence.set_modification_with_registry(position, name, self.modifications)?;
                }
            }
        }
        Ok(())
    }

    /// `<pep_var_mod_pos>`, three `.`-separated fields: the N-terminal slot,
    /// one digit per residue and the C-terminal slot. A non-`0` digit is the
    /// one-based index into `variable_modifications`, which is why that list
    /// must not be expanded before the document ends.
    fn variable_modification_positions(&mut self, text: &str) -> Result<()> {
        let parts = split(text, '.');
        if parts.len() != 3 {
            return Ok(());
        }
        let mut sequence = self.peptide_hit.sequence.clone();
        let length = sequence.len();
        for (position, digit) in parts[1].chars().enumerate() {
            if digit == '0' {
                continue;
            }
            let name = self.variable_modification(digit)?;
            if position >= length {
                return Err(parse(
                    "<pep_var_mod_pos> marks a modification beyond the peptide sequence",
                ));
            }
            match name {
                Some(name) => {
                    sequence.set_modification_with_registry(position, &name, self.modifications)?;
                }
                None => self
                    .warnings
                    .push("Cannot parse variable modification".to_owned()),
            }
        }
        if let Some(digit) = parts[0].chars().next().filter(|&c| c != '0') {
            match self.variable_modification(digit)? {
                Some(name) => {
                    sequence.set_n_terminal_modification_with_registry(&name, self.modifications)?
                }
                None => self
                    .warnings
                    .push("Cannot parse variable N-term modification".to_owned()),
            }
        }
        if let Some(digit) = parts[2].chars().next().filter(|&c| c != '0') {
            match self.variable_modification(digit)? {
                Some(name) => {
                    sequence.set_c_terminal_modification_with_registry(&name, self.modifications)?
                }
                None => self
                    .warnings
                    .push("Cannot parse variable C-term modification".to_owned()),
            }
        }
        self.peptide_hit.sequence = sequence;
        Ok(())
    }

    /// The bare modification name for a `<pep_var_mod_pos>` digit, or `None`
    /// when the configured entry has no specificity part to strip.
    fn variable_modification(&self, digit: char) -> Result<Option<String>> {
        let position = digit
            .to_digit(10)
            .and_then(|value| value.checked_sub(1))
            .and_then(|value| usize::try_from(value).ok())
            .ok_or_else(|| parse("invalid <pep_var_mod_pos> digit"))?;
        // The source uses vector::at, which throws when the index is too large.
        let entry = self
            .search
            .variable_modifications
            .get(position)
            .ok_or_else(|| {
                parse("<pep_var_mod_pos> index exceeds the declared variable modifications")
            })?;
        let parts = split(entry, ' ');
        Ok(if parts.len() >= 2 {
            Some(parts[0].to_owned())
        } else {
            None
        })
    }

    /// `<StringTitle>`: replace hit sequences from the modified-peptide map,
    /// then, if the query still has no retention time, read one from a
    /// `<m/z>_<RT>` title.
    fn string_title(&mut self, title: &str) -> Result<()> {
        if let Some(supplied) = self.peptides.get(title) {
            let identification = self.query_identification_mut()?;
            let mut hits = std::mem::take(&mut identification.hits);
            if supplied.len() != hits.len() {
                self.warnings
                    .push("pepXML hits and Mascot hits are not the same".to_owned());
            }
            // pepXML can hold more hits than Mascot XML, so every Mascot hit is
            // matched against every supplied sequence. Quadratic in the hit
            // count, which is always small.
            for hit in &mut hits {
                for candidate in supplied {
                    if candidate.is_modified() && candidate.as_str() == hit.sequence.as_str() {
                        hit.sequence = candidate.clone();
                        break;
                    }
                }
            }
            self.query_identification_mut()?.hits = hits;
        }
        if self.query_identification_mut()?.rt.is_none() {
            let parts = split(title, '_');
            if parts.len() == 2 {
                let rt = number(parts[1], "StringTitle retention time")?;
                self.query_identification_mut()?.rt = Some(rt);
            }
        }
        Ok(())
    }

    /// A `<name>` element inside `<fixed_mods>` or `<variable_mods>`.
    ///
    /// Mascot XML 1.x has `<name>` only inside `<variable_mods>`; from 2.1 both
    /// sections have one, so the parent element decides.
    fn modification_name(&mut self, text: &str) -> Result<()> {
        // The closing `</name>` has already been popped, so the enclosing
        // section is the second entry from the top, as `tags_open_.size() - 2`
        // is in the source.
        let parent = self
            .tags_open
            .len()
            .checked_sub(2)
            .and_then(|index| self.tags_open.get(index))
            .map(String::as_str)
            .unwrap_or("");
        if self.major_version == "1" || parent == "variable_mods" {
            // Phospho (ST) cannot be split here: the order of variable
            // modifications is the index space <pep_var_mod_pos> refers to.
            self.push_modifications(vec![text.to_owned()], false)?;
        } else if parent == "fixed_mods" {
            if self.remove_fixed_mods.iter().any(|m| m == text) {
                self.warnings.push(format!(
                    "Modification removed as fixed modification: '{text}'"
                ));
            } else {
                let split_mods = self.split_modification(text)?;
                self.push_modifications(split_mods, true)?;
            }
        }
        Ok(())
    }

    /// A `<warning>` element. Mascot reports a fixed modification it can only
    /// apply as a variable one; that entry must not enter the fixed list.
    fn mascot_warning(&mut self, body: &str) {
        self.warnings
            .push(format!("Warnings were present: '{body}'"));
        let trimmed = trim(body);
        if !trimmed.contains("can only be used as a variable modification") {
            return;
        }
        let parts = split(trimmed, ';');
        let Some(first) = parts.first() else {
            return;
        };
        let Some(rest) = first.strip_prefix('\'') else {
            return;
        };
        if let Some(end) = rest.find('\'') {
            if let Some(name) = rest.get(..end) {
                self.remove_fixed_mods.push(name.to_owned());
            }
        }
    }

    /// `</peptide>`: a hit whose sequence is already present only contributes
    /// another protein accession as peptide evidence.
    fn close_peptide(&mut self) -> Result<()> {
        let accession = self.protein_hit.accession.clone();
        let identifier = self.identifier.clone();
        let sequence = self.peptide_hit.sequence.clone();
        let limit = self.limits.max_hits;
        let mut evidence = std::mem::take(&mut self.evidence);
        evidence.protein_accession = accession;
        let hit = std::mem::take(&mut self.peptide_hit);
        let identification = self.identification_mut()?;
        match identification
            .hits
            .iter_mut()
            .find(|stored| stored.sequence == sequence)
        {
            Some(stored) => {
                if stored.evidences.len() >= limit {
                    return Err(bad("Mascot XML peptide evidence limit exceeded"));
                }
                stored.evidences.push(evidence);
            }
            None => {
                identification.identifier = identifier;
                identification.score_type = SCORE_TYPE.to_owned();
                if identification.hits.len() >= limit {
                    return Err(bad("Mascot XML peptide hit limit exceeded"));
                }
                let mut hit = hit;
                hit.evidences.push(evidence);
                identification.hits.push(hit);
            }
        }
        Ok(())
    }

    /// `splitModificationBySpecifiedAA`: expand `Phospho (ST)` into
    /// `Phospho (S)` and `Phospho (T)`, leaving terminal specifications and
    /// anything that is not exactly "name (residues)" untouched.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when an expanded
    /// identifier is not in the modification database, which is where the
    /// source throws `Exception::ElementNotFound`.
    fn split_modification(&self, modification: &str) -> Result<Vec<String>> {
        let parts = split(modification, ' ');
        if parts.len() != 2 {
            return Ok(vec![modification.to_owned()]);
        }
        if parts[1].starts_with("(N-term") || parts[1].starts_with("(C-term") {
            return Ok(vec![modification.to_owned()]);
        }
        let residues = parts[1].replace([')', '('], "");
        let mut expanded = Vec::new();
        for residue in residues.chars() {
            let candidate = format!("{} ({residue})", parts[0]);
            // `ModificationsDB::has` matches any registered identifier, which
            // is what `find` with no residue or terminus restriction does.
            if self.modifications.find(&candidate, None, None).is_empty() {
                return Err(bad(format!("modification '{candidate}' not found")));
            }
            expanded.push(candidate);
        }
        Ok(expanded)
    }

    fn push_modifications(&mut self, values: Vec<String>, fixed: bool) -> Result<()> {
        let target = if fixed {
            &mut self.search.fixed_modifications
        } else {
            &mut self.search.variable_modifications
        };
        if target.len() + values.len() > self.limits.max_modifications {
            return Err(bad("Mascot XML modification limit exceeded"));
        }
        target
            .try_reserve(values.len())
            .map_err(|_| bad("modification allocation failed"))?;
        target.extend(values);
        Ok(())
    }

    /// `MascotXMLFile::load`'s post-processing of the handler's output.
    fn finish(&mut self) -> Result<MascotXmlResult> {
        // The source reserves the unfiltered count, so a document declaring
        // five million queries and holding none still commits a five-million
        // entry vector — and the empty result keeps that capacity. Counting
        // the survivors first costs one pass and reserves what is needed.
        let keeps = |identification: &PeptideIdentification| {
            let hits = &identification.hits;
            !hits.is_empty() && (hits.len() > 1 || !hits[0].sequence.is_empty())
        };
        let mut kept: Vec<PeptideIdentification> = Vec::new();
        kept.try_reserve_exact(self.ids.iter().filter(|id| keeps(id)).count())
            .map_err(|_| bad("identification allocation failed"))?;
        let mut missing_sequence = 0usize;
        for identification in std::mem::take(&mut self.ids) {
            if keeps(&identification) {
                kept.push(identification);
            } else if !identification.hits.is_empty() {
                missing_sequence += 1;
            }
        }
        if missing_sequence > 0 {
            self.warnings.push(format!(
                "Warning: Removed {missing_sequence} peptide identifications without sequence."
            ));
        }
        let no_rt = kept.iter().filter(|id| id.rt.is_none()).count();
        if no_rt > 0 {
            self.warnings.push(format!(
                "Warning: {no_rt} (of {}) peptide identifications have no retention time value.",
                kept.len()
            ));
        }
        // A supplied mapping that resolved nothing is an error, not a warning.
        // The source compares the two counts without excluding the empty case,
        // so a non-empty lookup over a document with no surviving
        // identifications is refused as well; that is reproduced here.
        if !self.lookup.is_empty() && no_rt == kept.len() {
            return Err(Error::MissingInformation(
                "No retention time information for peptide identifications found".into(),
            ));
        }
        // Mascot 2.2 repeats the first hit, so one of the two is dropped. The
        // source comment says "erase first hit" but the code erases the second,
        // which is the behaviour reproduced here.
        for identification in &mut kept {
            let duplicate = identification.hits.len() > 1
                && identification.hits[0].score == identification.hits[1].score
                && identification.hits[0].sequence == identification.hits[1].sequence
                && identification.hits[0].charge == identification.hits[1].charge;
            if duplicate {
                identification.hits.remove(1);
            }
        }
        Ok(MascotXmlResult {
            protein_identification: std::mem::take(&mut self.protein),
            peptide_identifications: kept,
            warnings: std::mem::take(&mut self.warnings),
        })
    }
}

/// One `<pep_res_before>` or `<pep_res_after>` character.
///
/// Mascot writes `-` for a protein terminus. The source stores that character
/// verbatim on the peptide evidence, and `IdXMLFile` then writes
/// `aa_before="-"`, which is neither of OpenMS's own `[`/`]` terminus markers
/// and which nothing downstream interprets. This port maps `-` to the marker it
/// means — `[` before the peptide, `]` after it — because `FlankingResidue`
/// has no verbatim variant and the terminus is the information the character
/// carries. `docs/MASCOT_XML_SUPPORT.md` records the divergence.
fn flanking_residue(code: char, before: bool) -> Result<FlankingResidue> {
    if code == '-' {
        return Ok(if before {
            FlankingResidue::NTerminus
        } else {
            FlankingResidue::CTerminus
        });
    }
    FlankingResidue::from_code(code)
}

fn tolerance(value: f64, ppm: bool) -> crate::comparison::Tolerance {
    if ppm {
        crate::comparison::Tolerance::Ppm(value)
    } else {
        crate::comparison::Tolerance::Absolute(value)
    }
}
fn current(value: crate::comparison::Tolerance) -> f64 {
    let (crate::comparison::Tolerance::Absolute(v) | crate::comparison::Tolerance::Ppm(v)) = value;
    v
}
fn first_residue(sequence: &AASequence) -> Option<String> {
    sequence.as_str().chars().next().map(|c| c.to_string())
}
fn last_residue(sequence: &AASequence) -> Option<String> {
    sequence.as_str().chars().next_back().map(|c| c.to_string())
}

/// Expand one XML reference that quick-xml emitted as its own event.
///
/// quick-xml does not expand references inside text: it splits the text at
/// every `&...;` and hands the reference back separately. The five predefined
/// entities and numeric character references are resolvable without a DTD and
/// are expanded here, because they are ordinary in Mascot output. Anything else
/// would need a DTD, which this reader refuses outright, so it is an error
/// rather than a silent deletion — dropping it concatenated the surrounding
/// fragments and turned `1&#46;5` into the score `15`.
fn expand_reference(raw: &[u8]) -> Result<String> {
    let name = std::str::from_utf8(raw)
        .map_err(|_| parse("non-UTF-8 XML entity reference in Mascot XML"))?;
    match name {
        "amp" => return Ok("&".into()),
        "lt" => return Ok("<".into()),
        "gt" => return Ok(">".into()),
        "quot" => return Ok("\"".into()),
        "apos" => return Ok("'".into()),
        _ => {}
    }
    let digits = name.strip_prefix('#').ok_or_else(|| {
        Error::Unsupported(format!(
            "XML entity reference `&{name};` needs a DTD, which is not supported"
        ))
    })?;
    let code = match digits
        .strip_prefix('x')
        .or_else(|| digits.strip_prefix('X'))
    {
        Some(hex) => u32::from_str_radix(hex, 16),
        None => digits.parse::<u32>(),
    }
    .map_err(|_| parse("malformed numeric character reference in Mascot XML"))?;
    let character = char::from_u32(code)
        .ok_or_else(|| parse("numeric character reference is not a Unicode scalar value"))?;
    Ok(character.to_string())
}

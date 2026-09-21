// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Mascot generic format (MGF) reader, writer and search-submission header,
//! ported from `FORMAT/MascotGenericFile.h` and its implementation.
//!
//! A file is a sequence of `BEGIN IONS` … `END IONS` blocks. Each block carries
//! `KEY=value` header lines followed by whitespace-separated m/z–intensity
//! pairs. For the format itself see the Matrix Science description of the
//! generic file format (`data_file_help.html#GEN`).
//!
//! The source reader is deliberately lenient in some places and strict in
//! others, and this module reproduces both; every case is listed in
//! `docs/MASCOT_GENERIC_SUPPORT.md`.
//!
//! Entry points are the free functions
//! [`read`](crate::format::mascot_generic::read),
//! [`load`](crate::format::mascot_generic::load) and
//! [`consume`](crate::format::mascot_generic::consume) for the streaming
//! reader, and
//! [`MascotGenericFile`](crate::format::mascot_generic::MascotGenericFile) for
//! the writer, which owns the Mascot search parameters that become the MGF
//! parameter header.
//!
//! The source class derives from `ProgressLogger`; [`load_with_progress`](crate::format::mascot_generic::load_with_progress),
//! [`MascotGenericFile::store_with_progress`](crate::format::mascot_generic::MascotGenericFile::store_with_progress)
//! and [`MascotGenericFile::store_to_with_progress`](crate::format::mascot_generic::MascotGenericFile::store_to_with_progress)
//! make the progress calls of its `load` and `store` on a caller's logger;
//! every other entry point reports nothing.
//!
//! This module is independent of [`crate::format::mgf`], a stricter native MGF
//! interchange adapter that rejects much of the malformed input the Mascot
//! reader tolerates. Use this module when source fidelity matters.

pub use super::ms2::Limits;
use super::ms2::{
    Counter, TextInput, increment, invalid, is_space, push_peak, push_spectrum, trim, unsupported,
};
use super::parse_error;
use crate::chemistry::ModificationsDB;
use crate::concept::constants::user_param::{
    MSM_INCHI_STRING, MSM_METABOLITE_NAME, MSM_SMILES_STRING,
};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::interfaces::MSDataConsumer;
use crate::kernel::SpectrumType;
use crate::metadata::MetaValue;
use crate::param::value::{format_float, format_float32};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{MSExperiment, MSSpectrum, Peak1D, Precursor, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::{ControlFlow, Range};
use std::path::{Path, PathBuf};

/// Upper bound on peaks per written spectrum, from `MascotGenericFile.cpp`'s
/// `spec.size() >= 10000` guard. The source message explains the limit as a
/// profile-data guard: MGF is centroided by definition.
pub const MAX_WRITTEN_PEAKS: usize = 10_000;
/// Bound on `SEQ=` lines accumulated for one query.
pub const MAX_SEQ_ENTRIES: usize = 100_000;
/// Bound on identifiers in one modification parameter list.
pub const MAX_MODIFICATIONS: usize = 100_000;

/// Meta key under which a parsed `TITLE=` line is stored, as in the source.
pub const TITLE_KEY: &str = "TITLE";
/// Meta key under which `SEQ=` lines accumulate, always as a string list.
pub const SEQ_KEY: &str = "SEQ";
/// Meta key for a `SCANS=` line. The source deliberately renames it.
pub const SCAN_ID_KEY: &str = "Scan_ID";
/// Meta key for a `SPECTRUMID=` line, the GNPS library accession.
pub const GNPS_SPECTRUM_ID_KEY: &str = "GNPS_Spectrum_ID";
/// Sentinel the source substitutes for a missing or empty native-ID accession.
pub const UNKNOWN_NATIVE_ID_TYPE: &str = "UNKNOWN";

/// Per-spectrum state the source reader carries from one block to the next.
///
/// `MascotGenericFile::load` declares one spectrum outside its read loop, and
/// `getNextSpectrum_` clears only the peak list, the native ID, `TITLE` and
/// `SEQ`. Retention time, precursor m/z, precursor intensity, precursor
/// charge, MS level and every other meta value therefore survive into the next
/// block, so a block that omits `CHARGE=` silently inherits the previous
/// block's charge. The upstream class test asserts the non-inheritance of
/// `SEQ` only, which is the one field the source explicitly resets.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CarryOver {
    /// Start every block from a fresh record. This is the native default,
    /// because an inherited precursor or retention time is not recoverable by
    /// a caller and is almost never what the file meant.
    #[default]
    Reset,
    /// Reproduce the source: clear only peaks, native ID, `TITLE` and `SEQ`.
    Source,
}

/// Reader configuration: three range filters, the carry-over policy and the
/// shared text-adapter resource ceilings.
///
/// Ranges are half-open, including the minimum and excluding the maximum, as in
/// OpenMS `DRange`. `MascotGenericFile` is not a `PeakFileOptions` consumer —
/// the source reader applies no filtering at all — so these three filters are a
/// native addition mirroring [`crate::format::dta2d`], where the source does
/// consume them. [`MascotGenericReader`] is a streaming iterator, so a filtered
/// peak is never retained.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadOptions {
    /// Spectra whose retention time falls outside this range are dropped.
    pub rt_range: Option<Range<f64>>,
    /// Peaks whose m/z falls outside this range are dropped.
    pub mz_range: Option<Range<f64>>,
    /// Peaks whose intensity falls outside this range are dropped.
    pub intensity_range: Option<Range<f64>>,
    /// Whether unset fields inherit the previous block's values.
    pub carry_over: CarryOver,
    /// Byte, line, spectrum and peak ceilings, shared with the other text
    /// adapters. Counts include filtered and discarded input.
    pub limits: Limits,
}
impl ReadOptions {
    fn validate(&self) -> Result<()> {
        self.limits.validate()?;
        for range in [&self.rt_range, &self.mz_range, &self.intensity_range]
            .into_iter()
            .flatten()
        {
            if !range.start.is_finite() || !range.end.is_finite() || range.start > range.end {
                return Err(invalid("MGF ranges require finite ordered endpoints"));
            }
        }
        Ok(())
    }
}
fn includes(range: &Option<Range<f64>>, value: f64) -> bool {
    range.as_ref().is_none_or(|r| r.contains(&value))
}

/// Writer resource ceilings. The complete output is measured against
/// `max_output_bytes` before any byte is written or any file is created.
pub type WriteOptions = Limits;

/// `StringUtils::simplify`: every maximal run of space, tab, CR or LF collapses
/// to a single space. The source then substitutes tab for space redundantly.
fn simplify(text: &str) -> String {
    let mut result = String::with_capacity(text.len());
    let mut last_ws = false;
    for c in text.chars() {
        if is_space(c) {
            if !last_ws {
                result.push(' ');
            }
            last_ws = true;
        } else {
            result.push(c);
            last_ws = false;
        }
    }
    result
}

/// `StringUtils::split(s, c, out)` semantics: an empty input yields no chunks,
/// an input without the separator yields the whole string as one chunk, and
/// otherwise the separator partitions the input, keeping empty chunks.
fn source_split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// `StringUtils::toDouble`: leading and trailing space, tab, CR and LF are
/// skipped, one leading `+` is consumed, and the remainder must be a complete
/// decimal, `inf` or `nan` token. The source's `nan(payload)` extension is not
/// accepted, because Rust's parser has no such form and no caller writes one.
///
/// The remainder then goes to `std::from_chars`, which refuses a *second* `+`
/// and reports an overflowing decimal literal as `result_out_of_range` — a
/// conversion error, not an infinity. Rust's own parser accepts both, so `++5`
/// and `1e999` are rejected here explicitly; only the `inf`/`nan` words
/// themselves convert to a non-finite value, exactly as they do upstream.
fn source_double(text: &str) -> Option<f64> {
    let text = trim(text);
    let body = text.strip_prefix('+').unwrap_or(text);
    if body.starts_with('+') {
        return None;
    }
    let value = body.parse::<f64>().ok()?;
    if value.is_finite() {
        return Some(value);
    }
    let word = body.strip_prefix('-').unwrap_or(body);
    let head = word.get(..3).unwrap_or("");
    (head.eq_ignore_ascii_case("inf") || head.eq_ignore_ascii_case("nan")).then_some(value)
}

/// [`source_double`] restricted to finite values.
///
/// The source stores an infinity or a NaN silently; every later consumer of an
/// [`MSSpectrum`] rejects it, so it is refused at the point of parsing instead.
fn finite_double(text: &str, line: usize, label: &str) -> Result<f64> {
    source_double(text)
        .filter(|v| v.is_finite())
        .ok_or_else(|| parse_error(line, format!("invalid finite {label}: {text:?}")))
}

/// `StringUtils::toInt32`: leading and trailing space, tab, CR and LF, one
/// optional `+`, then a complete `i32`. Trailing characters and range overflow
/// are conversion errors, and so is a second `+`, which `std::from_chars`
/// refuses even though Rust's own parser accepts it.
fn source_int32(text: &str) -> Option<i32> {
    let text = trim(text);
    let body = text.strip_prefix('+').unwrap_or(text);
    if body.starts_with('+') {
        return None;
    }
    body.parse::<i32>().ok()
}

/// `std::stoi`, used only by the source's `MSLEVEL=` branch: leading
/// whitespace, an optional sign and the longest leading digit run, ignoring
/// whatever follows. `None` stands for the source's `std::invalid_argument`
/// and `Some(Err(()))` for its `std::out_of_range`.
fn source_stoi(text: &str) -> Option<std::result::Result<i32, ()>> {
    let bytes = text.as_bytes();
    let mut index = 0;
    while index < bytes.len() && matches!(bytes[index], b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
    {
        index += 1;
    }
    let start = index;
    if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
        index += 1;
    }
    let digits = index;
    while index < bytes.len() && bytes[index].is_ascii_digit() {
        index += 1;
    }
    if index == digits {
        return None;
    }
    // Every byte in start..index is ASCII, so the range is a character boundary.
    Some(text.get(start..index)?.parse::<i32>().map_err(|_| ()))
}

/// The suffix of `text` starting at byte `offset`, as `StringUtils::substr`
/// would return it.
///
/// `StringUtils::substr` clamps the start position to the string length, so a
/// header line shorter than its own key — a bare `NAME` or `MSLEVEL` — yields
/// an empty value rather than an error; that clamp is reproduced here.
///
/// Within the string the source slices raw bytes, so an offset landing inside a
/// multi-byte character yields an ill-formed string. This port refuses instead:
/// Rust string slicing at a non-boundary aborts the process, and half a
/// character is not a value any caller can use.
fn value_after<'a>(text: &'a str, offset: usize, line: usize, key: &str) -> Result<&'a str> {
    if offset >= text.len() {
        return Ok("");
    }
    text.get(offset..).ok_or_else(|| {
        parse_error(
            line,
            format!("MGF {key} value starts inside a multi-byte character"),
        )
    })
}

fn set_meta(
    spectrum: &mut MSSpectrum,
    key: &str,
    text: &str,
    offset: usize,
    line: usize,
    name: &str,
) -> Result<()> {
    let value = value_after(text, offset, line, name)?;
    spectrum
        .metadata
        .insert(key.to_owned(), MetaValue::from(value));
    Ok(())
}

/// The per-block starting point: MS level 2, centroided and one precursor,
/// exactly as `MascotGenericFile::load` initialises its reused spectrum. MGF is
/// always centroided by definition, which is why the type is not inferred.
fn base_spectrum() -> MSSpectrum {
    MSSpectrum {
        ms_level: 2,
        spectrum_type: SpectrumType::Centroid,
        precursors: vec![Precursor::default()],
        ..MSSpectrum::default()
    }
}

/// The single precursor the source indexes as `getPrecursors()[0]`.
fn precursor_mut(spectrum: &mut MSSpectrum) -> &mut Precursor {
    if spectrum.precursors.is_empty() {
        spectrum.precursors.push(Precursor::default());
    }
    &mut spectrum.precursors[0]
}

/// The `TITLE=` branch, which has two quite different halves.
///
/// When the line contains `min` anywhere, the source treats it as a Bruker
/// export like `TITLE= Cmpd 1, +MSn(595.3), 10.9 min`, splits on `,` and sets
/// the retention time from the first whitespace token of every chunk containing
/// `min`, read as minutes — so the last such chunk wins and no `TITLE` meta
/// value is stored at all. If one of those conversions fails, the source falls
/// back to storing the text between the first and second `=` as `TITLE`, which
/// truncates a title containing a second `=`.
///
/// Otherwise the value is the text after the first `=` at or after index 4,
/// suffixed with `_<native ID>` to keep titles unique — unless the already
/// stored `TITLE` contains the native ID, in which case a second `TITLE=` line
/// in the same block replaces it without the suffix.
///
/// The source calls `spectrum.setRT` *inside* the chunk loop and wraps the
/// whole loop in one `try`, so every conversion that succeeded before a later
/// one failed stays applied: `TITLE=run, 2 min, bad min, 3 min` keeps the
/// retention time 120 s from the second chunk *and* stores the fallback title.
fn read_title(spectrum: &mut MSSpectrum, text: &str, line: usize) -> Result<()> {
    if text.contains("min") {
        let mut failed = false;
        for chunk in source_split(text, ',') {
            if !chunk.contains("min") {
                continue;
            }
            let first = source_split(trim(chunk), ' ')
                .first()
                .copied()
                .unwrap_or("");
            match source_double(trim(first)) {
                Some(minutes) => {
                    let value = minutes * 60.0;
                    if !value.is_finite() {
                        return Err(parse_error(line, "MGF title retention time overflows"));
                    }
                    spectrum.rt = value;
                }
                None => {
                    failed = true;
                    break;
                }
            }
        }
        if failed {
            let parts = source_split(text, '=');
            if parts.len() >= 2 && !parts[1].is_empty() {
                spectrum
                    .metadata
                    .insert(TITLE_KEY.to_owned(), MetaValue::from(parts[1]));
            }
        }
        return Ok(());
    }
    let Some(equals) = text
        .as_bytes()
        .iter()
        .skip(4)
        .position(|&b| b == b'=')
        .map(|index| index + 4)
    else {
        return Ok(());
    };
    let value = value_after(text, equals + 1, line, "TITLE")?;
    let existing = spectrum
        .metadata
        .get(TITLE_KEY)
        .and_then(|v| v.as_str().ok())
        .unwrap_or("");
    let title = if existing.contains(&spectrum.native_id) {
        value.to_owned()
    } else {
        format!("{value}_{}", spectrum.native_id)
    };
    spectrum
        .metadata
        .insert(TITLE_KEY.to_owned(), MetaValue::from(title));
    Ok(())
}

/// Streaming Mascot-generic reader. One block is decoded per [`Iterator::next`].
///
/// Source `MascotGenericFile::getNextSpectrum_` is reproduced block for block,
/// including every case documented on [`read_with_options`].
pub struct MascotGenericReader<R> {
    input: TextInput<R>,
    options: ReadOptions,
    index: u32,
    spectra: usize,
    peaks: usize,
    template: MSSpectrum,
    finished: bool,
    line: String,
    /// Bytes of input read so far, line terminators included.
    consumed: usize,
}

impl<R: BufRead> MascotGenericReader<R> {
    /// Wrap a buffered stream without reading from it yet.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when the options
    /// exceed the hard ceilings or carry a non-finite or inverted range.
    pub fn new(reader: R, options: &ReadOptions) -> Result<Self> {
        options.validate()?;
        Ok(Self {
            input: TextInput::new(reader, &options.limits)?,
            options: options.clone(),
            index: 0,
            spectra: 0,
            peaks: 0,
            template: base_spectrum(),
            finished: false,
            line: String::new(),
            consumed: 0,
        })
    }

    /// The next input line into `self.line`, counting its bytes.
    fn next_line(&mut self) -> Result<bool> {
        let read = self.input.next_line(&mut self.line)?;
        if read {
            self.consumed = self.consumed.saturating_add(self.line.len());
        }
        Ok(read)
    }

    /// What the source's `is.tellg()` reports after the last line was read
    /// with `std::getline` (`MascotGenericFile.h:99`): the bytes consumed, or
    /// -1 once that line ended at the end of input without a newline, because
    /// `getline` then set `eofbit` and `tellg` fails.
    fn source_position(&self) -> Result<i64> {
        if !self.line.ends_with('\n') {
            return Ok(-1);
        }
        progress_value(self.consumed)
    }

    fn next_block(&mut self) -> Result<Option<MSSpectrum>> {
        let max_spectra = self.options.limits.max_spectra;
        let source_carry = self.options.carry_over == CarryOver::Source;
        let mut spectrum = if source_carry {
            self.template.clone()
        } else {
            base_spectrum()
        };
        spectrum.peaks.clear();
        spectrum.native_id = format!("index={}", self.index);
        spectrum.metadata.remove(TITLE_KEY);
        spectrum.metadata.remove(SEQ_KEY);
        // `SEQ=` lines accumulate here rather than through the meta value, so a
        // query with many of them costs one push each instead of copying the
        // whole list back and forth on every line as the source does.
        let mut sequences: Vec<String> = Vec::new();
        loop {
            if !self.next_line()? {
                return Ok(None);
            }
            if trim(&self.line) != "BEGIN IONS" {
                continue;
            }
            match self.read_block(&mut spectrum, &mut sequences)? {
                Some(()) => {
                    if !sequences.is_empty() {
                        spectrum
                            .metadata
                            .insert(SEQ_KEY.to_owned(), MetaValue::from(sequences));
                    }
                    increment(
                        &mut self.spectra,
                        max_spectra,
                        "MGF spectrum limit exceeded",
                    )?;
                    self.index = self.index.saturating_add(1);
                    if source_carry {
                        self.template = spectrum.clone();
                    }
                    return Ok(Some(spectrum));
                }
                // Source: the inner header loop ran to end of file, the outer
                // `getline` then failed too, and `getNextSpectrum_` returned
                // false. A block truncated before its first peak line is
                // dropped without an error.
                None => return Ok(None),
            }
        }
    }

    /// Read the body of one block. `Ok(None)` reports end of input before any
    /// peak line, which the source treats as a clean end of file.
    fn read_block(
        &mut self,
        spectrum: &mut MSSpectrum,
        sequences: &mut Vec<String>,
    ) -> Result<Option<()>> {
        loop {
            if !self.next_line()? {
                return Ok(None);
            }
            let text = trim(&self.line);
            if text.is_empty() {
                continue;
            }
            // Source `isdigit(line[0])`: only an ASCII digit starts the peak
            // list. The source passes a possibly negative `char` to `isdigit`,
            // which is undefined behaviour for non-ASCII input; this compares
            // the byte instead.
            let peak_line = text.as_bytes().first().is_some_and(u8::is_ascii_digit);
            if peak_line {
                self.read_peaks(spectrum)?;
                return Ok(Some(()));
            }
            self.read_header_line(spectrum, sequences)?;
        }
    }

    /// The peak sub-loop, entered with `self.line` holding a digit-leading
    /// line and returning after the block's `END IONS` is consumed.
    fn read_peaks(&mut self, spectrum: &mut MSSpectrum) -> Result<()> {
        let max_peaks = self.options.limits.max_peaks;
        let mz_range = self.options.mz_range.clone();
        let intensity_range = self.options.intensity_range.clone();
        loop {
            let line_number = self.input.line;
            let text = simplify(trim(&self.line));
            if !text.is_empty() {
                // Source `split(line, ' ', split, false)` returns false when the
                // line holds no separator at all, the only case the source
                // reports as "does not contain m/z and intensity".
                if !text.contains(' ') {
                    return Err(parse_error(
                        line_number,
                        format!(
                            "The content {text:?} does not contain m/z and intensity values separated by whitespace (space or tab)!"
                        ),
                    ));
                }
                let mut fields = text.split(' ');
                let mz_text = fields.next().unwrap_or("");
                let intensity_text = fields.next().unwrap_or("");
                // A third field is the optional per-peak charge, which the
                // source parses and then discards; further fields are ignored.
                let convert = |value: &str| -> Result<f64> {
                    source_double(value).ok_or_else(|| {
                        parse_error(
                            line_number,
                            format!(
                                "The content {text:?} could not be converted to a number! Expected two (m/z int) or three (m/z int charge) numbers separated by whitespace (space or tab)."
                            ),
                        )
                    })
                };
                let mz = convert(mz_text)?;
                let intensity = convert(intensity_text)? as f32;
                if !mz.is_finite() || !intensity.is_finite() {
                    return Err(parse_error(
                        line_number,
                        "MGF peak m/z and intensity must be finite",
                    ));
                }
                increment(&mut self.peaks, max_peaks, "MGF peak limit exceeded")?;
                if includes(&mz_range, mz) && includes(&intensity_range, f64::from(intensity)) {
                    push_peak(spectrum, Peak1D::new(mz, intensity))?;
                }
            }
            if !self.next_line()? {
                return Err(parse_error(
                    self.input.line,
                    "Reached end of file. Found \"BEGIN IONS\" but not the corresponding \"END IONS\"!",
                ));
            }
            if trim(&self.line) == "END IONS" {
                return Ok(());
            }
        }
    }

    /// One `KEY=value` header line inside a block.
    ///
    /// Prefix tests run in the source's order and each assumes `=` immediately
    /// after the key, so `ADDUCT=`, `ION_MODE=` and any other key are silently
    /// ignored. An `END IONS` line reaching here — a block with no peak line —
    /// is ignored too, so such a block merges into the following one.
    fn read_header_line(
        &self,
        spectrum: &mut MSSpectrum,
        sequences: &mut Vec<String>,
    ) -> Result<()> {
        let line = self.input.line;
        let text = trim(&self.line);
        if text.starts_with("PEPMASS") {
            let value = value_after(text, 8, line, "PEPMASS")?;
            // The source substitutes tab for space here but does *not*
            // simplify, so `PEPMASS=500  10` splits into three fields — the
            // middle one empty — and is a parse error. Only the peak lines get
            // the whitespace collapsing.
            let substituted = value.replace('\t', " ");
            let fields = source_split(&substituted, ' ');
            match fields.len() {
                1 => {
                    let mz = finite_double(fields[0], line, "precursor m/z")?;
                    precursor_mut(spectrum).mz = mz;
                }
                2 => {
                    let mz = finite_double(fields[0], line, "precursor m/z")?;
                    let intensity = finite_double(fields[1], line, "precursor intensity")? as f32;
                    if !intensity.is_finite() {
                        return Err(parse_error(line, "MGF precursor intensity exceeds f32"));
                    }
                    let precursor = precursor_mut(spectrum);
                    precursor.mz = mz;
                    precursor.intensity = intensity;
                }
                n => {
                    return Err(parse_error(
                        line,
                        format!(
                            "Cannot parse PEPMASS in {text:?} (expected 1 or 2 entries, but {n} were present)!"
                        ),
                    ));
                }
            }
        } else if text.starts_with("CHARGE") {
            // The source removes every '+' and then requires a complete i32, so
            // "2+" and "+2" parse while a charge list ("1,2,3"), a trailing '-'
            // ("2-") or any other trailing text raises a conversion error that
            // propagates out of `load`.
            let value = value_after(text, 7, line, "CHARGE")?.replace('+', "");
            let charge = source_int32(&value).ok_or_else(|| {
                parse_error(
                    line,
                    format!("Could not convert string {value:?} to an integer value"),
                )
            })?;
            precursor_mut(spectrum).charge = charge;
        } else if text.starts_with("RTINSECONDS") {
            let value = value_after(text, 12, line, "RTINSECONDS")?;
            spectrum.rt = finite_double(value, line, "retention time")?;
        } else if text.starts_with("TITLE") {
            read_title(spectrum, text, line)?;
        } else if text.starts_with("NAME") {
            set_meta(spectrum, MSM_METABOLITE_NAME, text, 5, line, "NAME")?;
        } else if text.starts_with("COMPOUND_NAME") {
            set_meta(
                spectrum,
                MSM_METABOLITE_NAME,
                text,
                14,
                line,
                "COMPOUND_NAME",
            )?;
        } else if text.starts_with("INCHI=") {
            set_meta(spectrum, MSM_INCHI_STRING, text, 6, line, "INCHI")?;
        } else if text.starts_with("SMILES") {
            set_meta(spectrum, MSM_SMILES_STRING, text, 7, line, "SMILES")?;
        } else if text.starts_with("IONMODE") {
            set_meta(spectrum, "IONMODE", text, 8, line, "IONMODE")?;
        } else if text.starts_with("MSLEVEL") {
            let value = value_after(text, 8, line, "MSLEVEL")?;
            match source_stoi(value) {
                // The source assigns whatever `std::stoi` returned. A
                // non-positive MS level makes the record invalid for every
                // consumer, so it is refused here instead of stored.
                Some(Ok(level)) => {
                    spectrum.ms_level =
                        u32::try_from(level)
                            .ok()
                            .filter(|&v| v > 0)
                            .ok_or_else(|| {
                                parse_error(line, format!("invalid MGF MS level {level}"))
                            })?;
                }
                // `std::invalid_argument`: the source falls back to MS2 and
                // records the fallback as a meta value.
                None => {
                    spectrum.ms_level = 2;
                    spectrum
                        .metadata
                        .insert("MSLEVEL".to_owned(), MetaValue::from("2"));
                }
                // `std::out_of_range`: the source falls back to MS2 silently.
                Some(Err(())) => spectrum.ms_level = 2,
            }
        } else if text.starts_with("SOURCE_INSTRUMENT") {
            set_meta(
                spectrum,
                "SOURCE_INSTRUMENT",
                text,
                18,
                line,
                "SOURCE_INSTRUMENT",
            )?;
        } else if text.starts_with("ORGANISM") {
            set_meta(spectrum, "ORGANISM", text, 9, line, "ORGANISM")?;
        } else if text.starts_with("PI") {
            set_meta(spectrum, "PI", text, 3, line, "PI")?;
        } else if text.starts_with("DATACOLLECTOR") {
            set_meta(spectrum, "DATACOLLECTOR", text, 14, line, "DATACOLLECTOR")?;
        } else if text.starts_with("LIBRARYQUALITY") {
            set_meta(spectrum, "LIBRARYQUALITY", text, 15, line, "LIBRARYQUALITY")?;
        } else if text.starts_with("SPECTRUMID") {
            set_meta(spectrum, GNPS_SPECTRUM_ID_KEY, text, 11, line, "SPECTRUMID")?;
        } else if text.starts_with("SCANS=") {
            set_meta(spectrum, SCAN_ID_KEY, text, 6, line, "SCANS")?;
        } else if text.starts_with("SEQ=") {
            // Per the Mascot specification a query may carry several SEQ lines,
            // each an independent sequence filter, so the value is always a
            // string list even when only one line was present. The source
            // round-trips the whole list through the meta value on every line,
            // which is quadratic; the list is accumulated here instead and
            // stored once when the block ends.
            let value = value_after(text, 4, line, "SEQ")?.to_owned();
            if sequences.len() >= MAX_SEQ_ENTRIES {
                return Err(invalid("MGF SEQ list limit exceeded"));
            }
            sequences.push(value);
        }
        Ok(())
    }
}

impl<R: BufRead> Iterator for MascotGenericReader<R> {
    type Item = Result<MSSpectrum>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        loop {
            match self.next_block() {
                Ok(Some(spectrum)) => {
                    if includes(&self.options.rt_range, spectrum.rt) {
                        return Some(Ok(spectrum));
                    }
                }
                Ok(None) => {
                    self.finished = true;
                    return None;
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }
    }
}

/// Read every block into an experiment, with default options.
///
/// ```
/// use openms::format::mascot_generic;
/// let text = "BEGIN IONS\nTITLE=demo\nPEPMASS=500.25 1200\nCHARGE=2+\n\
///             RTINSECONDS=42.5\n100.0 10.0\n200.0 20.0\nEND IONS\n";
/// let experiment = mascot_generic::read(text.as_bytes())?;
/// assert_eq!(experiment.spectra.len(), 1);
/// let spectrum = &experiment.spectra[0];
/// assert_eq!(spectrum.precursors[0].charge, 2);
/// assert_eq!(spectrum.rt, 42.5);
/// assert_eq!(spectrum.peaks.len(), 2);
/// // The source appends the native ID so titles stay unique.
/// assert_eq!(spectrum.metadata["TITLE"].as_str()?, "demo_index=0");
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// See [`read_with_options`].
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read every block into an experiment.
///
/// Source `MascotGenericFile::load` throws `FileNotFound` before opening a
/// file; [`load_with_options`] is the path-taking entry point that reports
/// that. The source also calls `updateRanges()` afterwards, which this crate
/// computes on demand rather than caching.
///
/// The reader reproduces these source behaviours, all of them surprising:
///
/// - Every line outside a `BEGIN IONS` … `END IONS` block is skipped, so an
///   MGF parameter header written by [`MascotGenericFile::store`] and any
///   global `CHARGE=` default are not read back.
/// - A block with no peak line does not end at its `END IONS`: that line
///   matches no header prefix and is ignored, so the block continues into the
///   following one and the two merge into a single spectrum.
/// - A block truncated before its first peak line is dropped silently; one
///   truncated after its first peak line is a [`Parse`](crate::Error::Parse)
///   error.
/// - Once a peak line has been seen, every following line must be a peak or
///   blank, so a header line after the peak list is a parse error.
/// - A third whitespace field on a peak line is the optional per-peak charge;
///   the source parses and discards it, and so does this reader.
///
/// # Errors
///
/// Returns [`Parse`](crate::Error::Parse) for input the source rejects — a
/// peak line with no whitespace or an unparsable number, a `PEPMASS=` with
/// more than two fields, a `CHARGE=` that is not one signed integer, and a
/// block truncated after its peak list started — and
/// [`InvalidValue`](crate::Error::InvalidValue) when a
/// [`ReadOptions::limits`] ceiling is reached. Nothing is appended to the
/// returned experiment before the whole file has been read, so a failed read
/// yields no partial result.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    read_reporting(reader, options, &mut ProgressReporter::silent())
}

/// The reader, with the source's `setProgress(is.tellg())` after each block
/// is added (`MascotGenericFile.h:96-101`). A block outside
/// [`ReadOptions::rt_range`], a native filter, is not returned and makes no
/// call.
fn read_reporting(
    reader: impl BufRead,
    options: &ReadOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<MSExperiment> {
    let mut result = MSExperiment::default();
    let mut blocks = MascotGenericReader::new(reader, options)?;
    while let Some(spectrum) = blocks.next() {
        push_spectrum(&mut result, spectrum?)?;
        if progress.is_reporting() {
            progress.set(blocks.source_position()?)?;
        }
    }
    Ok(result)
}

/// Replace `destination` only on success; a failure leaves it unchanged.
///
/// # Errors
///
/// See [`read_with_options`].
pub fn read_into(reader: impl BufRead, destination: &mut MSExperiment) -> Result<()> {
    read_into_with_options(reader, destination, &ReadOptions::default())
}

/// Replace `destination` only on success; a failure leaves it unchanged.
///
/// # Errors
///
/// See [`read_with_options`].
pub fn read_into_with_options(
    reader: impl BufRead,
    destination: &mut MSExperiment,
    options: &ReadOptions,
) -> Result<()> {
    *destination = read_with_options(reader, options)?;
    Ok(())
}

/// Read a file, with default options.
///
/// # Errors
///
/// Returns [`Io`](crate::Error::Io) when the file cannot be opened, which is
/// where the source throws `Exception::FileNotFound`; otherwise see
/// [`read_with_options`].
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    load_with_options(path, &ReadOptions::default())
}

/// Read a file.
///
/// # Errors
///
/// See [`load`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<MSExperiment> {
    load_reporting(path, options, &mut ProgressReporter::silent())
}

/// Read a file, reporting progress to `logger` as source
/// `MascotGenericFile::load` does (`MascotGenericFile.h:74-104`).
///
/// The calls are the source's: once the file is open,
/// `startProgress(0, file size in bytes, "loading MGF")`; after each block is
/// added, `setProgress` with the source's `is.tellg()`, the bytes consumed so
/// far, which is **-1** after a last `END IONS` line with no newline, because
/// `std::getline` then set `eofbit` and `tellg` fails (the command backend
/// prints its `Invalid progress value '-1'` diagnostic there, as the Release
/// build does); and `endProgress()` after the last block. The result is the
/// one [`load_with_options`] returns, and so is every error: both run the same
/// code, whose calls go nowhere for [`load_with_options`].
///
/// A file that cannot be opened makes no call, as the source's
/// `FileNotFound` makes none. (The source opens an existing but unreadable
/// file as an empty map, with a range ending at -1; this refuses it, as
/// [`load`] does.) A failure after the start, including invalid `options`, a
/// native check, leaves the section open, as in the source, where the
/// exception bypasses `endProgress`: no `-- done` line is printed, the nesting
/// depth stays one level deeper, and a command backend of `logger` refuses its
/// next start.
///
/// # Errors
///
/// As [`load`], plus [`InvalidValue`](crate::Error::InvalidValue) for a file
/// larger than `i64::MAX` bytes, and the errors of the progress calls
/// ([`ProgressLogger::start_progress`] and its siblings).
pub fn load_with_progress(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    logger: &mut ProgressLogger,
) -> Result<MSExperiment> {
    load_reporting(path, options, &mut ProgressReporter::new(Some(logger)))
}

/// The source's `load`: its section, sized by the file, around the reader.
fn load_reporting(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    progress: &mut ProgressReporter<'_>,
) -> Result<MSExperiment> {
    let file = File::open(path)?;
    if progress.is_reporting() {
        let size = i64::try_from(file.metadata()?.len())
            .map_err(|_| invalid("MGF file size exceeds the progress range"))?;
        progress.start(0, size, "loading MGF")?;
    }
    let experiment = read_reporting(BufReader::new(file), options, progress)?;
    progress.end()?;
    Ok(experiment)
}

/// Stream blocks into a consumer, with default options.
///
/// # Errors
///
/// See [`consume_with_options`].
pub fn consume(reader: impl BufRead, consumer: &mut impl MSDataConsumer) -> Result<()> {
    consume_with_options(reader, consumer, &ReadOptions::default())
}

/// Stream blocks into an [`MSDataConsumer`], retaining one spectrum at a time.
///
/// Native addition: the source MGF adapter has no consumer interface, only the
/// whole-file `load`. The block count is unknown before the file is read, so
/// `set_expected_size` is called once with `(0, 0)` and
/// `set_experimental_settings` with the default settings, because MGF carries
/// none. A consumer returning [`ControlFlow::Break`] stops the read
/// successfully; earlier callback effects are not rolled back.
///
/// # Errors
///
/// See [`read_with_options`], plus any error the consumer itself returns.
pub fn consume_with_options(
    reader: impl BufRead,
    consumer: &mut impl MSDataConsumer,
    options: &ReadOptions,
) -> Result<()> {
    consumer.set_expected_size(0, 0)?;
    consumer.set_experimental_settings(&Default::default())?;
    for spectrum in MascotGenericReader::new(reader, options)? {
        let mut spectrum = spectrum?;
        if consumer.consume_spectrum(&mut spectrum)?.is_break() {
            return Ok(());
        }
    }
    Ok(())
}

/// What a whole store call did, replacing the source's console output.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct WriteReport {
    /// Spectra actually written as a `BEGIN IONS` block.
    pub written: usize,
    /// MS2 spectra skipped because their precursor m/z was exactly zero. The
    /// source prints "No precursor m/z information … skipping spectrum!" to
    /// standard output; this port counts and reports it instead.
    pub skipped_without_precursor: usize,
    /// Spectra skipped because their MS level is not 2. The source warns only
    /// for level 0 and drops every other non-MS2 level in silence.
    pub skipped_ms_level: usize,
    /// Messages the source sends to `cerr`, `cout` or `OPENMS_LOG_WARN`.
    pub warnings: Vec<String>,
}

/// Result of offering one spectrum to [`MascotGenericFile::write_spectrum`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpectrumOutcome {
    /// Whether a block was written. A precursor m/z of exactly zero skips it.
    pub written: bool,
    /// Messages the source writes to `cerr` or `cout`.
    pub warnings: Vec<String>,
}

/// C++ default `ostream <<` formatting of a `double`: `%g` with precision 6.
///
/// `writeHeader_` streams the tolerance parameters without changing the stream
/// flags, so `3.0` is written as `3` and `0.3` as `0.3`.
fn ostream_double(value: f64) -> String {
    ostream_g(value, 6)
}

/// `%g` with an explicit precision, i.e. `ostream <<` after `setprecision(n)`
/// while the stream is still in its default float format.
///
/// The compact spectrum writer needs precisions 5 and 3: the source sets
/// `fixed` only while composing a *generated* `TITLE=` line, so a spectrum that
/// already carries a `TITLE` meta value has its `PEPMASS=` and `RTINSECONDS=`
/// written in this significant-digit form instead.
fn ostream_g(value: f64, precision: u32) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    if value == 0.0 {
        return "0".into();
    }
    let digits = i32::try_from(precision.max(1)).unwrap_or(6);
    let mantissa_decimals = usize::try_from(digits - 1).unwrap_or(0);
    let scientific = format!("{value:.mantissa_decimals$e}");
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        return scientific;
    };
    let Ok(exponent) = exponent.parse::<i32>() else {
        return scientific;
    };
    if !(-4..digits).contains(&exponent) {
        let mantissa = strip_trailing_zeros(mantissa);
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.unsigned_abs())
    } else {
        let decimals = usize::try_from(digits - 1 - exponent).unwrap_or(0);
        strip_trailing_zeros(&format!("{value:.decimals$}"))
    }
}
fn strip_trailing_zeros(text: &str) -> String {
    if !text.contains('.') {
        return text.to_owned();
    }
    let trimmed = text.trim_end_matches('0');
    trimmed.strip_suffix('.').unwrap_or(trimmed).to_owned()
}

/// `precisionWrapper`: `StringUtils::toStr(value, true)`, 15 significant digits
/// for a `double`, always keeping at least one fractional digit.
fn precision_wrapper(value: f64) -> String {
    format_float(value, true)
}
/// `precisionWrapper` for a `float`: six significant digits.
fn precision_wrapper_f32(value: f32) -> String {
    format_float32(value, true)
}

/// CV accessions whose native-ID format is `scan=<integer>`.
const SCAN_ACCESSIONS: [&str; 6] = [
    "MS:1000768",
    "MS:1000769",
    "MS:1000771",
    "MS:1000772",
    "MS:1000776",
    "MS:1002818",
];
/// CV accessions whose native-ID format is `file=<integer>`.
const FILE_ACCESSIONS: [&str; 2] = ["MS:1000773", "MS:1000775"];

/// The digits of the last `<key><digits>` match in `text`.
///
/// The source's regex token iterator collects every match and then takes
/// `matches.back()` *before* converting it, so the last match wins even when
/// its digits do not fit an `Int` — an earlier, convertible match is not a
/// fallback.
fn last_keyed_digits<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    let bytes = text.as_bytes();
    let mut found = None;
    let mut start = 0;
    while let Some(offset) = text.get(start..).and_then(|rest| rest.find(key)) {
        let at = start + offset;
        let digits_at = at + key.len();
        let mut end = digits_at;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
        if end > digits_at {
            // ASCII digits only, so this range is a character boundary.
            found = text.get(digits_at..end);
        }
        start = at + 1;
    }
    found
}

/// The last maximal ASCII digit run in `text`, unconverted.
fn last_digits(text: &str) -> Option<&str> {
    let bytes = text.as_bytes();
    let mut found = None;
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index].is_ascii_digit() {
            let start = index;
            while index < bytes.len() && bytes[index].is_ascii_digit() {
                index += 1;
            }
            found = text.get(start..index);
        } else {
            index += 1;
        }
    }
    found
}

/// WIFF native IDs encode the scan number as `cycle * 1000 + experiment`. The
/// source refuses an experiment of 1000 or more, because the encoding collides.
///
/// `cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+)` is collected with *two*
/// subgroups, and only the final match's pair is examined: an earlier pair with
/// an experiment of 1000 or more is never seen. Boost's `\s` covers vertical
/// tab and form feed as well as the four common blanks.
///
/// # Errors
///
/// The experiment ceiling raises `Exception::InvalidValue` upstream, which
/// `extractScanNumber` does *not* catch — it catches only `ConversionError` —
/// so it aborts the whole store rather than writing a sentinel.
fn wiff_scan_number(native_id: &str) -> Result<(i32, Option<String>)> {
    let bytes = native_id.as_bytes();
    let mut last: Option<(&str, &str)> = None;
    let mut start = 0;
    while let Some(offset) = native_id
        .get(start..)
        .and_then(|rest| rest.find("cycle="))
        .map(|index| start + index)
    {
        let mut position = offset + "cycle=".len();
        let cycle_start = position;
        while position < bytes.len() && bytes[position].is_ascii_digit() {
            position += 1;
        }
        let cycle = native_id.get(cycle_start..position).unwrap_or("");
        let space_start = position;
        while position < bytes.len()
            && matches!(bytes[position], b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
        {
            position += 1;
        }
        let separated = position > space_start;
        let rest = native_id.get(position..).unwrap_or("");
        let experiment = rest.strip_prefix("experiment=").map_or("", |digits| {
            let end = digits
                .as_bytes()
                .iter()
                .take_while(|b| b.is_ascii_digit())
                .count();
            digits.get(..end).unwrap_or("")
        });
        if !cycle.is_empty() && separated && !experiment.is_empty() {
            last = Some((cycle, experiment));
        }
        start = offset + 1;
    }
    let Some((cycle, experiment)) = last else {
        return Ok((-1, Some(no_match_warning(native_id))));
    };
    let (Ok(cycle), Ok(experiment)) = (cycle.parse::<i32>(), experiment.parse::<i32>()) else {
        return Ok((
            -1,
            Some(format!(
                "Values: '{cycle}', '{experiment}' could not be converted to int in string. Native ID='{native_id}' accession='MS:1000770'"
            )),
        ));
    };
    if experiment >= 1000 {
        return Err(invalid(&format!(
            "The value of experiment is too large and can not be handled properly.: '{experiment}'"
        )));
    }
    // The source computes `cycle * 1000 + experiment` in `int`, which is signed
    // overflow for a large cycle; the sentinel is written instead of wrapping.
    match cycle
        .checked_mul(1000)
        .and_then(|value| value.checked_add(experiment))
    {
        Some(value) => Ok((value, None)),
        None => Ok((
            -1,
            Some(format!(
                "native_id '{native_id}' encodes a scan number that overflows a 32-bit integer."
            )),
        )),
    }
}
fn no_match_warning(native_id: &str) -> String {
    format!("native_id '{native_id}' is invalid. Could not extract scan number.")
}

/// Scan number for `native_id` under `accession`, plus any message the source
/// would have logged.
///
/// Returns the source's `-1` sentinel when nothing matches or the last match
/// does not convert, because `writeSpectrum` streams the returned integer
/// verbatim into `SCANS=`. The conversion is `toInt32`, so the result is an
/// `i32`, as `extractScanNumber`'s return type is.
///
/// `METADATA/SpectrumLookup.h` and `METADATA/SpectrumNativeIDParser.h` are
/// separate unported headers; only the accession table the MGF writer reaches
/// is reproduced, and `docs/MASCOT_GENERIC_SUPPORT.md` records that.
///
/// # Errors
///
/// See [`wiff_scan_number`]: a WIFF experiment of 1000 or more aborts the store.
fn extract_scan_number(native_id: &str, accession: &str) -> Result<(i32, Option<String>)> {
    let digits = if SCAN_ACCESSIONS.contains(&accession) {
        last_keyed_digits(native_id, "scan=")
    } else if FILE_ACCESSIONS.contains(&accession) {
        last_keyed_digits(native_id, "file=")
    } else if accession == "MS:1000774" {
        last_keyed_digits(native_id, "index=")
    } else if accession == "MS:1001508" {
        last_keyed_digits(native_id, "scanId=")
    } else if accession == "MS:1000777" {
        last_keyed_digits(native_id, "spectrum=")
    } else if accession == "MS:1001530" {
        last_digits(native_id)
    } else if accession == "MS:1000770" {
        return wiff_scan_number(native_id);
    } else {
        return Ok((
            -1,
            Some(format!(
                "native_id: {native_id} accession: {accession} Could not extract scan number - no valid native_id_type_accession was provided"
            )),
        ));
    };
    let Some(digits) = digits else {
        return Ok((-1, Some(no_match_warning(native_id))));
    };
    let Ok(value) = digits.parse::<i32>() else {
        return Ok((
            -1,
            Some(format!(
                "Value: '{digits}' could not be converted to int in string. Native ID='{native_id}'"
            )),
        ));
    };
    if accession != "MS:1000774" {
        return Ok((value, None));
    }
    // An `index=` native ID is one less than the scan number consumers such as
    // pepXML expect, so the source adds one — in `int` arithmetic, which is
    // signed overflow at `INT_MAX`. The sentinel is written instead.
    match value.checked_add(1) {
        Some(value) => Ok((value, None)),
        None => Ok((
            -1,
            Some(format!(
                "native_id '{native_id}' scan number {value} + 1 overflows a 32-bit integer."
            )),
        )),
    }
}

/// The accession the source reads from the experiment's first source file, or
/// the [`UNKNOWN_NATIVE_ID_TYPE`] sentinel when there is none or it is empty.
fn native_id_accession(experiment: &MSExperiment) -> (String, Option<String>) {
    match experiment.settings.source_files.first() {
        None => (
            UNKNOWN_NATIVE_ID_TYPE.to_owned(),
            Some("MascotGenericFile: no native ID accession.".to_owned()),
        ),
        Some(file) if file.native_id_type_accession.is_empty() => (
            UNKNOWN_NATIVE_ID_TYPE.to_owned(),
            Some("MascotGenericFile: empty native ID accession.".to_owned()),
        ),
        Some(file) => (file.native_id_type_accession.clone(), None),
    }
}

/// `std::regex_replace(path.stem(), "[^a-zA-Z0-9]", "")`: the file-name stem
/// with every character that is not an ASCII letter or digit removed. A
/// multi-byte character consists only of non-ASCII bytes, so dropping the bytes
/// and dropping the characters give the same result.
fn filtered_filename(filename: &str) -> String {
    let stem = PathBuf::from(filename)
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    stem.chars().filter(char::is_ascii_alphanumeric).collect()
}

/// Every UniMod-backed modification identifier, as
/// `ModificationsDB::getAllSearchModifications` returns it: the full IDs of
/// records carrying a UniMod record ID, sorted case-insensitively with the
/// shorter of two otherwise-equal names first.
fn all_search_modifications(database: &ModificationsDB) -> Vec<String> {
    let mut names: Vec<String> = database
        .entries()
        .iter()
        .filter(|record| record.record_id().is_some_and(|id| id > 0))
        .map(|record| record.full_id().to_owned())
        .collect();
    // The source comparator lower-cases character by character and falls back
    // to the shorter string; this key is the same order and a total one, which
    // Rust's sort requires.
    names.sort_by_cached_key(|name| (name.to_lowercase(), name.len()));
    names
}

/// `updateMembers_`: expand `special_modifications` into the per-residue map.
fn special_modification_groups(parameters: &Param) -> Result<BTreeMap<String, String>> {
    let special = match parameters.value("special_modifications")? {
        ParamValue::String(value) => value.clone(),
        ParamValue::Empty => String::new(),
        _ => return Err(invalid("special_modifications must be a string")),
    };
    let groups = crate::data_structures::list::create::<String>(&special, b',')?;
    if groups.len() > MAX_MODIFICATIONS {
        return Err(invalid("MGF special modification list limit exceeded"));
    }
    let mut map = BTreeMap::new();
    for group in &groups {
        // `prefix(group, ' ')` is the whole string when there is no space, and
        // `prefix(suffix(group, '('), ')')` the residues between the last '('
        // and the first ')' after it, or the whole suffix when either is absent.
        let name = group.split(' ').next().unwrap_or("");
        let after = match group.rfind('(') {
            Some(position) => group.get(position + 1..).unwrap_or(""),
            None => group.as_str(),
        };
        let residues = match after.find(')') {
            Some(position) => after.get(..position).unwrap_or(""),
            None => after,
        };
        for residue in residues.chars() {
            map.insert(format!("{name} ({residue})"), group.clone());
        }
    }
    Ok(map)
}

/// Read/write Mascot generic files (MGF).
///
/// This is the port of the source class: a `DefaultParamHandler` whose
/// parameters are the Mascot search settings written as the MGF parameter
/// header, plus the peak-list writer. Reading needs none of that state and is
/// available as the free [`read`] and [`load`] functions; [`Self::load`] exists
/// so the source API maps one to one.
///
/// # Notes
///
/// The source class also derives from `ProgressLogger` and reports progress
/// while loading and storing. Progress reporting is not threaded through this
/// port; see `docs/MASCOT_GENERIC_SUPPORT.md`.
#[derive(Clone, Debug)]
pub struct MascotGenericFile {
    handler: DefaultParamHandler,
    store_compact: bool,
    mod_group_map: BTreeMap<String, String>,
}

impl MascotGenericFile {
    /// A writer with the source defaults, taking modification identifiers from
    /// the shared [`ModificationsDB`].
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) if the parameter
    /// tree exceeds its resource ceilings, which the modification valid-string
    /// lists dominate, or if a modification identifier contains a comma, which
    /// a `Param` valid-string list cannot represent.
    pub fn new() -> Result<Self> {
        Self::with_modifications(ModificationsDB::global())
    }

    /// A writer whose modification valid-string lists come from `database`.
    ///
    /// The source reads the process-wide `ModificationsDB` singleton; passing
    /// the registry in keeps a caller-owned override possible, as elsewhere in
    /// this crate.
    ///
    /// # Errors
    ///
    /// See [`Self::new`].
    pub fn with_modifications(database: &ModificationsDB) -> Result<Self> {
        let mut handler = DefaultParamHandler::new("MascotGenericFile")?;
        handler.set_defaults(defaults(database)?)?;
        let (groups, _) = handler.defaults_to_parameters_with(special_modification_groups)?;
        Ok(Self {
            handler,
            store_compact: false,
            mod_group_map: groups,
        })
    }

    /// Current parameters.
    pub fn parameters(&self) -> &Param {
        self.handler.parameters()
    }

    /// The source defaults, as `DefaultParamHandler::getDefaults` returns them.
    pub fn defaults(&self) -> &Param {
        self.handler.defaults()
    }

    /// Replace the parameters and rebuild the specificity-group map.
    ///
    /// This is `setParameters` together with the source's `updateMembers_`
    /// override. Unknown keys come back as warnings; a type or restriction
    /// violation is an error that leaves the previous parameters in place.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) for a parameter
    /// that violates its declared type or restriction, or when a resource
    /// ceiling is reached.
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (groups, warnings) = self
            .handler
            .set_parameters_with(parameters, special_modification_groups)?;
        self.mod_group_map = groups;
        Ok(warnings)
    }

    /// Specificity groups Mascot needs spelled as a group, keyed by the
    /// per-residue identifier OpenMS uses.
    ///
    /// `updateMembers_` expands `special_modifications`, so `Deamidated (NQ)`
    /// becomes the two entries `Deamidated (N)` and `Deamidated (Q)`, both
    /// mapping back to `Deamidated (NQ)`. [`Self::write_header_to`] rewrites
    /// each configured modification through this map.
    pub fn special_modification_groups(&self) -> &BTreeMap<String, String> {
        &self.mod_group_map
    }

    /// Whether the last store call requested the compact format.
    ///
    /// The source keeps this as the protected member `store_compact_`, assigned
    /// by `store` and read by `writeSpectrum`.
    pub fn store_compact(&self) -> bool {
        self.store_compact
    }

    /// Read an MGF file, as the source's template `load` member does.
    ///
    /// # Errors
    ///
    /// See [`load`].
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MSExperiment> {
        load(path)
    }

    /// [`Self::load`], reporting progress to `logger` as the source's `load`
    /// does; see the free [`load_with_progress`].
    ///
    /// # Errors
    ///
    /// See [`load_with_progress`].
    pub fn load_with_progress(
        &self,
        path: impl AsRef<Path>,
        logger: &mut ProgressLogger,
    ) -> Result<MSExperiment> {
        load_with_progress(path, &ReadOptions::default(), logger)
    }

    /// Write `experiment` to `path`.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when `path` does
    /// not end in `.mgf`, which is where the source throws
    /// `Exception::UnableToCreateFile`, or when a spectrum holds
    /// [`MAX_WRITTEN_PEAKS`] peaks or more; [`Io`](crate::Error::Io) when the
    /// file cannot be created, where the source throws
    /// `Exception::FileNotWritable`. The whole output is measured before the
    /// file is created, so a rejected store leaves no truncated file behind.
    pub fn store(
        &mut self,
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
        compact: bool,
    ) -> Result<WriteReport> {
        self.store_with_options(path, experiment, compact, &WriteOptions::default())
    }

    /// Write `experiment` to `path` under explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// See [`Self::store`].
    pub fn store_with_options(
        &mut self,
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
        compact: bool,
        options: &WriteOptions,
    ) -> Result<WriteReport> {
        self.store_reporting(
            path.as_ref(),
            experiment,
            compact,
            options,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`Self::store_with_options`], reporting progress to `logger` as the
    /// source's `store` does (`MascotGenericFile.cpp:458-476`).
    ///
    /// The calls are the source's: after the parameter header, the HTTP
    /// opening and the native-ID accession are written or determined,
    /// `startProgress(0, spectra, "storing mascot generic file")`;
    /// `setProgress(i)` before spectrum `i`, whether or not it is written; and
    /// `endProgress()` after the HTTP closing, before the file is flushed. The
    /// written bytes, the report and every error are those of
    /// [`Self::store_with_options`], which runs the same code with the calls
    /// going nowhere; its checks, including the measured dry run, precede the
    /// file and every call. With `internal:content` set to `header_only` the
    /// source writes no peak list and makes no call, and neither does this. A
    /// failure after the start leaves the section open, as described at
    /// [`load_with_progress`].
    ///
    /// # Errors
    ///
    /// See [`Self::store`], plus the errors of the progress calls.
    pub fn store_with_progress(
        &mut self,
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
        compact: bool,
        options: &WriteOptions,
        logger: &mut ProgressLogger,
    ) -> Result<WriteReport> {
        self.store_reporting(
            path.as_ref(),
            experiment,
            compact,
            options,
            &mut ProgressReporter::new(Some(logger)),
        )
    }

    /// The path `store`: checks, a measured dry run that reports nothing, then
    /// the file, whose peak lists report to `progress`.
    fn store_reporting(
        &mut self,
        path: &Path,
        experiment: &MSExperiment,
        compact: bool,
        options: &WriteOptions,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<WriteReport> {
        let name = path.to_string_lossy().into_owned();
        if !super::file_types::has_valid_extension(&name, super::FileType::Mgf) {
            return Err(invalid(
                "invalid file extension, expected 'mgf' for a Mascot generic file",
            ));
        }
        options.validate()?;
        let mut counter = Counter(options.max_output_bytes);
        self.render(
            &mut counter,
            &name,
            experiment,
            compact,
            &mut ProgressReporter::silent(),
        )?;
        let mut writer = BufWriter::new(File::create(path)?);
        let report = self.render(&mut writer, &name, experiment, compact, progress)?;
        writer.flush()?;
        Ok(report)
    }

    /// Write `experiment` to a stream, noting `filename` inside the file.
    ///
    /// `filename` is what the default `TITLE=` line and the HTTP peak-list
    /// enclosure record; Mascot echoes it back in its response. Only the
    /// file-name stem, stripped of every non-alphanumeric character, reaches
    /// the titles.
    ///
    /// # Errors
    ///
    /// See [`Self::store`]. No extension or writability check applies to a
    /// caller-provided stream, as in the source.
    pub fn store_to(
        &mut self,
        writer: impl Write,
        filename: &str,
        experiment: &MSExperiment,
        compact: bool,
    ) -> Result<WriteReport> {
        let mut writer = writer;
        self.render(
            &mut writer,
            filename,
            experiment,
            compact,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`Self::store_to`], reporting progress to `logger` with the calls
    /// [`Self::store_with_progress`] describes; the source's stream `store`
    /// makes them too.
    ///
    /// # Errors
    ///
    /// See [`Self::store_to`], plus the errors of the progress calls.
    pub fn store_to_with_progress(
        &mut self,
        writer: impl Write,
        filename: &str,
        experiment: &MSExperiment,
        compact: bool,
        logger: &mut ProgressLogger,
    ) -> Result<WriteReport> {
        let mut writer = writer;
        self.render(
            &mut writer,
            filename,
            experiment,
            compact,
            &mut ProgressReporter::new(Some(logger)),
        )
    }

    fn render(
        &mut self,
        writer: &mut impl Write,
        filename: &str,
        experiment: &MSExperiment,
        compact: bool,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<WriteReport> {
        self.store_compact = compact;
        let content = self.string_parameter("internal:content")?;
        let mut report = WriteReport::default();
        if content != "peaklist_only" {
            self.write_header_to(writer)?;
        }
        if content != "header_only" {
            self.write_experiment(writer, filename, experiment, &mut report, progress)?;
        }
        Ok(report)
    }

    /// The peak lists, with the source's progress section around them
    /// (`MascotGenericFile.cpp:458-476`).
    fn write_experiment(
        &self,
        writer: &mut impl Write,
        filename: &str,
        experiment: &MSExperiment,
        report: &mut WriteReport,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<()> {
        let enclosure = self.http_peak_list_enclosure(filename)?;
        let http = self.string_parameter("internal:HTTP_format")? == "true";
        if http {
            write!(writer, "{}", enclosure.0)?;
        }
        let stem = filtered_filename(filename);
        let (accession, warning) = native_id_accession(experiment);
        report.warnings.extend(warning);
        // The `fixed` flag the compact writer sets lives in the C++ ostream, so
        // it is sticky for the rest of the file once any spectrum has set it.
        let mut fixed = false;
        progress.start_count(experiment.spectra.len(), "storing mascot generic file")?;
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            progress.set_count(index)?;
            match spectrum.ms_level {
                2 => {
                    let outcome =
                        self.write_spectrum_with(writer, spectrum, &stem, &accession, &mut fixed)?;
                    report.warnings.extend(outcome.warnings);
                    if outcome.written {
                        report.written += 1;
                    } else {
                        report.skipped_without_precursor += 1;
                    }
                }
                0 => {
                    report.skipped_ms_level += 1;
                    report.warnings.push(
                        "MascotGenericFile: MSLevel is set to 0, ignoring this spectrum!"
                            .to_owned(),
                    );
                }
                // The source silently drops every other MS level.
                _ => report.skipped_ms_level += 1,
            }
        }
        if http {
            write!(writer, "{}", enclosure.1)?;
        }
        progress.end()
    }

    /// Strings that enclose the peak-list body for an HTTP submission.
    ///
    /// Returns the MIME part opening and the closing boundary, so custom peak
    /// content in another format — mzXML, say — can be embedded between them
    /// while this class writes only the parameter header. `filename` can later
    /// be found in the Mascot response.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when
    /// `internal:boundary` is missing or is not a string.
    pub fn http_peak_list_enclosure(&self, filename: &str) -> Result<(String, String)> {
        let boundary = self.string_parameter("internal:boundary")?;
        Ok((
            format!(
                "--{boundary}\nContent-Disposition: form-data; name=\"FILE\"; filename=\"{filename}\"\n\n"
            ),
            format!("\n\n--{boundary}--\n"),
        ))
    }

    /// Write one spectrum as a `BEGIN IONS` … `END IONS` block.
    ///
    /// `filename` is the already-filtered stem the default `TITLE=` line
    /// carries, and `native_id_type_accession` selects how `SCANS=` is derived
    /// from the native ID: the [`UNKNOWN_NATIVE_ID_TYPE`] sentinel takes the
    /// text after the native ID's last `=`, and anything else goes through the
    /// accession table `SpectrumLookup::extractScanNumber` uses, whose failure
    /// sentinel `-1` is written verbatim.
    ///
    /// A spectrum whose precursor m/z is exactly zero is skipped and
    /// [`SpectrumOutcome::written`] is then `false`. Only the first precursor
    /// is used; further precursors produce a warning. A stored `TITLE` meta
    /// value is written as-is, because it was either parsed from an MGF or set
    /// to be written to one; otherwise the title is composed from precursor
    /// m/z, retention time, native ID and `filename`.
    ///
    /// In compact form peak m/z values carry five fixed decimals and peak
    /// intensities three, and zero-intensity peaks are omitted. Otherwise every
    /// value is written at full precision.
    ///
    /// `PEPMASS=` and `RTINSECONDS=` are the exception: the source sets the
    /// stream's `fixed` flag only in the branch that *generates* a `TITLE=`
    /// line, so a compact spectrum that already carries a `TITLE` meta value
    /// gets five and three *significant* digits — `901.23` and `235`, not
    /// `901.23457` and `234.568`. Because the flag lives in the stream it stays
    /// set afterwards, so within one [`Self::store`] only the spectra before the
    /// first generated title or written peak line are affected. This entry
    /// point is one call on a fresh stream, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when the spectrum
    /// holds [`MAX_WRITTEN_PEAKS`] peaks or more — the source's guard against
    /// profile data — when precursor m/z or retention time is not finite, when
    /// a stored `TITLE` or `SEQ` meta value is not a string or string list, or
    /// when a WIFF native ID encodes an experiment of 1000 or more, which the
    /// source raises as an uncaught `InvalidValue`; [`Io`](crate::Error::Io) on
    /// a write failure.
    pub fn write_spectrum(
        &self,
        writer: &mut impl Write,
        spectrum: &MSSpectrum,
        filename: &str,
        native_id_type_accession: &str,
    ) -> Result<SpectrumOutcome> {
        let mut fixed = false;
        self.write_spectrum_with(
            writer,
            spectrum,
            filename,
            native_id_type_accession,
            &mut fixed,
        )
    }

    /// [`Self::write_spectrum`] with the caller's stream `fixed` flag, which
    /// the compact branch reads and then sets.
    fn write_spectrum_with(
        &self,
        writer: &mut impl Write,
        spectrum: &MSSpectrum,
        filename: &str,
        native_id_type_accession: &str,
        fixed: &mut bool,
    ) -> Result<SpectrumOutcome> {
        let mut warnings = Vec::new();
        if spectrum.precursors.len() > 1 {
            warnings.push(
                "Warning: The spectrum written to Mascot file has more than one precursor. The first precursor is used!"
                    .to_owned(),
            );
        }
        if spectrum.peaks.len() >= MAX_WRITTEN_PEAKS {
            return Err(invalid(&format!(
                "Spectrum to be written as MGF has {} peaks; the upper limit is 10,000. Only centroided data is allowed - this is most likely profile data.",
                spectrum.peaks.len()
            )));
        }
        let default_precursor = Precursor::default();
        let precursor = spectrum.precursors.first().unwrap_or(&default_precursor);
        let mz = precursor.mz;
        let rt = spectrum.rt;
        if mz == 0.0 {
            warnings.push(format!(
                "No precursor m/z information for spectrum with rt {} present, skipping spectrum!",
                ostream_double(rt)
            ));
            return Ok(SpectrumOutcome {
                written: false,
                warnings,
            });
        }
        if !mz.is_finite() || !rt.is_finite() {
            return Err(invalid("nonfinite MGF precursor m/z or retention time"));
        }
        writeln!(writer)?;
        writeln!(writer, "BEGIN IONS")?;
        let title = spectrum.metadata.get(TITLE_KEY);
        let (mz_text, rt_text) = if self.store_compact {
            if title.is_none() {
                // The source streams `fixed` while composing the generated
                // title, so the flag is already set for this spectrum's
                // PEPMASS and RTINSECONDS.
                *fixed = true;
            }
            if *fixed {
                (format!("{mz:.5}"), format!("{rt:.3}"))
            } else {
                (ostream_g(mz, 5), ostream_g(rt, 3))
            }
        } else {
            (precision_wrapper(mz), precision_wrapper(rt))
        };
        match title {
            Some(value) => writeln!(writer, "TITLE={}", value.as_str()?)?,
            None => writeln!(
                writer,
                "TITLE={mz_text}_{rt_text}_{}_{filename}",
                spectrum.native_id
            )?,
        }
        writeln!(writer, "PEPMASS={mz_text}")?;
        writeln!(writer, "RTINSECONDS={rt_text}")?;
        if native_id_type_accession == UNKNOWN_NATIVE_ID_TYPE {
            let scans = match spectrum.native_id.rfind('=') {
                Some(position) => spectrum.native_id.get(position + 1..).unwrap_or(""),
                None => spectrum.native_id.as_str(),
            };
            writeln!(writer, "SCANS={scans}")?;
        } else {
            let (scans, warning) =
                extract_scan_number(&spectrum.native_id, native_id_type_accession)?;
            warnings.extend(warning);
            writeln!(writer, "SCANS={scans}")?;
        }
        if precursor.charge != 0 && self.string_parameter("skip_spectrum_charges")? != "true" {
            let sign = if precursor.charge < 0 { '-' } else { '+' };
            writeln!(writer, "CHARGE={}{sign}", precursor.charge)?;
        }
        if let Some(value) = spectrum.metadata.get(SEQ_KEY) {
            for sequence in value.as_string_list()? {
                writeln!(writer, "SEQ={sequence}")?;
            }
        }
        for peak in &spectrum.peaks {
            if self.store_compact {
                if peak.intensity == 0.0 {
                    continue;
                }
                // The source streams `fixed` here too, on every peak line.
                *fixed = true;
                writeln!(writer, "{:.5} {:.3}", peak.mz, peak.intensity)?;
            } else {
                writeln!(
                    writer,
                    "{} {}",
                    precision_wrapper(peak.mz),
                    precision_wrapper_f32(peak.intensity)
                )?;
            }
        }
        writeln!(writer, "END IONS")?;
        Ok(SpectrumOutcome {
            written: true,
            warnings,
        })
    }

    /// Write the full MGF parameter header.
    ///
    /// Every line is a `KEY=value` pair, or a MIME part when
    /// `internal:HTTP_format` is `true`. `FORMAT` deliberately stays within the
    /// first five lines: that is how OpenMS recognises its own MGF files when
    /// the suffix is not `.mgf`. `COM` is omitted for an empty `search_title`
    /// and `USEREMAIL` for an empty `email`; `DECOY` appears only when `decoy`
    /// is `true`; `REPORT` is `AUTO` when `number_of_hits` is zero.
    ///
    /// This is the source's protected `writeHeader_`. It is public here because
    /// `internal:content = header_only` makes it the whole output of
    /// [`Self::store_to`], so a caller writing its own peak list needs it.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidValue`](crate::Error::InvalidValue) when a parameter is
    /// missing or has the wrong type, and [`Io`](crate::Error::Io) on a write
    /// failure.
    pub fn write_header_to(&self, writer: &mut impl Write) -> Result<()> {
        let title = self.string_parameter("search_title")?;
        if !title.is_empty() {
            self.write_parameter_header(writer, "COM")?;
            writeln!(writer, "{title}")?;
        }
        self.write_parameter_header(writer, "USERNAME")?;
        writeln!(writer, "{}", self.string_parameter("username")?)?;
        let email = self.string_parameter("email")?;
        if !email.is_empty() {
            self.write_parameter_header(writer, "USEREMAIL")?;
            writeln!(writer, "{email}")?;
        }
        self.write_parameter_header(writer, "FORMAT")?;
        writeln!(writer, "{}", self.string_parameter("internal:format")?)?;
        self.write_parameter_header(writer, "TOLU")?;
        writeln!(
            writer,
            "{}",
            self.string_parameter("precursor_error_units")?
        )?;
        self.write_parameter_header(writer, "ITOLU")?;
        writeln!(writer, "{}", self.string_parameter("fragment_error_units")?)?;
        self.write_parameter_header(writer, "FORMVER")?;
        writeln!(writer, "1.01")?;
        self.write_parameter_header(writer, "DB")?;
        writeln!(writer, "{}", self.string_parameter("database")?)?;
        if self.string_parameter("decoy")? == "true" {
            self.write_parameter_header(writer, "DECOY")?;
            writeln!(writer, "1")?;
        }
        self.write_parameter_header(writer, "SEARCH")?;
        writeln!(writer, "{}", self.string_parameter("search_type")?)?;
        self.write_parameter_header(writer, "REPORT")?;
        let hits = self.integer_parameter("number_of_hits")?;
        if hits == 0 {
            writeln!(writer, "AUTO")?;
        } else {
            writeln!(writer, "{hits}")?;
        }
        self.write_parameter_header(writer, "CLE")?;
        writeln!(writer, "{}", self.string_parameter("enzyme")?)?;
        self.write_parameter_header(writer, "MASS")?;
        writeln!(writer, "{}", self.string_parameter("mass_type")?)?;
        self.write_modifications(writer, "fixed_modifications", "MODS")?;
        self.write_modifications(writer, "variable_modifications", "IT_MODS")?;
        self.write_parameter_header(writer, "INSTRUMENT")?;
        writeln!(writer, "{}", self.string_parameter("instrument")?)?;
        self.write_parameter_header(writer, "PFA")?;
        writeln!(writer, "{}", self.integer_parameter("missed_cleavages")?)?;
        self.write_parameter_header(writer, "TOL")?;
        let precursor_tolerance = self.float_parameter("precursor_mass_tolerance")?;
        writeln!(writer, "{}", ostream_double(precursor_tolerance))?;
        self.write_parameter_header(writer, "ITOL")?;
        let fragment_tolerance = self.float_parameter("fragment_mass_tolerance")?;
        writeln!(writer, "{}", ostream_double(fragment_tolerance))?;
        self.write_parameter_header(writer, "TAXONOMY")?;
        writeln!(writer, "{}", self.string_parameter("taxonomy")?)?;
        self.write_parameter_header(writer, "CHARGE")?;
        writeln!(writer, "{}", self.string_parameter("charges")?)?;
        Ok(())
    }

    /// One parameter name, as `NAME=` or as a MIME part header.
    fn write_parameter_header(&self, writer: &mut impl Write, name: &str) -> Result<()> {
        if self.string_parameter("internal:HTTP_format")? == "true" {
            let boundary = self.string_parameter("internal:boundary")?;
            write!(
                writer,
                "--{boundary}\nContent-Disposition: form-data; name=\"{name}\"\n\n"
            )?;
        } else {
            write!(writer, "{name}=")?;
        }
        Ok(())
    }

    /// One `MODS=`/`IT_MODS=` line per modification, after rewriting each
    /// through [`Self::special_modification_groups`] and de-duplicating.
    ///
    /// The source collects into a `std::set<std::string>`, so the output is
    /// sorted by identifier and a group named twice appears once.
    fn write_modifications(&self, writer: &mut impl Write, key: &str, tag: &str) -> Result<()> {
        let list = match self.parameters().value(key)? {
            ParamValue::StringList(values) => values.clone(),
            ParamValue::String(value) => {
                crate::data_structures::list::create::<String>(value, b',')?
            }
            ParamValue::Empty => Vec::new(),
            _ => return Err(invalid("modification parameter must be a string list")),
        };
        if list.len() > MAX_MODIFICATIONS {
            return Err(invalid("MGF modification list limit exceeded"));
        }
        let mut filtered = BTreeSet::new();
        for modification in &list {
            filtered.insert(
                self.mod_group_map
                    .get(modification)
                    .cloned()
                    .unwrap_or_else(|| modification.clone()),
            );
        }
        for modification in &filtered {
            self.write_parameter_header(writer, tag)?;
            writeln!(writer, "{modification}")?;
        }
        Ok(())
    }

    fn string_parameter(&self, key: &str) -> Result<String> {
        match self.parameters().value(key)? {
            ParamValue::String(value) => Ok(value.clone()),
            ParamValue::Empty => Ok(String::new()),
            _ => Err(invalid("MGF parameter must be a string")),
        }
    }
    fn integer_parameter(&self, key: &str) -> Result<i64> {
        match self.parameters().value(key)? {
            ParamValue::Integer(value) => Ok(*value),
            _ => Err(invalid("MGF parameter must be an integer")),
        }
    }
    fn float_parameter(&self, key: &str) -> Result<f64> {
        match self.parameters().value(key)? {
            ParamValue::Float(value) => Ok(*value),
            #[allow(clippy::cast_precision_loss)]
            ParamValue::Integer(value) => Ok(*value as f64),
            _ => Err(invalid("MGF parameter must be a floating-point value")),
        }
    }
}

/// The source constructor's `defaults_`, including the `internal:` section that
/// is deliberately not shown to TOPP users.
fn defaults(database: &ModificationsDB) -> Result<Param> {
    let modifications = all_search_modifications(database);
    let advanced = [String::from("advanced")];
    let none: [String; 0] = [];
    let mut param = Param::new();
    let entries: [(&str, ParamValue, &str, &[String]); 24] = [
        (
            "database",
            "MSDB".into(),
            "Name of the sequence database",
            &none,
        ),
        (
            "search_type",
            "MIS".into(),
            "Name of the search type for the query",
            &advanced,
        ),
        (
            "enzyme",
            "Trypsin".into(),
            "The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).",
            &none,
        ),
        (
            "instrument",
            "Default".into(),
            "Instrument definition which specifies the fragmentation rules",
            &none,
        ),
        (
            "missed_cleavages",
            ParamValue::Integer(1),
            "Number of missed cleavages allowed for the enzyme",
            &none,
        ),
        (
            "precursor_mass_tolerance",
            ParamValue::Float(3.0),
            "Tolerance of the precursor peaks",
            &none,
        ),
        (
            "precursor_error_units",
            "Da".into(),
            "Units of the precursor mass tolerance",
            &none,
        ),
        (
            "fragment_mass_tolerance",
            ParamValue::Float(0.3),
            "Tolerance of the peaks in the fragment spectrum",
            &none,
        ),
        (
            "fragment_error_units",
            "Da".into(),
            "Units of the fragment peaks tolerance",
            &none,
        ),
        (
            "charges",
            "1,2,3".into(),
            "Charge states to consider, given as a comma separated list of integers (only used for spectra without precursor charge information)",
            &none,
        ),
        (
            "taxonomy",
            "All entries".into(),
            "Taxonomy specification of the sequences",
            &none,
        ),
        (
            "fixed_modifications",
            ParamValue::StringList(Vec::new()),
            "List of fixed modifications, according to UniMod definitions.",
            &none,
        ),
        (
            "variable_modifications",
            ParamValue::StringList(Vec::new()),
            "Variable modifications given as UniMod definitions.",
            &none,
        ),
        (
            "special_modifications",
            "Cation:Na (DE),Deamidated (NQ),Oxidation (HW),Phospho (ST),Sulfo (ST)".into(),
            "Modifications with specificity groups that are used by Mascot and have to be treated specially",
            &advanced,
        ),
        (
            "mass_type",
            "monoisotopic".into(),
            "Defines the mass type, either monoisotopic or average",
            &none,
        ),
        (
            "number_of_hits",
            ParamValue::Integer(0),
            "Number of hits which should be returned, if 0 AUTO mode is enabled.",
            &none,
        ),
        (
            "skip_spectrum_charges",
            "false".into(),
            "Sometimes precursor charges are given for each spectrum but are wrong, setting this to 'true' does not write any charge information to the spectrum, the general charge information is however kept.",
            &none,
        ),
        (
            "decoy",
            "false".into(),
            "Set to true if mascot should generate the decoy database.",
            &none,
        ),
        (
            "search_title",
            "OpenMS_search".into(),
            "Sets the title of the search.",
            &advanced,
        ),
        (
            "username",
            "OpenMS".into(),
            "Sets the username which is mentioned in the results file.",
            &advanced,
        ),
        (
            "email",
            "".into(),
            "Sets the email which is mentioned in the results file. Note: Some server require that a proper email is provided.",
            &none,
        ),
        (
            "internal:format",
            "Mascot generic".into(),
            "Sets the format type of the peak list, this should not be changed unless you write the header only.",
            &advanced,
        ),
        (
            "internal:boundary",
            "GZWgAaYKjHFeUaLOLEIOMq".into(),
            "MIME boundary for parameter header (if using HTTP format)",
            &advanced,
        ),
        (
            "internal:content",
            "all".into(),
            "Use parameter header + the peak lists with BEGIN IONS... or only one of them.",
            &advanced,
        ),
    ];
    for (key, value, description, tags) in entries {
        param.set_value(key, value, description, tags)?;
    }
    param.set_value(
        "internal:HTTP_format",
        "false".into(),
        "Write header with MIME boundaries instead of simple key-value pairs. For HTTP submission only.",
        &advanced,
    )?;
    let strings = |values: &[&str]| -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    };
    param.set_valid_strings("search_type", &strings(&["MIS", "SQ", "PMF"]))?;
    param.set_min_int("missed_cleavages", 0)?;
    param.set_min_float("precursor_mass_tolerance", 0.0)?;
    param.set_valid_strings(
        "precursor_error_units",
        &strings(&["%", "ppm", "mmu", "Da"]),
    )?;
    param.set_min_float("fragment_mass_tolerance", 0.0)?;
    param.set_valid_strings("fragment_error_units", &strings(&["mmu", "Da"]))?;
    param.set_valid_strings("fixed_modifications", &modifications)?;
    param.set_valid_strings("variable_modifications", &modifications)?;
    param.set_valid_strings("mass_type", &strings(&["monoisotopic", "average"]))?;
    param.set_min_int("number_of_hits", 0)?;
    param.set_valid_strings("skip_spectrum_charges", &strings(&["true", "false"]))?;
    param.set_valid_strings("decoy", &strings(&["true", "false"]))?;
    param.set_valid_strings(
        "internal:format",
        &strings(&["Mascot generic", "mzData (.XML)", "mzML (.mzML)"]),
    )?;
    param.set_valid_strings("internal:HTTP_format", &strings(&["true", "false"]))?;
    param.set_valid_strings(
        "internal:content",
        &strings(&["all", "peaklist_only", "header_only"]),
    )?;
    Ok(param)
}

/// Write an experiment with the source default parameters.
///
/// # Errors
///
/// See [`MascotGenericFile::store_to`].
pub fn write(
    writer: impl Write,
    filename: &str,
    experiment: &MSExperiment,
    compact: bool,
) -> Result<WriteReport> {
    MascotGenericFile::new()?.store_to(writer, filename, experiment, compact)
}

/// Store an experiment with the source default parameters.
///
/// # Errors
///
/// See [`MascotGenericFile::store`].
pub fn store(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    compact: bool,
) -> Result<WriteReport> {
    MascotGenericFile::new()?.store(path, experiment, compact)
}

/// A consumer that collects delivered spectra into an experiment.
///
/// Native helper for [`consume`]; the source has no consumer interface for MGF.
#[derive(Clone, Debug, Default)]
pub struct ExperimentCollector {
    /// Spectra delivered so far.
    pub experiment: MSExperiment,
}
impl MSDataConsumer for ExperimentCollector {
    fn set_expected_size(&mut self, spectra: usize, _chromatograms: usize) -> Result<()> {
        self.experiment
            .spectra
            .try_reserve(spectra)
            .map_err(|_| invalid("spectrum allocation failed"))
    }
    fn set_experimental_settings(
        &mut self,
        settings: &crate::metadata::ExperimentalSettings,
    ) -> Result<()> {
        self.experiment.settings = settings.clone();
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        push_spectrum(&mut self.experiment, std::mem::take(spectrum))?;
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(
        &mut self,
        _chromatogram: &mut crate::MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        Err(unsupported("MGF cannot store chromatograms"))
    }
}

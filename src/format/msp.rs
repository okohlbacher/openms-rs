// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! NIST peptide spectral libraries in MSP text, ported from
//! `FORMAT/MSPFile.h`.
//!
//! A library entry is a `Name:` line (`<peptide>/<charge>`, newer NIST
//! releases append `_<modifications>`), an optional `MW:` line, a `Comment:`
//! line of `key=value` tokens (among them `Mods=` and `Inst=`), a
//! `Num peaks:` (or SpectraST `NumPeaks:`) line and one line per peak:
//! m/z, intensity and a quoted annotation, separated by a tab, a carriage
//! return or two or more whitespace characters.
//!
//! ```text
//! Name: AAAAAAGAGPEMVR/2_2(0,A,Acetyl)(11,M,Oxidation)
//! MW: 1301.6400
//! Comment: Spec=Consensus Mods=2(0,A,Acetyl)(11,M,Oxidation) Inst=it
//! Num peaks: 2
//! 200.0   39      "a3/-0.12"
//! 300.1   258     "b3/-0.02,y2-H2O/-0.08"
//! ```
//!
//! Metabolite MSP libraries are read by [`crate::format::msp_generic`]
//! instead. See `docs/MSP_FILE_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::concept::constants::PROTON_MASS_U;
use crate::data_structures::string_utils::{to_double, to_float, to_int32};
use crate::format::file_info::text_format::ostream_g;
use crate::format::file_types::{FileType, has_valid_extension};
use crate::format::line_reader::{LineLimits, LineReader};
use crate::identification::{PeakAnnotation, PeptideHit, PeptideIdentification};
use crate::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor};
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Name of the source `DefaultParamHandler`.
const HANDLER_NAME: &str = "MSPFile";
/// Name of the string data array whose entries [`MSPFile::store`] writes as
/// peak annotations (source `Constants::UserParam::IonNames`).
pub const ION_NAMES: &str = "IonNames";
/// The source's fixed renaming of NIST modification names to OpenMS names
/// (`MSPFile.cpp:85-89`).
const MODNAME_TO_UNIMOD: [(&str, &str); 4] = [
    ("Pyro-glu", "Gln->pyro-Glu"),
    ("CAM", "Carbamidomethyl"),
    ("AB_old_ICATd8", "ICAT-D:2H(8)"),
    ("AB_old_ICATd0", "ICAT-D"),
];

/// Which spectra [`MSPFile::load`] keeps (parameter `instrument`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum InstrumentFilter {
    /// Every spectrum (the empty string).
    #[default]
    Any,
    /// Ion trap spectra, `Inst=it`.
    It,
    /// Q-TOF spectra, `Inst=qtof`.
    Qtof,
    /// TOF/TOF spectra, `Inst=toftof`.
    Toftof,
}

impl InstrumentFilter {
    /// The parameter spelling, also the `Inst=` value it keeps.
    pub fn name(self) -> &'static str {
        match self {
            Self::Any => "",
            Self::It => "it",
            Self::Qtof => "qtof",
            Self::Toftof => "toftof",
        }
    }

    fn from_name(name: &str) -> Result<Self> {
        match name {
            "" => Ok(Self::Any),
            "it" => Ok(Self::It),
            "qtof" => Ok(Self::Qtof),
            "toftof" => Ok(Self::Toftof),
            other => Err(Error::InvalidValue(format!(
                "instrument '{other}' is not one of '', 'it', 'qtof', 'toftof'"
            ))),
        }
    }
}

/// Resource ceilings for one MSP file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Line length and total input size.
    pub lines: LineLimits,
    /// Largest number of peaks summed over all kept spectra.
    pub max_peaks: usize,
    /// Largest number of `Name:` entries (peptide identifications) and of
    /// kept spectra.
    pub max_spectra: usize,
}

impl Default for ReadLimits {
    /// [`LineLimits::default`], 200 million peaks and 10 million spectra.
    fn default() -> Self {
        Self {
            lines: LineLimits::default(),
            max_peaks: 200_000_000,
            max_spectra: 10_000_000,
        }
    }
}

/// File adapter for NIST peptide spectral libraries in MSP text (source
/// `MSPFile`, a `DefaultParamHandler`).
///
/// The parameters are `parse_headers` (store the `Comment:` tokens as meta
/// values of each spectrum, default `false`), `parse_peakinfo` (store each
/// peak's annotation as a peak annotation of the peptide hit, default `true`),
/// `parse_firstpeakinfo_only` (only the first of several comma-separated peak
/// annotations, default `true`) and `instrument` (keep only spectra whose
/// `Inst=` matches; empty keeps all). The source stores the three flags as the
/// strings `"true"` and `"false"`; they are `bool` here and in
/// [`Self::to_param`] keep the source's string form.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MSPFile {
    /// Parameter `parse_headers`.
    pub parse_headers: bool,
    /// Parameter `parse_peakinfo`.
    pub parse_peakinfo: bool,
    /// Parameter `parse_firstpeakinfo_only`.
    pub parse_firstpeakinfo_only: bool,
    /// Parameter `instrument`.
    pub instrument: InstrumentFilter,
}

fn flag_entry(name: &str, value: bool, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::String(if value { "true" } else { "false" }.into()),
        valid_strings: vec!["true".into(), "false".into()],
        ..ParamEntry::default()
    }
}

impl MSPFile {
    /// A reader with the default parameters (source default constructor):
    /// headers not parsed, the first annotation of every peak parsed, every
    /// instrument kept.
    pub fn new() -> Self {
        Self {
            parse_headers: false,
            parse_peakinfo: true,
            parse_firstpeakinfo_only: true,
            instrument: InstrumentFilter::Any,
        }
    }

    /// The default parameters (installed by the source constructor,
    /// `MSPFile.cpp:25-39`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::new().to_param()
    }

    /// The current parameters in the layout of [`Self::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                flag_entry(
                    "parse_headers",
                    self.parse_headers,
                    "Flag whether header information should be parsed an stored for each spectrum",
                ),
                flag_entry(
                    "parse_peakinfo",
                    self.parse_peakinfo,
                    "Flag whether the peak annotation information should be parsed and stored for each peak",
                ),
                flag_entry(
                    "parse_firstpeakinfo_only",
                    self.parse_firstpeakinfo_only,
                    "Flag whether only the first or all peak annotation information should be parsed and stored for each peak.",
                ),
                ParamEntry {
                    name: "instrument".into(),
                    description: "If instrument given, only spectra of these type of instrument (Inst= in header) are parsed".into(),
                    value: ParamValue::String(self.instrument.name().into()),
                    valid_strings: vec![
                        String::new(),
                        "it".into(),
                        "qtof".into(),
                        "toftof".into(),
                    ],
                    ..ParamEntry::default()
                },
            ],
            nodes: Vec::new(),
        })
    }

    /// A reader from parameters (source `setParameters`). Missing parameters
    /// take their defaults; unknown ones are ignored, as the source only warns
    /// (see [`Self::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// is not one of its valid strings (the source's
    /// `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings
    /// the source logs.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(HANDLER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            let flag = |key: &str| -> Result<bool> { Ok(param.value(key)?.as_str()? == "true") };
            Ok(Self {
                parse_headers: flag("parse_headers")?,
                parse_peakinfo: flag("parse_peakinfo")?,
                parse_firstpeakinfo_only: flag("parse_firstpeakinfo_only")?,
                instrument: InstrumentFilter::from_name(param.value("instrument")?.as_str()?)?,
            })
        })
    }

    /// Load an MSP library (source `load(filename, ids, exp)`), with
    /// [`ReadLimits::default`].
    ///
    /// Returns one peptide identification per `Name:` entry that was kept and
    /// the experiment with one spectrum per kept entry, in file order; they
    /// replace whatever a previous call returned, as the source now resets
    /// both containers. The experiment records `filename` as its loaded file
    /// path and the content-detected type (MSP) as its loaded file type. The
    /// parsing rules are those of [`Self::read`].
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened or read (the source's
    /// `FileNotFound` and `FileNotReadable`), and the errors of [`Self::read`].
    pub fn load(
        &self,
        filename: impl AsRef<Path>,
    ) -> Result<(Vec<PeptideIdentification>, MSExperiment)> {
        self.load_with_limits(filename, ReadLimits::default())
    }

    /// [`Self::load`] under caller-chosen ceilings.
    ///
    /// # Errors
    ///
    /// As [`Self::load`].
    pub fn load_with_limits(
        &self,
        filename: impl AsRef<Path>,
        limits: ReadLimits,
    ) -> Result<(Vec<PeptideIdentification>, MSExperiment)> {
        let path = filename.as_ref();
        let file = File::open(path)?;
        let (ids, mut experiment) = self.read(BufReader::new(file), limits)?;
        let text = path
            .to_str()
            .ok_or_else(|| Error::InvalidValue("MSP path is not UTF-8".into()))?;
        experiment.settings.document.set_loaded_file_path(text)?;
        experiment.settings.document.set_loaded_file_type(path)?;
        Ok((ids, experiment))
    }

    /// Parse MSP text from `reader` (the body of source `load`).
    ///
    /// Lines end at LF only, as with `std::getline`; a CR before it stays in
    /// the line. Line by line:
    ///
    /// - `Name: <peptide>/<charge>[_...]` starts an entry: every innermost
    ///   `(...)` group is removed from the peptide, which the `Comment:`
    ///   line's `Mods=` describes instead, and a peptide identification with
    ///   one hit (score 0, rank 0, the charge) is appended.
    /// - `MW: <mass>` sets the pending spectrum's precursor m/z to
    ///   `(mass + charge * PROTON_MASS_U) / charge`.
    /// - `Comment:` is read token by token (single spaces). With an
    ///   `instrument` filter, the first `Inst=` token decides: a different
    ///   instrument drops the entry's identification and skips its peaks. A
    ///   `Mods=<n>/<pos>,<aa>,<name>/...` or `Mods=<n>(<pos>,<aa>,<name>)...`
    ///   token sets the modifications: position 0 is tried as an N-terminal,
    ///   the last position as a C-terminal modification, falling back to a
    ///   residue modification; `Pyro-glu`, `CAM`, `AB_old_ICATd8` and
    ///   `AB_old_ICATd0` are renamed to OpenMS names; a position outside the
    ///   peptide is ignored. With `parse_headers`, every `key=value` token
    ///   becomes a string meta value of the spectrum.
    /// - `Num peaks:` or `NumPeaks:` reads peak lines while they start with a
    ///   digit. m/z and intensity are read as `f32`. With `parse_peakinfo`, a
    ///   third, quoted field is required: `?...` becomes the annotation `?`
    ///   with charge 0, anything else is cut at its first space, split at
    ///   commas, each part cut at `/`, and `ion^charge` gives the charge
    ///   (default 1); `parse_firstpeakinfo_only` keeps only the first part.
    ///   The annotations replace those of the entry's hit, the spectrum gets
    ///   the native ID `index=<n>`, where `n` counts every `Num peaks:` block
    ///   including skipped ones, and is appended.
    ///
    /// The line that ends a peak block is consumed and not interpreted, as in
    /// the source: a `Name:` line directly after the peaks, with no empty line
    /// between, is lost. A spectrum whose instrument is filtered out keeps its
    /// pending precursor and meta values for the next entry, also as in the
    /// source.
    ///
    /// # Errors
    ///
    /// - [`Error::Parse`] when a `Name:` line lacks its space, `/` or charge,
    ///   a number does not convert (the source's `ConversionError`), a peak
    ///   line has no intensity, a requested annotation is missing (the
    ///   source's `ParseError`) or is not quoted, an annotation part has no ion
    ///   name, a `Mods=` token names more modifications than it lists or has
    ///   an entry without three fields, or an `MW:`, `Comment:` or peak block
    ///   precedes the first `Name:` line. Each of the last four is undefined
    ///   behaviour in the source, which indexes without checks.
    /// - [`Error::InvalidValue`] when an `MW:` line follows a charge below 1
    ///   (the source divides by it, or by its unsigned wrap-around), and the
    ///   errors of [`AASequence::parse`] and of setting a modification the
    ///   modification database does not know (the source's
    ///   `Exception::ElementNotFound`).
    /// - [`Error::InvalidValue`] when a limit is exceeded, and [`Error::Io`]
    ///   when reading fails.
    pub fn read(
        &self,
        reader: impl BufRead,
        limits: ReadLimits,
    ) -> Result<(Vec<PeptideIdentification>, MSExperiment)> {
        let mut lines = LineReader::lf_only(reader, limits.lines);
        let mut ids: Vec<PeptideIdentification> = Vec::new();
        let mut experiment = MSExperiment::new();
        let mut spectrum = MSSpectrum::default();
        let mut inst_type_correct = true;
        let mut spectrum_number = 0usize;
        let mut total_peaks = 0usize;
        let instrument = self.instrument.name();
        let mut line = String::new();

        while lines.next_line(&mut line)? {
            let number = lines.line_number();
            if line.starts_with("Name:") {
                if ids.len() >= limits.max_spectra {
                    return Err(Error::InvalidValue(format!(
                        "MSP entry limit of {} exceeded",
                        limits.max_spectra
                    )));
                }
                ids.push(name_entry(&line, number)?);
                inst_type_correct = true;
            } else if line.starts_with("MW:") {
                let fields = split(&line, ' ');
                if fields.len() == 2 {
                    let charge = first_hit(&mut ids, number)?.charge;
                    if charge < 1 {
                        return Err(Error::InvalidValue(format!(
                            "MSP line {number}: MW with charge {charge}"
                        )));
                    }
                    let mass = to_double(fields[1]).map_err(|message| parse(number, message))?;
                    let charge = f64::from(charge);
                    spectrum.precursors.resize_with(1, Precursor::default);
                    spectrum.precursors[0].mz = (mass + charge * PROTON_MASS_U) / charge;
                }
            } else if line.starts_with("Comment:") {
                for token in split(&line, ' ') {
                    if !inst_type_correct {
                        break;
                    }
                    if !instrument.is_empty() && token.starts_with("Inst=") {
                        if suffix(token, '=') != instrument {
                            inst_type_correct = false;
                            if ids.pop().is_none() {
                                return Err(before_name(number, "Comment:"));
                            }
                        }
                        break;
                    }
                    if token.starts_with("Mods=") && token != "Mods=0" {
                        apply_mods(token, first_hit(&mut ids, number)?, number)?;
                    }
                }
                if self.parse_headers && inst_type_correct {
                    parse_header(&line, &mut spectrum);
                }
            } else if line.starts_with("Num peaks:") || line.starts_with("NumPeaks:") {
                if inst_type_correct {
                    let mut annotations = Vec::new();
                    while lines.next_line(&mut line)? && starts_with_digit(&line) {
                        let number = lines.line_number();
                        let peak = self.peak_line(&line, number, &mut annotations)?;
                        total_peaks += 1;
                        if total_peaks > limits.max_peaks {
                            return Err(Error::InvalidValue(format!(
                                "MSP peak limit of {} exceeded",
                                limits.max_peaks
                            )));
                        }
                        spectrum.peaks.try_reserve(1).map_err(|_| {
                            Error::InvalidValue("MSP peak allocation failed".into())
                        })?;
                        spectrum.peaks.push(peak);
                    }
                    first_hit(&mut ids, number)?.peak_annotations = annotations;
                    if experiment.spectra.len() >= limits.max_spectra {
                        return Err(Error::InvalidValue(format!(
                            "MSP spectrum limit of {} exceeded",
                            limits.max_spectra
                        )));
                    }
                    let mut done = std::mem::take(&mut spectrum);
                    done.native_id = format!("index={spectrum_number}");
                    experiment.spectra.push(done);
                } else {
                    while lines.next_line(&mut line)? && starts_with_digit(&line) {}
                }
                spectrum_number += 1;
            }
        }
        Ok((ids, experiment))
    }

    /// One peak line: the peak, and its annotations appended to `annotations`.
    fn peak_line(
        &self,
        line: &str,
        number: usize,
        annotations: &mut Vec<PeakAnnotation>,
    ) -> Result<Peak1D> {
        let layout = || {
            parse(
                number,
                format!(
                    r#"'{line}' is not <mz><tab/spaces><intensity><tab/spaces>"<annotation>"<tab/spaces>"<comment>" in line {number}"#
                ),
            )
        };
        let fields = whitespace_fields(line);
        let mz = to_float(fields.first().ok_or_else(layout)?)
            .map_err(|message| parse(number, message))?;
        let intensity = to_float(fields.get(1).ok_or_else(layout)?)
            .map_err(|message| parse(number, message))?;
        let annotation = |name: &str, charge: i32| PeakAnnotation {
            annotation: name.into(),
            charge,
            mz: f64::from(mz),
            intensity: f64::from(intensity),
        };
        if self.parse_peakinfo {
            let Some(field) = fields.get(2) else {
                return Err(parse(
                    number,
                    format!(
                        "Requested reading peak info but no annotation found for line {number}"
                    ),
                ));
            };
            let text = unquote(field).ok_or_else(|| {
                parse(
                    number,
                    format!("'{field}' does not have the expected format of a quoted string"),
                )
            })?;
            if text.starts_with('?') {
                annotations.push(annotation("?", 0));
            } else {
                let text = prefix(&text, ' ');
                for part in split(text, ',') {
                    let ion = split(prefix(part, '/'), '^');
                    let Some(name) = ion.first() else {
                        return Err(parse(
                            number,
                            format!("peak annotation '{part}' has no ion name"),
                        ));
                    };
                    let charge = match ion.get(1) {
                        Some(charge) => {
                            to_int32(charge).map_err(|message| parse(number, message))?
                        }
                        None => 1,
                    };
                    annotations.push(annotation(name, charge));
                    if self.parse_firstpeakinfo_only {
                        break;
                    }
                }
            }
        }
        Ok(Peak1D::new(f64::from(mz), intensity))
    }

    /// Store spectra with their identifications as MSP text (source `store`,
    /// which takes an `AnnotatedMSRun`: `experiment.spectra[i]` pairs with
    /// `identifications[i]`). The layout is that of [`Self::write`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `filename` has a known extension other
    /// than `.msp` (the source's `UnableToCreateFile`), [`Error::Io`] when the
    /// file cannot be created or written (the source's `FileNotWritable`; the
    /// source then writes to an unchecked stream), and the errors of
    /// [`Self::write`].
    pub fn store(
        &self,
        filename: impl AsRef<Path>,
        experiment: &MSExperiment,
        identifications: &[PeptideIdentification],
    ) -> Result<()> {
        let path = filename.as_ref();
        let text = path.to_string_lossy();
        if !has_valid_extension(&text, FileType::Msp) {
            return Err(Error::InvalidValue(format!(
                "cannot create '{text}': invalid file extension, expected 'msp'"
            )));
        }
        let mut writer = BufWriter::new(File::create(path)?);
        self.write(experiment, identifications, &mut writer)?;
        writer.flush()?;
        Ok(())
    }

    /// Write spectra with their identifications as MSP text.
    ///
    /// Only pairs whose identification has hits are written, from its first
    /// hit, as:
    ///
    /// - `Name: <peptide>/<charge>`, where the peptide is the one-letter
    ///   sequence and only a modified `M` whose difference formula weighs
    ///   16 +- 0.01 Da is written, as `M(O)`;
    /// - `MW: <monoisotopic weight>`;
    /// - `Comment: Mods=<n>/<pos>,<aa>,<name>/... Inst=it` (or `Mods=0`),
    ///   where an N-terminal modification is position 0 and residue positions
    ///   count from 0; a C-terminal modification is not written;
    /// - `Num peaks: <n>` and one line per peak: m/z, the intensity scaled so
    ///   that the most intense peak is 10000 (unscaled when every intensity is
    ///   0 or below), and the quoted entry of the `IonNames` string data array,
    ///   or `"?"` without one; tab-separated;
    /// - an empty line.
    ///
    /// Numbers are written as by the source's `std::ostream`, with six
    /// significant digits in `%g` style (`1857.92`).
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when the two slices differ in length (the
    ///   source's paired iteration assumes they do not), or an `IonNames`
    ///   array is shorter than its spectrum (the source reads past it).
    /// - The errors of [`AASequence::mono_mass`] for a sequence without a
    ///   known weight, and [`Error::Io`] when writing fails.
    pub fn write(
        &self,
        experiment: &MSExperiment,
        identifications: &[PeptideIdentification],
        mut writer: impl Write,
    ) -> Result<()> {
        if experiment.spectra.len() != identifications.len() {
            return Err(Error::InvalidValue(format!(
                "{} spectra but {} peptide identifications",
                experiment.spectra.len(),
                identifications.len()
            )));
        }
        for (spectrum, identification) in experiment.spectra.iter().zip(identifications) {
            let Some(hit) = identification.hits.first() else {
                continue;
            };
            write_entry(spectrum, hit, &mut writer)?;
        }
        Ok(())
    }
}

/// The `Name:` line: the peptide with its charge.
fn name_entry(line: &str, number: usize) -> Result<PeptideIdentification> {
    let fields = split(line, ' ');
    let name = fields
        .get(1)
        .ok_or_else(|| parse(number, "Name line has no space-separated name".into()))?;
    let parts = split(name, '/');
    let (Some(peptide), Some(rest)) = (parts.first(), parts.get(1)) else {
        return Err(parse(
            number,
            format!("name '{name}' is not <peptide>/<charge>"),
        ));
    };
    let charge_text = split(rest, '_');
    let charge_text = charge_text
        .first()
        .ok_or_else(|| parse(number, format!("name '{name}' has no charge")))?;
    let charge = to_int32(charge_text).map_err(|message| parse(number, message))?;
    let peptide = remove_parenthesized(peptide);
    let sequence = AASequence::parse(&peptide)?;
    let mut identification = PeptideIdentification::default();
    identification
        .hits
        .push(PeptideHit::new(0.0, 0, charge, sequence)?);
    Ok(identification)
}

/// The first hit of the last identification (source
/// `ids.back().getHits()[0]`).
fn first_hit(ids: &mut [PeptideIdentification], number: usize) -> Result<&mut PeptideHit> {
    ids.last_mut()
        .and_then(|id| id.hits.first_mut())
        .ok_or_else(|| before_name(number, "an MW, Comment or peak block"))
}

fn before_name(number: usize, what: &str) -> Error {
    parse(number, format!("{what} precedes the first Name line"))
}

fn parse(line: usize, message: String) -> Error {
    Error::Parse { line, message }
}

/// A `Mods=` token applied to the hit's sequence.
fn apply_mods(token: &str, hit: &mut PeptideHit, number: usize) -> Result<()> {
    let mods = suffix(token, '=');
    let mut entries: Vec<&str> = split(mods, '/');
    if entries.len() <= 1 {
        // Mods=2(0,A,Acetyl)(11,M,Oxidation): the count, then groups.
        let count = prefix(mods, '(');
        let bound = to_int32(count).map_err(|message| parse(number, message))?;
        entries = vec![count];
        for group in parenthesized(mods) {
            // The source collects while `mod_split.size() - 1 <= sz`, with a
            // negative count wrapped to a huge `Size`.
            if bound >= 0 && entries.len() - 1 > bound as usize {
                break;
            }
            entries.push(&mods[group.0 + 1..group.1 - 1]);
        }
    }
    let count = to_int32(entries[0]).map_err(|message| parse(number, message))?;
    let count = usize::try_from(count)
        .map_err(|_| parse(number, format!("negative modification count in '{token}'")))?;
    let mut peptide = hit.sequence.clone();
    for index in 1..=count {
        let entry = entries.get(index).ok_or_else(|| {
            parse(
                number,
                format!("'{token}' names {count} modifications but lists fewer"),
            )
        })?;
        let fields = split(entry, ',');
        let (Some(position), Some(name)) = (fields.first(), fields.get(2)) else {
            return Err(parse(
                number,
                format!("modification '{entry}' is not <position>,<residue>,<name>"),
            ));
        };
        let name = MODNAME_TO_UNIMOD
            .iter()
            .find(|(nist, _)| nist == name)
            .map_or(*name, |(_, openms)| openms);
        let position = to_int32(position).map_err(|message| parse(number, message))?;
        // A negative position wraps to a huge `Size` in the source and is
        // ignored, as is any position past the last residue.
        let Ok(position) = usize::try_from(position) else {
            continue;
        };
        let last = peptide.len().wrapping_sub(1);
        if position > 0 && position < last {
            peptide.set_modification(position, name)?;
        } else if position == 0 {
            if peptide.set_n_terminal_modification(name).is_err() {
                peptide.set_modification(position, name)?;
            }
        } else if position == last && peptide.set_c_terminal_modification(name).is_err() {
            peptide.set_modification(position, name)?;
        }
    }
    hit.sequence = peptide;
    Ok(())
}

/// Source `parseHeader_`: every space-separated, trimmed `key=value` token
/// with exactly one `=` becomes a string meta value.
fn parse_header(header: &str, spectrum: &mut MSSpectrum) {
    for token in split(header, ' ') {
        let token = token.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r'));
        let parts = split(token, '=');
        if parts.len() == 2 {
            spectrum
                .metadata
                .insert(parts[0].to_owned(), MetaValue::from(parts[1]));
        }
    }
}

/// One MSP entry of [`MSPFile::write`].
fn write_entry(spectrum: &MSSpectrum, hit: &PeptideHit, writer: &mut impl Write) -> Result<()> {
    let sequence = &hit.sequence;
    let residues = sequence.as_str().as_bytes();
    let mut peptide = String::with_capacity(residues.len());
    let mut modifications = Vec::new();
    if let Some(modification) = sequence.n_terminal_modification() {
        let first = residues
            .first()
            .map_or(String::new(), |r| char::from(*r).to_string());
        modifications.push(format!("0,{first},{}", modification.name()));
    }
    for (index, &residue) in residues.iter().enumerate() {
        let residue = char::from(residue);
        match sequence.residue_modification(index)? {
            Some(modification) => {
                let oxidation = residue == 'M'
                    && modification
                        .diff_formula()
                        .is_ok_and(|formula| (formula.mono_mass() - 16.0).abs() < 0.01);
                if oxidation {
                    peptide.push_str("M(O)");
                } else {
                    peptide.push(residue);
                }
                modifications.push(format!("{index},{residue},{}", modification.name()));
            }
            None => peptide.push(residue),
        }
    }
    writeln!(writer, "Name: {peptide}/{}", hit.charge)?;
    writeln!(writer, "MW: {}", ostream_g(sequence.mono_mass()?, 6))?;
    write!(writer, "Comment:")?;
    if modifications.is_empty() {
        write!(writer, " Mods=0")?;
    } else {
        write!(
            writer,
            " Mods={}/{}",
            modifications.len(),
            modifications.join("/")
        )?;
    }
    writeln!(writer, " Inst=it")?;
    writeln!(writer, "Num peaks: {}", spectrum.peaks.len())?;

    let max_intensity = spectrum
        .peaks
        .iter()
        .map(|peak| f64::from(peak.intensity))
        .fold(0.0_f64, |max, value| if value > max { value } else { max });
    let ion_names = spectrum
        .string_data_arrays
        .iter()
        .find(|array| array.name == ION_NAMES);
    if let Some(names) = ion_names {
        if names.data.len() < spectrum.peaks.len() {
            return Err(Error::InvalidValue(format!(
                "IonNames array has {} entries for {} peaks",
                names.data.len(),
                spectrum.peaks.len()
            )));
        }
    }
    for (index, peak) in spectrum.peaks.iter().enumerate() {
        let intensity = if max_intensity == 0.0 {
            peak.intensity
        } else {
            (f64::from(peak.intensity) / max_intensity * 10000.0) as f32
        };
        let name = ion_names.map_or("?", |names| names.data[index].as_str());
        writeln!(
            writer,
            "{}\t{}\t\"{name}\"",
            ostream_g(peak.mz, 6),
            ostream_g(f64::from(intensity), 6)
        )?;
    }
    writeln!(writer)?;
    Ok(())
}

/// Source `StringUtils::split(s, char, out)`: an empty string yields no
/// fields; otherwise every separator ends a field, so consecutive separators
/// yield empty fields.
fn split(text: &str, separator: char) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// Source `StringUtils::prefix(s, char)`: the text before the first
/// `delimiter`, or all of it.
fn prefix(text: &str, delimiter: char) -> &str {
    text.split_once(delimiter).map_or(text, |(head, _)| head)
}

/// Source `StringUtils::suffix(s, char)`: the text after the last
/// `delimiter`, or all of it.
fn suffix(text: &str, delimiter: char) -> &str {
    text.rsplit_once(delimiter).map_or(text, |(_, tail)| tail)
}

/// Source `StringUtils::unquote(s)` with its default `"` and escape method:
/// `None` when `text` is not enclosed in double quotes, otherwise the inside
/// with `\"` and then `\\` unescaped.
fn unquote(text: &str) -> Option<String> {
    let inner = text
        .strip_prefix('"')
        .and_then(|rest| rest.strip_suffix('"'))
        .filter(|_| text.len() >= 2)?;
    Some(inner.replace("\\\"", "\"").replace("\\\\", "\\"))
}

/// `isdigit(line[0])` of the source's peak loop, on a non-empty line.
fn starts_with_digit(line: &str) -> bool {
    line.as_bytes().first().is_some_and(u8::is_ascii_digit)
}

/// The byte ranges of the source's `std::regex` `\((.*?)\)` matches, in order:
/// each `(` with the first `)` after it, where `.` matches neither LF nor CR.
fn parenthesized(text: &str) -> Vec<(usize, usize)> {
    let bytes = text.as_bytes();
    let mut matches = Vec::new();
    let mut start = 0;
    while start < bytes.len() {
        if bytes[start] == b'(' {
            let close = bytes[start + 1..]
                .iter()
                .position(|&b| matches!(b, b')' | b'\n' | b'\r'))
                .map(|offset| start + 1 + offset)
                .filter(|&end| bytes[end] == b')');
            if let Some(close) = close {
                matches.push((start, close + 1));
                start = close + 1;
                continue;
            }
        }
        start += 1;
    }
    matches
}

/// Source `std::regex_replace(peptide, rex, "")`: every
/// [`parenthesized`] match removed.
fn remove_parenthesized(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut from = 0;
    for (start, end) in parenthesized(text) {
        out.push_str(&text[from..start]);
        from = end;
    }
    out.push_str(&text[from..]);
    out
}

/// `std::isspace` in the "C" locale, the `\s` of the source's ECMAScript
/// expression.
fn c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The fields of a peak line: the source's `std::sregex_token_iterator` with
/// submatch -1 over `\s{2,}|\t|\r`. A separator is a run of two or more
/// whitespace characters (taken whole, the first alternative being greedy),
/// or a single tab or CR; the text before each separator is a field, and the
/// text after the last one is a field when it is not empty.
fn whitespace_fields(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut fields = Vec::new();
    let mut start = 0;
    let mut index = 0;
    while index < bytes.len() {
        let run = bytes[index..].iter().take_while(|&&b| c_space(b)).count();
        let separator = if run >= 2 {
            run
        } else if matches!(bytes[index], b'\t' | b'\r') {
            1
        } else {
            0
        };
        if separator > 0 {
            fields.push(&line[start..index]);
            index += separator;
            start = index;
        } else {
            index += 1;
        }
    }
    if start < bytes.len() {
        fields.push(&line[start..]);
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_fields_follow_the_source_expression() {
        assert_eq!(
            whitespace_fields("608.7\t974\t\"y6/0.39 2/2 2.4\""),
            vec!["608.7", "974", "\"y6/0.39 2/2 2.4\""]
        );
        assert_eq!(
            whitespace_fields("200.0   39      \"a3/-0.12\""),
            vec!["200.0", "39", "\"a3/-0.12\""]
        );
        // A single space is not a separator; a trailing CR is, and the empty
        // text after it is not a field.
        assert_eq!(whitespace_fields("1 2\t3\r"), vec!["1 2", "3"]);
        assert_eq!(whitespace_fields("1 \t2"), vec!["1", "2"]);
    }

    #[test]
    fn innermost_parentheses_are_removed() {
        assert_eq!(remove_parenthesized("AM(O)FDIK"), "AMFDIK");
        assert_eq!(remove_parenthesized("A(x)B(y)C"), "ABC");
        // Shortest match from the leftmost '(' that has one.
        assert_eq!(remove_parenthesized("((a))"), ")");
        assert_eq!(remove_parenthesized("A(B"), "A(B");
        assert_eq!(remove_parenthesized("(\r)(b)"), "(\r)");
    }

    #[test]
    fn unquote_requires_both_quotes() {
        assert_eq!(unquote("\"a\"").as_deref(), Some("a"));
        assert_eq!(unquote("\"\"").as_deref(), Some(""));
        assert_eq!(unquote("\"").as_deref(), None);
        assert_eq!(unquote("a").as_deref(), None);
        assert_eq!(unquote(r#""a\"b\\c""#).as_deref(), Some("a\"b\\c"));
    }

    #[test]
    fn string_helpers_follow_string_utils() {
        assert_eq!(split("", ' '), Vec::<&str>::new());
        assert_eq!(split("a  b", ' '), vec!["a", "", "b"]);
        assert_eq!(prefix("b8/-0.07 2/2", ' '), "b8/-0.07");
        assert_eq!(prefix("abc", '/'), "abc");
        assert_eq!(suffix("Mods=a=b", '='), "b");
        assert_eq!(suffix("abc", '='), "abc");
    }
}

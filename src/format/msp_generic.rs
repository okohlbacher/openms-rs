// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Metabolite spectral libraries in MSP text, ported from
//! `FORMAT/MSPGenericFile.h`.
//!
//! This reader is specialised for metabolite data. The required fields are
//! `Name`, `Num Peaks` and the peak data. The two numbers of a point may be
//! separated by a space, a tab or a colon; points may be separated by a space
//! or a semicolon. An example of the expected format:
//!
//! ```text
//! Name: foo
//! Num Peaks: 11
//! 35 310; 36 1230; 37 27; 38 303; 47 5240;
//! 66 203; 67 68; 68 77; 82 63; 83 240;
//! 136 350;
//! ```
//!
//! Another supported format:
//!
//! ```text
//! Name: bar
//! Num Peaks: 11
//! 35:310 36:1230 37:27 38:303 47:5240
//! 66:203 67:68 68:77 82:63 83:240
//! 136:350
//! ```
//!
//! See `docs/MSP_GENERIC_FILE_SUPPORT.md`.

use crate::concept::boost_regex::{BoostRegex, Captures, RegexOptions};
use crate::concept::constants::user_param::{
    MSM_CCS, MSM_INCHI_STRING, MSM_METABOLITE_NAME, MSM_PRECURSOR_ADDUCT, MSM_SMILES_STRING,
    MSM_SUM_FORMULA,
};
use crate::concept::fuzzy_string_comparator::format_g;
use crate::concept::string_conversion::{stod, stoul, to_double};
use crate::data_structures::list::ListFormat;
use crate::format::line_reader::{LineLimits, LineReader};
use crate::kernel::{MSExperiment, MSSpectrum, Peak1D, Precursor};
use crate::metadata::{MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Name of the source `DefaultParamHandler`.
const HANDLER_NAME: &str = "MSPGenericFile";
/// Meta value the source sets on the spectrum under construction: 1 while it
/// may be added to the library, 0 once added or while there is none.
pub const IS_VALID: &str = "is_valid";
/// Meta keys [`MSPGenericFile::store`] writes on lines of their own or not at
/// all.
const IGNORED_ON_STORE: [&str; 4] = ["Synon", "CAS#", "NIST#", "Num Peaks"];

/// Resource ceilings for one MSP file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Line length and total input size.
    pub lines: LineLimits,
    /// Largest number of peaks summed over all spectra, including those of
    /// spectra that are then rejected as duplicates.
    pub max_peaks: usize,
    /// Largest number of spectra in the library.
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

/// How [`MSPGenericFile::store_with_options`] spells numbers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NumberFormat {
    /// The shortest decimal text that reads back as the same value, never in
    /// exponent notation, so every written value can be read again. The
    /// native default.
    #[default]
    RoundTrip,
    /// The source's `std::ostream << double`: six significant digits in `%g`
    /// style. It loses precision beyond six digits, and a value from 1e6 up
    /// (or below 1e-4) is written in exponent notation, which the reader's
    /// point expression cannot read back: `100:1.5e+06` loads as the peak
    /// (100, 1.5) followed by the bogus peak (6, next m/z).
    Source,
}

/// Options for [`MSPGenericFile::store_with_options`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StoreOptions {
    /// Spelling of the retention time, the peak positions and intensities.
    pub numbers: NumberFormat,
}

impl StoreOptions {
    /// The source's `store`: numbers in six-digit `%g` style.
    pub fn source() -> Self {
        Self {
            numbers: NumberFormat::Source,
        }
    }
}

/// Load an MSP text file into an [`MSExperiment`], and store one (source
/// `MSPGenericFile`, a `DefaultParamHandler`).
///
/// The one parameter, `synonyms_separator` (default `|`), is the text placed
/// between the synonyms joined into a spectrum's `Synon` meta value, and the
/// text [`Self::store`] splits that value on.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MSPGenericFile {
    /// The source's `synonyms_separator_`.
    synonyms_separator: String,
    /// Names already added to the library during the current load (source
    /// `loaded_spectra_names_`), used to reject duplicates.
    loaded_spectra_names: BTreeSet<String>,
    /// Synonyms collected for the spectrum under construction (source
    /// `synonyms_`).
    synonyms: Vec<String>,
}

impl Default for MSPGenericFile {
    fn default() -> Self {
        Self {
            synonyms_separator: "|".into(),
            loaded_spectra_names: BTreeSet::new(),
            synonyms: Vec::new(),
        }
    }
}

/// The source's expressions, compiled once per load.
struct Expressions {
    name: BoostRegex,
    retention_time: BoostRegex,
    synon: BoostRegex,
    points_line: BoostRegex,
    point: BoostRegex,
    cas_nist: BoostRegex,
    precursor_mz: BoostRegex,
    num_peaks: BoostRegex,
    inchi: BoostRegex,
    smiles: BoostRegex,
    sum_formula: BoostRegex,
    precursor_type: BoostRegex,
    ccs: BoostRegex,
    metadatum: BoostRegex,
}

fn regex(pattern: &str, no_mod_s: bool, icase: bool) -> Result<BoostRegex> {
    BoostRegex::with_options(
        pattern,
        RegexOptions {
            icase,
            no_mod_s,
            ..RegexOptions::default()
        },
    )
}

impl Expressions {
    /// `MSPGenericFile.cpp:74-89`, with the same construction flags.
    fn new() -> Result<Self> {
        Ok(Self {
            name: regex("(?:^Name|^NAME): (.+)", true, false)?,
            retention_time: regex("(?:^Retention Time|^RETENTIONTIME): (.+)", true, false)?,
            synon: regex("^synon(?:yms?)?: (.+)", true, true)?,
            points_line: regex(r"^\d", false, false)?,
            point: regex(r"(\d+(?:\.\d+)?)[: \t](\d+(?:\.\d+)?);? ?", false, false)?,
            cas_nist: regex(r"^CAS#: ([\d-]+);  NIST#: (\d+)", false, false)?,
            precursor_mz: regex("^PRECURSORMZ: (.+)", false, false)?,
            num_peaks: regex("(?:^Num Peaks|^Num peaks): (\\d+)", false, false)?,
            inchi: regex("^INCHIKEY: (.+)", false, false)?,
            smiles: regex("^SMILES: (.+)", false, false)?,
            sum_formula: regex("^FORMULA: (.+)", false, false)?,
            precursor_type: regex("^PRECURSORTYPE: (.+)", false, false)?,
            ccs: regex("^CCS: (.+)", false, false)?,
            metadatum: regex("^(.+): (.+)", true, false)?,
        })
    }
}

/// One capture group as UTF-8 text; an unmatched group is empty, as a Boost
/// `sub_match` converts to an empty string.
fn group<'h>(line: &'h [u8], captures: &Captures, index: usize, number: usize) -> Result<&'h str> {
    let range = captures.get(index).unwrap_or(0..0);
    std::str::from_utf8(&line[range]).map_err(|_| Error::Parse {
        line: number,
        message: "MSP field is not valid UTF-8".into(),
    })
}

fn string_value(text: &str) -> MetaValue {
    MetaValue::from(text)
}

/// `StringUtils::split(s, std::string splitter, out)`: an empty string gives
/// no parts, an empty splitter gives one part per character.
fn split_on(text: &str, splitter: &str) -> Vec<String> {
    if text.is_empty() {
        Vec::new()
    } else if splitter.is_empty() {
        text.chars().map(String::from).collect()
    } else {
        text.split(splitter).map(str::to_owned).collect()
    }
}

impl MSPGenericFile {
    /// A reader with the default parameters (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// The class's default parameters (source `getDefaultParameters`): the
    /// string `synonyms_separator`, default `|`.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
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
            entries: vec![ParamEntry {
                name: "synonyms_separator".into(),
                description:
                    "The character that will separate the synonyms in the Synon metaValue.".into(),
                value: ParamValue::String(self.synonyms_separator.clone()),
                ..ParamEntry::default()
            }],
            nodes: Vec::new(),
        })
    }

    /// A reader from parameters, as source `setParameters` followed by
    /// `updateMembers_`. A missing `synonyms_separator` takes its default;
    /// unknown parameters are ignored, as the source only warns (see
    /// [`Self::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `synonyms_separator` is not a
    /// string (the source's `Exception::InvalidParameter`).
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
            Ok(Self {
                synonyms_separator: param.value("synonyms_separator")?.as_str()?.to_owned(),
                ..Self::default()
            })
        })
    }

    /// The separator placed between synonyms (source `synonyms_separator_`).
    pub fn synonyms_separator(&self) -> &str {
        &self.synonyms_separator
    }

    /// Load the MSP file at `filename` into a new library (source `load`).
    ///
    /// Equivalent to [`Self::load_with_limits`] with [`ReadLimits::default`].
    ///
    /// # Errors
    ///
    /// See [`Self::read_with_limits`]; also [`Error::Io`] when the file cannot
    /// be opened, where the source throws `Exception::FileNotFound`,
    /// `Exception::FileNotReadable` or `Exception::IOException` (the
    /// `NotFound`, `PermissionDenied` and other I/O error kinds here).
    pub fn load(&mut self, filename: impl AsRef<Path>) -> Result<MSExperiment> {
        self.load_with_limits(filename, ReadLimits::default())
    }

    /// [`Self::load`] under explicit `limits`.
    ///
    /// # Errors
    ///
    /// As [`Self::load`].
    pub fn load_with_limits(
        &mut self,
        filename: impl AsRef<Path>,
        limits: ReadLimits,
    ) -> Result<MSExperiment> {
        let file = File::open(filename.as_ref())?;
        self.read_with_limits(BufReader::new(file), limits)
    }

    /// Read MSP text from `reader` into a new library.
    ///
    /// # Errors
    ///
    /// See [`Self::read_with_limits`].
    pub fn read(&mut self, reader: impl BufRead) -> Result<MSExperiment> {
        self.read_with_limits(reader, ReadLimits::default())
    }

    /// Read MSP text from `reader` into a new library, under `limits`.
    ///
    /// Lines end at LF only, as `std::getline` reads them, so a CRLF file
    /// leaves a CR at the end of every line. Each line is tested against the
    /// source's expressions in this order, the first match deciding:
    ///
    /// 1. A line that starts with a digit is a peak line. Every point is
    ///    matched in turn, each search continuing after the previous match,
    ///    and becomes a peak (intensity narrowed to `f32`).
    /// 2. `synon:`, `synonym:` or `synonyms:` (any case): a synonym.
    /// 3. `Name: ` or `NAME: `: the spectrum built so far is offered to the
    ///    library (see [`Self::add_spectrum_to_library`]) and a new one starts
    ///    with this name, also stored as the `Metabolite_Name` meta value.
    /// 4. `Num Peaks: ` or `Num peaks: ` and digits: the `Num Peaks` meta
    ///    value, as a string.
    /// 5. `Retention Time: ` or `RETENTIONTIME: `: the retention time, read by
    ///    `std::stod` (the longest numeric prefix).
    /// 6. `PRECURSORMZ: `: one precursor at that m/z, read by `std::stod`.
    /// 7. `CAS#: <digits and dashes>;  NIST#: <digits>`: the `CAS#` and
    ///    `NIST#` meta values.
    /// 8. `INCHIKEY: `, `SMILES: `, `FORMULA: `, `PRECURSORTYPE: `: the
    ///    `Inchi_String`, `SMILES_String`, `Sum_Formula` and `Precursor_Ion`
    ///    meta values.
    /// 9. `CCS: `: the `CCS` meta value, trimmed, as a float when the source's
    ///    `StringUtils::toDouble` converts it and as the text otherwise.
    /// 10. `<key>: <value>`: any other meta value, as a string. The key runs
    ///     to the last `: ` of the line.
    ///
    /// The expressions of lines 2, 3, 5 and 10 use `no_mod_s`, so their value
    /// stops before a CR; the others do not, so on a CRLF file the `INCHIKEY`,
    /// `SMILES`, `FORMULA` and `PRECURSORTYPE` values keep a trailing CR. This
    /// is the source's behaviour and is kept.
    ///
    /// Anything before the first name belongs to no spectrum and is discarded.
    /// At the end of input the last spectrum is offered to the library. A
    /// spectrum without a retention time gets its index in the library.
    /// Both the names seen and the synonyms are reset at the start.
    ///
    /// # Errors
    ///
    /// - The errors of [`Self::add_spectrum_to_library`], which abort the load.
    /// - [`Error::Parse`] when a peak line holds no point or a number does not
    ///   convert: the source lets `std::stod`'s `std::invalid_argument` or
    ///   `std::out_of_range` escape, which is not an OpenMS exception.
    /// - [`Error::Parse`], a native check, when a peak, retention time or
    ///   precursor m/z is NaN or infinite, or an intensity overflows `f32`;
    ///   and when a stored field is not UTF-8, where the source keeps bytes.
    /// - [`Error::InvalidValue`] when the input exceeds `limits`.
    pub fn read_with_limits(
        &mut self,
        reader: impl BufRead,
        limits: ReadLimits,
    ) -> Result<MSExperiment> {
        let expressions = Expressions::new()?;
        self.loaded_spectra_names.clear();
        self.synonyms.clear();
        let mut library = MSExperiment::new();
        let mut spectrum = MSSpectrum::default();
        // To avoid adding invalid spectra to the library.
        spectrum
            .metadata
            .insert(IS_VALID.into(), MetaValue::from(0i64));
        let mut lines = LineReader::lf_only(reader, limits.lines);
        let mut line = Vec::new();
        let mut peaks = 0usize;

        while lines.next_bytes(&mut line)? {
            let number = lines.line_number();
            let x = &expressions;
            if x.points_line.is_search_match(&line)? {
                let mut offset = 0;
                let mut first = true;
                loop {
                    let Some(captures) = x.point.search(&line[offset..])? else {
                        if first {
                            return Err(Error::Parse {
                                line: number,
                                message: "stod: no conversion of '' (a peak line without a point)"
                                    .into(),
                            });
                        }
                        break;
                    };
                    first = false;
                    let slice = &line[offset..];
                    let position =
                        stod(group(slice, &captures, 1, number)?).map_err(|message| {
                            Error::Parse {
                                line: number,
                                message,
                            }
                        })?;
                    let intensity =
                        stod(group(slice, &captures, 2, number)?).map_err(|message| {
                            Error::Parse {
                                line: number,
                                message,
                            }
                        })? as f32;
                    if !position.is_finite() || !intensity.is_finite() {
                        return Err(Error::Parse {
                            line: number,
                            message: "MSP peak is NaN, infinite or outside the f32 intensity range"
                                .into(),
                        });
                    }
                    peaks += 1;
                    if peaks > limits.max_peaks {
                        return Err(Error::InvalidValue(format!(
                            "MSP input exceeds {} peaks",
                            limits.max_peaks
                        )));
                    }
                    spectrum
                        .peaks
                        .try_reserve(1)
                        .map_err(|_| Error::InvalidValue("MSP peak allocation failed".into()))?;
                    spectrum.peaks.push(Peak1D::new(position, intensity));
                    let end = captures.range().end;
                    // An empty match cannot occur (a point has two digits), so
                    // the search always advances.
                    offset += end.max(1);
                    if offset >= line.len() {
                        break;
                    }
                }
            } else if let Some(m) = x.synon.search(&line)? {
                let synonym = group(&line, &m, 1, number)?.to_owned();
                self.synonyms.push(synonym);
            } else if let Some(m) = x.name.search(&line)? {
                let name = group(&line, &m, 1, number)?.to_owned();
                if library.spectra.len() >= limits.max_spectra {
                    return Err(Error::InvalidValue(format!(
                        "MSP input exceeds {} spectra",
                        limits.max_spectra
                    )));
                }
                self.add_spectrum_to_library(&mut spectrum, &mut library)?;
                spectrum = MSSpectrum::default();
                self.synonyms.clear();
                spectrum.name = name.clone();
                spectrum
                    .metadata
                    .insert(MSM_METABOLITE_NAME.into(), string_value(&name));
                spectrum
                    .metadata
                    .insert(IS_VALID.into(), MetaValue::from(1i64));
            } else if let Some(m) = x.num_peaks.search(&line)? {
                let value = group(&line, &m, 1, number)?;
                spectrum
                    .metadata
                    .insert("Num Peaks".into(), string_value(value));
            } else if let Some(m) = x.retention_time.search(&line)? {
                let rt = stod(group(&line, &m, 1, number)?).map_err(|message| Error::Parse {
                    line: number,
                    message,
                })?;
                if !rt.is_finite() {
                    return Err(Error::Parse {
                        line: number,
                        message: "MSP retention time is NaN or infinite".into(),
                    });
                }
                spectrum.rt = rt;
            } else if let Some(m) = x.precursor_mz.search(&line)? {
                let mz = stod(group(&line, &m, 1, number)?).map_err(|message| Error::Parse {
                    line: number,
                    message,
                })?;
                if !mz.is_finite() {
                    return Err(Error::Parse {
                        line: number,
                        message: "MSP precursor m/z is NaN or infinite".into(),
                    });
                }
                spectrum.precursors = vec![Precursor::new(mz, 0)];
            } else if let Some(m) = x.cas_nist.search(&line)? {
                let cas = group(&line, &m, 1, number)?;
                let nist = group(&line, &m, 2, number)?;
                spectrum.metadata.insert("CAS#".into(), string_value(cas));
                spectrum.metadata.insert("NIST#".into(), string_value(nist));
            } else if let Some(m) = x.inchi.search(&line)? {
                let value = group(&line, &m, 1, number)?;
                spectrum
                    .metadata
                    .insert(MSM_INCHI_STRING.into(), string_value(value));
            } else if let Some(m) = x.smiles.search(&line)? {
                let value = group(&line, &m, 1, number)?;
                spectrum
                    .metadata
                    .insert(MSM_SMILES_STRING.into(), string_value(value));
            } else if let Some(m) = x.sum_formula.search(&line)? {
                let value = group(&line, &m, 1, number)?;
                spectrum
                    .metadata
                    .insert(MSM_SUM_FORMULA.into(), string_value(value));
            } else if let Some(m) = x.precursor_type.search(&line)? {
                let value = group(&line, &m, 1, number)?;
                spectrum
                    .metadata
                    .insert(MSM_PRECURSOR_ADDUCT.into(), string_value(value));
            } else if let Some(m) = x.ccs.search(&line)? {
                // Collision cross section, typically a plain number in square
                // angstrom; a non-numeric value is kept verbatim as a string.
                let value = crate::format::text::trim(group(&line, &m, 1, number)?);
                let stored = match to_double(value) {
                    Ok(ccs) => MetaValue::new(MetaValueData::Float(ccs))
                        .unwrap_or_else(|_| string_value(value)),
                    Err(_) => string_value(value),
                };
                spectrum.metadata.insert(MSM_CCS.into(), stored);
            } else if let Some(m) = x.metadatum.search(&line)? {
                // Other metadata; this matches everything and so comes last.
                let key = group(&line, &m, 1, number)?.to_owned();
                let value = group(&line, &m, 2, number)?;
                spectrum.metadata.insert(key, string_value(value));
            }
        }
        // To make sure a spectrum is added even if no empty line is present
        // before the end of input.
        self.add_spectrum_to_library(&mut spectrum, &mut library)?;
        Ok(library)
    }

    /// Validate `spectrum` and add a copy to `library` (the source's private
    /// `addSpectrumToLibrary`, reached from its class test through
    /// `MSPGenericFile_friend`).
    ///
    /// Nothing happens unless the spectrum's `is_valid` meta value is non-zero.
    /// A spectrum whose name was already added during this load is a duplicate
    /// and is skipped (the source logs `DUPLICATE: <name>`). Otherwise the
    /// synonyms collected since its name are joined with the separator into its
    /// `Synon` meta value, `is_valid` is removed, a negative retention time is
    /// replaced by the number of spectra already in `library`, and a copy is
    /// appended. These changes are made to `spectrum` itself, as in the
    /// source, so a caller that offers the same spectrum again sees them. In
    /// every case that gets past the first check, `is_valid` of `spectrum` is
    /// then set to 0.
    ///
    /// The source removes the last *character* after joining, which is only
    /// correct for a one-character separator: with `||` it leaves a trailing
    /// `|`, and with an empty separator it cuts the last synonym. This port
    /// joins with the separator, which is identical for the default `|`.
    ///
    /// # Errors
    ///
    /// - [`Error::MissingInformation`] when the name is empty ("The current
    ///   spectrum misses the Name information.") or, for a new name, the
    ///   `Num Peaks` meta value is missing ("The current spectrum misses the
    ///   Num Peaks information."); the source's `Exception::MissingInformation`.
    /// - [`Error::Parse`] when the number of peaks differs from `Num Peaks`
    ///   ("The number of points parsed does not coincide with `Num Peaks`."),
    ///   the source's `Exception::ParseError`, and when `Num Peaks` does not
    ///   convert by `std::stoul`.
    /// - [`Error::InvalidValue`] when `is_valid` is absent or not an integer,
    ///   where the source's `DataValue` conversion throws
    ///   `Exception::ConversionError`.
    ///
    /// The source checks the name before the duplicate test and `Num Peaks`
    /// after it, and so does this port. On an error neither `spectrum` nor
    /// `library` has changed, and the checks all precede the changes.
    pub fn add_spectrum_to_library(
        &mut self,
        spectrum: &mut MSSpectrum,
        library: &mut MSExperiment,
    ) -> Result<()> {
        let valid = match spectrum.metadata.get(IS_VALID).map(MetaValue::data) {
            Some(MetaValueData::Integer(value)) => *value,
            _ => {
                return Err(Error::InvalidValue(format!(
                    "the spectrum's '{IS_VALID}' meta value is not an integer"
                )));
            }
        };
        if valid as i32 == 0 {
            return Ok(());
        }
        // Check that the required metadata (Name, Num Peaks) is present; Num
        // Peaks is checked with the number of points parsed.
        if spectrum.name.is_empty() {
            return Err(Error::MissingInformation(
                "The current spectrum misses the Name information.".into(),
            ));
        }
        if !self.loaded_spectra_names.contains(&spectrum.name) {
            let num_peaks = spectrum.metadata.get("Num Peaks").ok_or_else(|| {
                Error::MissingInformation(
                    "The current spectrum misses the Num Peaks information.".into(),
                )
            })?;
            let num_peaks = num_peaks.to_list_text()?.into_owned();
            let expected =
                stoul(&num_peaks).map_err(|message| Error::Parse { line: 0, message })?;
            if spectrum.peaks.len() as u64 != expected {
                return Err(Error::Parse {
                    line: 0,
                    message: format!(
                        "{num_peaks}: The number of points parsed does not coincide with `Num Peaks`."
                    ),
                });
            }
            library
                .spectra
                .try_reserve(1)
                .map_err(|_| Error::InvalidValue("MSP library allocation failed".into()))?;
            // As in the source, these changes are made to `spectrum` itself
            // and stay visible to the caller after the copy is added.
            if !self.synonyms.is_empty() {
                spectrum.metadata.insert(
                    "Synon".into(),
                    string_value(&self.synonyms.join(&self.synonyms_separator)),
                );
            }
            spectrum.metadata.remove(IS_VALID);
            if spectrum.rt < 0.0 {
                // Set the retention time to the spectrum index.
                spectrum.rt = library.spectra.len() as f64;
            }
            self.loaded_spectra_names.insert(spectrum.name.clone());
            library.spectra.push(spectrum.clone());
        }
        spectrum
            .metadata
            .insert(IS_VALID.into(), MetaValue::from(0i64));
        Ok(())
    }

    /// Save the spectra and their metadata in MSP text, with the native
    /// [`StoreOptions::default`] (numbers that read back exactly).
    ///
    /// # Errors
    ///
    /// See [`Self::store_with_options`].
    pub fn store(&self, filename: impl AsRef<Path>, library: &MSExperiment) -> Result<()> {
        self.store_with_options(filename, library, StoreOptions::default())
    }

    /// Save `library` to `filename` in MSP text.
    ///
    /// Spectra without peaks are skipped. Each other spectrum is written as:
    /// `Name: <name>`; `Retention Time: <rt>`; one `Synon: <s>` line per part
    /// of a string `Synon` meta value split on the separator; `CAS#: <c>;
    /// NIST#: <n>` when both meta values exist; `<key>: <value>` for every
    /// other meta value except `Num Peaks`, in key order; `Num Peaks: <n>`;
    /// the peaks as `<mz>:<intensity> `, five to a line; and an empty line.
    /// Meta values are written as the source's `StringUtils::toStr` spells a
    /// `DataValue`.
    ///
    /// The source writes meta values in the order their names were first
    /// registered in the process-wide `MetaInfoRegistry`, which depends on the
    /// process history; this port writes them sorted by key.
    ///
    /// # Errors
    ///
    /// - [`Error::MissingInformation`] when any spectrum, with or without
    ///   peaks, has no name ("The current spectrum misses the Name
    ///   information."). All spectra are checked before the file is created;
    ///   the source throws after writing the spectra before the unnamed one.
    /// - [`Error::InvalidValue`], native checks, when a name, key or written
    ///   value contains a line break, which would split the record, and when
    ///   a retention time or peak is NaN or infinite.
    /// - [`Error::Io`] when the file cannot be created or written, the
    ///   source's `Exception::FileNotWritable`.
    pub fn store_with_options(
        &self,
        filename: impl AsRef<Path>,
        library: &MSExperiment,
        options: StoreOptions,
    ) -> Result<()> {
        let text = self.render(library, options)?;
        let mut writer = BufWriter::new(File::create(filename.as_ref())?);
        writer.write_all(text.as_bytes())?;
        writer.flush()?;
        Ok(())
    }

    /// Write `library` in the layout of [`Self::store_with_options`] to
    /// `writer`.
    ///
    /// # Errors
    ///
    /// As [`Self::store_with_options`].
    pub fn write_with_options(
        &self,
        mut writer: impl Write,
        library: &MSExperiment,
        options: StoreOptions,
    ) -> Result<()> {
        let text = self.render(library, options)?;
        writer.write_all(text.as_bytes())?;
        writer.flush()?;
        Ok(())
    }

    fn render(&self, library: &MSExperiment, options: StoreOptions) -> Result<String> {
        if library
            .spectra
            .iter()
            .any(|spectrum| spectrum.name.is_empty())
        {
            return Err(Error::MissingInformation(
                "The current spectrum misses the Name information.".into(),
            ));
        }
        let f64_text = |value: f64| -> Result<String> {
            if !value.is_finite() {
                return Err(Error::InvalidValue(
                    "MSP output value is NaN or infinite".into(),
                ));
            }
            Ok(match options.numbers {
                NumberFormat::RoundTrip => format!("{value}"),
                NumberFormat::Source => format_g(value),
            })
        };
        let f32_text = |value: f32| -> Result<String> {
            if !value.is_finite() {
                return Err(Error::InvalidValue(
                    "MSP output intensity is NaN or infinite".into(),
                ));
            }
            Ok(match options.numbers {
                NumberFormat::RoundTrip => format!("{value}"),
                NumberFormat::Source => format_g(f64::from(value)),
            })
        };
        let line = |text: &mut String, content: &str| -> Result<()> {
            if content.contains(['\n', '\r']) {
                return Err(Error::InvalidValue(format!(
                    "MSP output line would contain a line break: {content:?}"
                )));
            }
            text.push_str(content);
            text.push('\n');
            Ok(())
        };

        let mut text = String::new();
        for spectrum in &library.spectra {
            // Spectra without peaks are not stored.
            if spectrum.peaks.is_empty() {
                continue;
            }
            line(&mut text, &format!("Name: {}", spectrum.name))?;
            line(
                &mut text,
                &format!("Retention Time: {}", f64_text(spectrum.rt)?),
            )?;
            if let Some(synonyms) = spectrum.metadata.get("Synon") {
                if let MetaValueData::String(value) = synonyms.data() {
                    for synonym in split_on(value, &self.synonyms_separator) {
                        line(&mut text, &format!("Synon: {synonym}"))?;
                    }
                }
            }
            if let (Some(cas), Some(nist)) = (
                spectrum.metadata.get("CAS#"),
                spectrum.metadata.get("NIST#"),
            ) {
                line(
                    &mut text,
                    &format!(
                        "CAS#: {};  NIST#: {}",
                        cas.to_list_text()?,
                        nist.to_list_text()?
                    ),
                )?;
            }
            for (key, value) in &spectrum.metadata {
                if IGNORED_ON_STORE.contains(&key.as_str()) {
                    continue;
                }
                line(&mut text, &format!("{key}: {}", value.to_list_text()?))?;
            }
            line(&mut text, &format!("Num Peaks: {}", spectrum.peaks.len()))?;
            for (index, peak) in spectrum.peaks.iter().enumerate() {
                text.push_str(&f64_text(peak.mz)?);
                text.push(':');
                text.push_str(&f32_text(peak.intensity)?);
                text.push(' ');
                if (index + 1) % 5 == 0 {
                    text.push('\n');
                }
            }
            if spectrum.peaks.len() % 5 != 0 {
                text.push('\n');
            }
            // separator
            text.push('\n');
        }
        Ok(text)
    }
}

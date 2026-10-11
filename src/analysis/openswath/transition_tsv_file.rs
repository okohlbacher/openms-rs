// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reading and writing of OpenSWATH transition lists.
//!
//! Ports `ANALYSIS/OPENSWATH/TransitionTSVFile.h` (with its `.cpp`). The
//! transition lists can be either comma- or tab-separated plain text files
//! (CSV or TSV format); a semicolon is recognised as well. Modifications
//! should be provided in UniMod format (`PEPT(UniMod:21)IDE`), but can also be
//! provided in TPP format or as names OpenMS understands
//! (`PEPT(Phosphorylation)IDE(UniMod:27)A`). For another format that stores
//! transitions see `TransitionPQPFile`.
//!
//! The following columns are required (a later synonym is tried only when the
//! earlier names are absent or empty):
//!
//! | Column | Type | Meaning |
//! |---|---|---|
//! | `PrecursorMz` | float | precursor ion m/z |
//! | `ProductMz` | float; synonym `FragmentMz` | product (fragment) ion m/z |
//! | `LibraryIntensity` | float; synonyms `RelativeIntensity`, `RelativeFragmentIntensity` | relative intensity of the fragment ion |
//! | `NormalizedRetentionTime` | float; synonyms `RetentionTimeCalculatorScore`, `iRT`, `RetentionTime`, `Tr_recalibrated` | expected (normalized) retention time |
//!
//! For targeted proteomics these are recommended: `GeneName`; `ProteinId`
//! (synonym `ProteinName`, several separated by `;`); `PeptideSequence`
//! (synonyms `Sequence`, `StrippedSequence`; the sequence without
//! modifications); `ModifiedPeptideSequence` (synonyms
//! `FullUniModPeptideName`, `FullPeptideName`, `ModifiedSequence`);
//! `PrecursorCharge` (synonym `Charge`); `ProductCharge` (synonym
//! `FragmentCharge`); `FragmentType` (e.g. `b` or `y`; synonym
//! `FragmentIonType`); `FragmentSeriesNumber` (e.g. `7` for y7; synonyms
//! `FragmentNumber`, `FragmentIonOrdinal`).
//!
//! OpenSWATH groups transitions to detect candidate analyte signals. The
//! groups are derived from the input unless given: `TransitionGroupId`
//! (synonyms `TransitionGroupName`, `transition_group_id`); `TransitionId`
//! (unique per transition; synonyms `TransitionName`, `transition_name`);
//! `Decoy` (1 decoy, 0 target; synonyms `decoy`, `IsDecoy`);
//! `PeptideGroupLabel` (the label group, PSI-MS `MS:1000893`, of isotopically
//! labelled forms of one peptide); `DetectingTransition`,
//! `IdentifyingTransition`, `QuantifyingTransition` (0 or 1; also in snake
//! case; `Quantitative` is a synonym of the last).
//!
//! Optional columns that OpenSWATH does not use: `CollisionEnergy` (synonym
//! `CE`), `Annotation`, `UniprotId` (synonym `UniprotID`), `LabelType`,
//! `PrecursorIonMobility`, `FragmentMzDelta`, `FragmentModification` and
//! `Peptidoforms` (separated by `|`). For targeted metabolomics:
//! `CompoundName` (synonym `CompoundId` in the source's documentation only;
//! a row with a compound name that is not empty or `NA` is a compound),
//! `SMILES`, `SumFormula` and `Adducts`.
//!
//! Peptide label groups designate isotopically modified forms of the same
//! peptide: `PEPTIDEAK` and `PEPTIDEAK[+8]` share `PEPTIDEAK_gr1`, the
//! phosphorylated forms share `PEPTIDEAK_gr2`. A group whose members have
//! different sequences is most likely an error and, unless
//! `override_group_label_check` is set, each member of it is given its own
//! transition group id as label.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TRANSITION_TSV_SUPPORT.md)
//! for the member mapping and the evidence. The source contains no
//! `#pragma omp`; the port is serial. The source's `ProgressLogger` base is not
//! ported, and the `OpenSwath::LightTargetedExperiment` overloads wait for
//! `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`.

use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    Compound, Interpretation, Peptide, PeptideModification, Protein, RTType, RTUnit, RetentionTime,
    UnresolvedModificationPolicy, get_aa_sequence_with,
};
use crate::chemistry::{AASequence, ModificationsDB, PeptideFragmentType, SequenceModification};
use crate::format::FileType;
use crate::metadata::{CVTerm, MetaValue, MetaValueData};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, Read, Write};
use std::path::Path;

const PARAM_NAME: &str = "TransitionTSVFile";

/// How the retention-time column is interpreted
/// (`retentionTimeInterpretation`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RetentionTimeInterpretation {
    /// Normalized iRT units, stored without a unit (`iRT`, the default).
    #[default]
    Irt,
    /// Local retention time in seconds (`seconds`).
    Seconds,
    /// Local retention time in minutes (`minutes`).
    Minutes,
}

impl RetentionTimeInterpretation {
    /// The parameter spelling.
    pub fn name(self) -> &'static str {
        match self {
            Self::Irt => "iRT",
            Self::Seconds => "seconds",
            Self::Minutes => "minutes",
        }
    }
}

/// One line of a transition list (the source's protected `TSVTransition`).
///
/// Numeric fields use `-1` for "not given", as in the source.
#[derive(Clone, Debug, PartialEq)]
pub struct TsvTransition {
    /// Precursor m/z.
    pub precursor: f64,
    /// Product (fragment ion) m/z.
    pub product: f64,
    /// Normalized retention time.
    pub rt_calibrated: f64,
    /// Unique transition name.
    pub transition_name: String,
    /// Collision energy.
    pub ce: f64,
    /// Relative library intensity of the fragment ion.
    pub library_intensity: f64,
    /// Transition group identifier (transitions of the same analyte).
    pub group_id: String,
    /// Whether the transition is a decoy.
    pub decoy: bool,
    /// Peptide sequence, amino acids only.
    pub peptide_sequence: String,
    /// Protein identifiers.
    pub protein_name: Vec<String>,
    /// Gene identifier.
    pub gene_name: String,
    /// Fragment ion annotation.
    pub annotation: String,
    /// Full peptide sequence with UniMod modifications.
    pub full_peptide_name: String,
    /// Compound name (metabolomics).
    pub compound_name: String,
    /// SMILES (metabolomics).
    pub smiles: String,
    /// Molecular formula (metabolomics).
    pub sum_formula: String,
    /// Adducts (metabolomics).
    pub adducts: String,
    /// Precursor charge, as text.
    pub precursor_charge: String,
    /// Peptide group label (isotopically labelled forms of one peptide).
    pub peptide_group_label: String,
    /// Label type, e.g. `heavy` or `light`.
    pub label_type: String,
    /// Fragment charge, as text; `NA` by default.
    pub fragment_charge: String,
    /// Fragment number (7 for a y7 ion).
    pub fragment_nr: i32,
    /// Fragment m/z delta to the theoretical ion.
    pub fragment_mzdelta: f64,
    /// Ion mobility drift time.
    pub drift_time: f64,
    /// Fragment modification (a negative value is a neutral loss).
    pub fragment_modification: i32,
    /// Fragment type (`y` for a y7 ion).
    pub fragment_type: String,
    /// UniProt identifiers of the proteins.
    pub uniprot_id: Vec<String>,
    /// Use the transition to detect the peak group.
    pub detecting_transition: bool,
    /// Use the transition for peptidoform inference (IPF).
    pub identifying_transition: bool,
    /// Use the transition to quantify the peak group.
    pub quantifying_transition: bool,
    /// Peptidoforms.
    pub peptidoforms: Vec<String>,
}

impl Default for TsvTransition {
    fn default() -> Self {
        Self {
            precursor: -1.0,
            product: -1.0,
            rt_calibrated: -1.0,
            transition_name: String::new(),
            ce: -1.0,
            library_intensity: -1.0,
            group_id: String::new(),
            decoy: false,
            peptide_sequence: String::new(),
            protein_name: Vec::new(),
            gene_name: String::new(),
            annotation: String::new(),
            full_peptide_name: String::new(),
            compound_name: String::new(),
            smiles: String::new(),
            sum_formula: String::new(),
            adducts: String::new(),
            precursor_charge: String::new(),
            peptide_group_label: String::new(),
            label_type: String::new(),
            fragment_charge: "NA".into(),
            fragment_nr: -1,
            fragment_mzdelta: -1.0,
            drift_time: -1.0,
            fragment_modification: 0,
            fragment_type: String::new(),
            uniprot_id: Vec::new(),
            detecting_transition: true,
            identifying_transition: false,
            quantifying_transition: true,
            peptidoforms: Vec::new(),
        }
    }
}

impl TsvTransition {
    /// Whether the line is a peptide: by convention, when the compound name
    /// is empty or `NA`.
    pub fn is_peptide(&self) -> bool {
        self.compound_name.is_empty() || self.compound_name == "NA"
    }
}

/// The header the writer emits, in order (`header_names_`).
pub const HEADER_NAMES: [&str; 29] = [
    "PrecursorMz",
    "ProductMz",
    "PrecursorCharge",
    "ProductCharge",
    "LibraryIntensity",
    "NormalizedRetentionTime",
    "PeptideSequence",
    "ModifiedPeptideSequence",
    "PeptideGroupLabel",
    "LabelType",
    "CompoundName",
    "SumFormula",
    "SMILES",
    "Adducts",
    "ProteinId",
    "UniprotId",
    "GeneName",
    "FragmentType",
    "FragmentSeriesNumber",
    "Annotation",
    "CollisionEnergy",
    "PrecursorIonMobility",
    "TransitionGroupId",
    "TransitionId",
    "Decoy",
    "DetectingTransition",
    "IdentifyingTransition",
    "QuantifyingTransition",
    "Peptidoforms",
];

/// Resource ceilings for reading a transition list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Largest accepted input, in bytes. Default 4 GiB.
    pub max_bytes: u64,
    /// Largest number of transitions. Default 50 million.
    pub max_transitions: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 4 * 1024 * 1024 * 1024,
            max_transitions: 50_000_000,
        }
    }
}

/// The transition list reader and writer, with its three parameters.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TransitionTSVFile {
    /// How the retention-time column is interpreted.
    pub retention_time_interpretation: RetentionTimeInterpretation,
    /// Keep a peptide label group whose members have different sequences
    /// instead of splitting it (`override_group_label_check`).
    pub override_group_label_check: bool,
    /// Read a peptide whose modified sequence does not parse from its plain
    /// sequence instead of refusing it (`force_invalid_mods`).
    pub force_invalid_mods: bool,
    /// Resource ceilings for reading.
    pub limits: ReadLimits,
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

fn string_entry(name: &str, value: &str, description: &str, valid: &[&str]) -> ParamEntry {
    let mut entry = ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::String(value.into()),
        ..ParamEntry::default()
    };
    entry.tags.insert("advanced".into());
    entry.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    entry
}

/// `StringUtils::toDouble`: surrounding whitespace ignored.
fn to_double(text: &str) -> Result<f64> {
    text.trim_matches([' ', '\t', '\n', '\r'])
        .parse::<f64>()
        .map_err(|_| {
            bad(format!(
                "Could not convert string '{text}' to a double value"
            ))
        })
}

/// `StringUtils::toInt32`: surrounding whitespace ignored.
fn to_int32(text: &str) -> Result<i32> {
    text.trim_matches([' ', '\t', '\n', '\r'])
        .parse::<i32>()
        .map_err(|_| bad(format!("Could not convert string '{text}' to an int value")))
}

/// `StringUtils::split(s, c)`: empty input gives no parts.
pub(super) fn split(text: &str, separator: char) -> Vec<String> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).map(str::to_owned).collect()
    }
}

/// `StringUtils::split(s, "/")` on a multi-character-capable splitter.
fn split_str(text: &str, separator: &str) -> Vec<String> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).map(str::to_owned).collect()
    }
}

/// Remove `"`, `'` and `,`, as the reader cleans every text field.
fn clean(text: &str) -> String {
    text.chars()
        .filter(|c| !matches!(c, '"' | '\'' | ','))
        .collect()
}

/// `DataValue::toString()` of a metadata value (full precision).
pub(super) fn value_text(value: &MetaValue) -> String {
    use crate::param::value::format_float;
    fn list<T>(values: &[T], f: impl Fn(&T) -> String) -> String {
        format!("[{}]", values.iter().map(f).collect::<Vec<_>>().join(", "))
    }
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(s) => s.clone(),
        MetaValueData::Integer(i) => i.to_string(),
        MetaValueData::Float(x) => format_float(*x, true),
        MetaValueData::StringList(v) => list(v, Clone::clone),
        MetaValueData::IntegerList(v) => list(v, ToString::to_string),
        MetaValueData::FloatList(v) => list(v, |x| format_float(*x, true)),
    }
}

/// One parsed line: the header-name to column map and the fields.
struct Row<'a> {
    header: &'a BTreeMap<String, usize>,
    fields: &'a [String],
}

impl Row<'_> {
    fn raw(&self, name: &str) -> Option<&str> {
        self.header
            .get(name)
            .and_then(|&i| self.fields.get(i))
            .map(String::as_str)
    }

    /// `extractName<std::string>`: present (even empty) sets the value,
    /// cleaned of quotes and commas.
    fn text(&self, value: &mut String, names: &[&str]) -> bool {
        for name in names {
            if let Some(raw) = self.raw(name) {
                *value = clean(raw);
                return true;
            }
        }
        false
    }

    /// `extractName<double>`: an empty field counts as absent.
    fn double(&self, value: &mut f64, names: &[&str]) -> Result<bool> {
        for name in names {
            if let Some(raw) = self.raw(name).filter(|r| !r.is_empty()) {
                *value = to_double(raw)?;
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// `extractName<int>`: an empty field counts as absent.
    fn int(&self, value: &mut i32, names: &[&str]) -> Result<bool> {
        for name in names {
            if let Some(raw) = self.raw(name).filter(|r| !r.is_empty()) {
                *value = to_int32(raw)?;
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// `extractName<bool>`: `1`/`TRUE` and `0`/`FALSE` (any case); any other
    /// non-empty value counts as absent, so the next synonym is tried.
    fn flag(&self, value: &mut bool, names: &[&str]) -> bool {
        for name in names {
            if let Some(raw) = self.raw(name).filter(|r| !r.is_empty()) {
                if raw == "1" || raw.eq_ignore_ascii_case("true") {
                    *value = true;
                    return true;
                }
                if raw == "0" || raw.eq_ignore_ascii_case("false") {
                    *value = false;
                    return true;
                }
            }
        }
        false
    }
}

/// `TextFile::getLine`: lines end at `\n`, `\r` or `\r\n`; a final line
/// without a terminator is returned, an empty remainder is not.
fn read_lines(mut input: impl BufRead, max_bytes: u64) -> Result<Vec<String>> {
    let mut bytes = Vec::new();
    let read = (&mut input)
        .take(max_bytes.saturating_add(1))
        .read_to_end(&mut bytes)?;
    if u64::try_from(read).unwrap_or(u64::MAX) > max_bytes {
        return Err(bad(format!(
            "transition list exceeds the {max_bytes}-byte input limit"
        )));
    }
    let text =
        String::from_utf8(bytes).map_err(|e| bad(format!("transition list is not UTF-8: {e}")))?;
    let mut lines = Vec::new();
    let mut current = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\n' => lines.push(std::mem::take(&mut current)),
            '\r' => {
                if chars.peek() == Some(&'\n') {
                    chars.next();
                }
                lines.push(std::mem::take(&mut current));
            }
            _ => current.push(c),
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }
    Ok(lines)
}

impl TransitionTSVFile {
    /// A reader/writer with the default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters (source `getDefaults`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let flag = |b: bool| if b { "true" } else { "false" };
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                string_entry(
                    "retentionTimeInterpretation",
                    self.retention_time_interpretation.name(),
                    "How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)",
                    &["iRT", "seconds", "minutes"],
                ),
                string_entry(
                    "override_group_label_check",
                    flag(self.override_group_label_check),
                    "Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.",
                    &["true", "false"],
                ),
                string_entry(
                    "force_invalid_mods",
                    flag(self.force_invalid_mods),
                    "Force reading even if invalid modifications are encountered (OpenMS may not recognize the modification)",
                    &["true", "false"],
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build from parameters (source `setParameters` + `updateMembers_`);
    /// missing parameters take their defaults and unknown ones are ignored.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a wrong type or a value outside the valid
    /// strings (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(PARAM_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let interpretation = match param.value("retentionTimeInterpretation")?.as_str()? {
            "seconds" => RetentionTimeInterpretation::Seconds,
            "minutes" => RetentionTimeInterpretation::Minutes,
            _ => RetentionTimeInterpretation::Irt,
        };
        Ok(Self {
            retention_time_interpretation: interpretation,
            override_group_label_check: param.value("override_group_label_check")?.to_bool()?,
            force_invalid_mods: param.value("force_invalid_mods")?.to_bool()?,
            limits: ReadLimits::default(),
        })
    }

    /// Determine the separator of a header line and map each (cleaned)
    /// column name to its index (`getTSVHeader_`).
    ///
    /// Tab, semicolon and comma are tried in that order; the first that
    /// splits the header into at least four fields wins. A repeated column
    /// name maps to its last index.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no separator gives four fields (the
    /// source `Exception::IllegalArgument`, with its message).
    pub fn tsv_header(line: &str) -> Result<(char, BTreeMap<String, usize>)> {
        const MIN_HEADER_SIZE: usize = 4;
        let mut header = Vec::new();
        let mut delimiter = '\t';
        for candidate in ['\t', ';', ','] {
            delimiter = candidate;
            // `std::getline` drops a final empty field.
            header = if line.is_empty() {
                Vec::new()
            } else {
                let mut parts: Vec<String> = line.split(candidate).map(clean).collect();
                if line.ends_with(candidate) {
                    parts.pop();
                }
                parts
            };
            if header.len() >= MIN_HEADER_SIZE {
                break;
            }
        }
        if header.len() < MIN_HEADER_SIZE {
            let shown = if delimiter == '\t' { 't' } else { delimiter };
            return Err(bad(format!(
                "Determined your csv/tsv file to have delimiter '{shown}', but the parsed header has only {} fields instead of the minimal {MIN_HEADER_SIZE}. Please check your input file.",
                header.len()
            )));
        }
        let mut map = BTreeMap::new();
        for (index, name) in header.into_iter().enumerate() {
            map.insert(name, index);
        }
        Ok((delimiter, map))
    }

    /// Read the lines of a transition list (`readUnstructuredTSVInput_`).
    ///
    /// `file_type` [`FileType::Mrm`] reads a headerless SpectraST MRM list;
    /// anything else reads a list with a header. Returns the transitions and
    /// the warnings the source prints.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a header without four fields, a line whose
    /// field count differs from the number of distinct column names, a
    /// missing required column, a number that does not convert, or a
    /// sequence that does not parse (all `Exception::IllegalArgument` or
    /// `ConversionError` in the source); [`Error::Io`] from `input`.
    pub fn read_transitions(
        &self,
        input: impl BufRead,
        file_type: FileType,
    ) -> Result<(Vec<TsvTransition>, Vec<String>)> {
        let lines = read_lines(input, self.limits.max_bytes)?;
        let mut lines = lines.into_iter();
        let mrm = file_type == FileType::Mrm;
        let (delimiter, header) = if mrm {
            let names = [
                "SpectraSTBestSample",
                "SpectraSTmaxNumUsed/totNumUsed",
                "SpectraSTpI",
                "PrecursorMz",
                "SpectraSTRetentionTime",
                "ProductMz",
                "LibraryIntensity",
                "SpectraSTAnnotation",
                "FragmentCharge",
                "SpectraSTFullPeptideName",
                "SpectraSTUnknown",
                "SpectraSTNumberOfProteinsMappedTo",
                "ProteinName",
            ];
            (
                '\t',
                names
                    .iter()
                    .enumerate()
                    .map(|(i, n)| ((*n).to_owned(), i))
                    .collect::<BTreeMap<_, _>>(),
            )
        } else {
            Self::tsv_header(&lines.next().unwrap_or_default())?
        };
        let mut transitions = Vec::new();
        let mut warnings = Vec::new();
        let mut spectrast_legacy = false;
        for (index, line) in lines.enumerate() {
            let count = index + 1;
            let fields: Vec<String> = line.split(delimiter).map(str::to_owned).collect();
            if fields.len() != header.len() {
                return Err(bad(format!(
                    "Error reading the file on line {count}: length of the header and length of the line do not match: {} != {}",
                    fields.len(),
                    header.len()
                )));
            }
            if transitions.len() >= self.limits.max_transitions {
                return Err(bad("transition list exceeds the transition limit"));
            }
            let row = Row {
                header: &header,
                fields: &fields,
            };
            let mut t = TsvTransition::default();
            let mut skip = false;
            t.precursor = to_double(
                row.raw("PrecursorMz")
                    .ok_or_else(|| bad("Expected a header named PrecursorMz but found none"))?,
            )?;
            if !row.double(&mut t.product, &["ProductMz", "FragmentMz"])? {
                return Err(bad(
                    "Expected a header named ProductMz or FragmentMz but found none",
                ));
            }
            if !row.double(
                &mut t.library_intensity,
                &[
                    "LibraryIntensity",
                    "RelativeIntensity",
                    "RelativeFragmentIntensity",
                ],
            )? {
                return Err(bad(
                    "Expected a header named LibraryIntensity or RelativeFragmentIntensity but found none",
                ));
            }
            if !row.double(
                &mut t.rt_calibrated,
                &[
                    "RetentionTimeCalculatorScore",
                    "iRT",
                    "NormalizedRetentionTime",
                    "RetentionTime",
                    "Tr_recalibrated",
                ],
            )? {
                if let Some(raw) = row.raw("SpectraSTRetentionTime") {
                    spectrast_rt_extract(raw, &mut t.rt_calibrated, &mut spectrast_legacy)?;
                } else {
                    return Err(bad(
                        "Expected a header named RetentionTime, NormalizedRetentionTime, iRT, RetentionTimeCalculatorScore, Tr_recalibrated or SpectraSTRetentionTime but found none",
                    ));
                }
            }
            row.text(&mut t.precursor_charge, &["PrecursorCharge", "Charge"]);
            row.text(&mut t.fragment_type, &["FragmentType", "FragmentIonType"]);
            row.text(&mut t.fragment_charge, &["FragmentCharge", "ProductCharge"]);
            row.int(
                &mut t.fragment_nr,
                &[
                    "FragmentSeriesNumber",
                    "FragmentNumber",
                    "FragmentIonOrdinal",
                ],
            )?;
            row.double(&mut t.drift_time, &["PrecursorIonMobility"])?;
            row.double(&mut t.fragment_mzdelta, &["FragmentMzDelta"])?;
            row.int(&mut t.fragment_modification, &["FragmentModification"])?;
            row.text(&mut t.gene_name, &["GeneName"]);
            let mut proteins = String::new();
            row.text(&mut proteins, &["ProteinName", "ProteinId"]);
            if proteins != "NA" && !proteins.is_empty() {
                t.protein_name = split(&proteins, ';');
            }
            row.text(&mut t.peptide_group_label, &["PeptideGroupLabel"]);
            row.text(&mut t.label_type, &["LabelType"]);
            row.text(
                &mut t.peptide_sequence,
                &["PeptideSequence", "Sequence", "StrippedSequence"],
            );
            row.text(
                &mut t.full_peptide_name,
                &[
                    "FullUniModPeptideName",
                    "FullPeptideName",
                    "ModifiedSequence",
                    "ModifiedPeptideSequence",
                ],
            );
            row.flag(
                &mut t.detecting_transition,
                &["detecting_transition", "DetectingTransition"],
            );
            row.flag(
                &mut t.identifying_transition,
                &["identifying_transition", "IdentifyingTransition"],
            );
            row.flag(
                &mut t.quantifying_transition,
                &[
                    "quantifying_transition",
                    "QuantifyingTransition",
                    "Quantitative",
                ],
            );
            let mut peptidoforms = String::new();
            row.text(&mut peptidoforms, &["Peptidoforms"]);
            t.peptidoforms = split(&peptidoforms, '|');
            row.text(&mut t.compound_name, &["CompoundName"]);
            row.text(&mut t.sum_formula, &["SumFormula"]);
            row.text(&mut t.smiles, &["SMILES"]);
            row.text(&mut t.adducts, &["Adducts"]);
            row.text(&mut t.annotation, &["Annotation"]);
            let mut uniprot = String::new();
            row.text(&mut uniprot, &["UniprotId", "UniprotID"]);
            if uniprot != "NA" && !uniprot.is_empty() {
                t.uniprot_id = split(&uniprot, ';');
            }
            row.double(&mut t.ce, &["CE", "CollisionEnergy"])?;
            row.flag(&mut t.decoy, &["decoy", "Decoy", "IsDecoy"]);
            if let Some(raw) = row.raw("SpectraSTAnnotation") {
                skip = spectrast_annotation_extract(raw, &mut t)?;
            }
            if mrm {
                let name = row.raw("SpectraSTFullPeptideName").unwrap_or("");
                let parts = split_str(name, "/");
                let sequence = AASequence::parse(parts.first().map_or("", String::as_str))?;
                t.full_peptide_name = sequence.to_string();
                t.peptide_sequence = sequence.as_str().to_owned();
                t.precursor_charge = parts
                    .get(1)
                    .cloned()
                    .ok_or_else(|| bad(format!("SpectraST peptide name '{name}' has no charge")))?;
                t.transition_name = count.to_string();
                t.group_id = format!("{}_{}", t.full_peptide_name, t.precursor_charge);
            } else {
                if !row.text(
                    &mut t.transition_name,
                    &["transition_name", "TransitionName", "TransitionId"],
                ) {
                    t.transition_name = count.to_string();
                }
                if !row.text(
                    &mut t.group_id,
                    &[
                        "transition_group_id",
                        "TransitionGroupId",
                        "TransitionGroupName",
                    ],
                ) {
                    let sequence = AASequence::parse(&t.full_peptide_name)?;
                    t.group_id = format!("{sequence}_{}", t.precursor_charge);
                }
            }
            cleanup_transition(&mut t);
            if !skip {
                transitions.push(t);
            }
        }
        if spectrast_legacy
            && self.retention_time_interpretation == RetentionTimeInterpretation::Irt
        {
            warnings.push("Warning: SpectraST was not run in RT normalization mode but the converted list was interpreted to have iRT units. Check whether you need to adapt the parameter -algorithm:retentionTimeInterpretation. You can ignore this warning if you used a legacy SpectraST 4.0 file.".into());
        }
        Ok((transitions, warnings))
    }

    /// Read a transition list file into a [`TargetedExperiment`]
    /// (`convertTSVToTargetedExperiment`).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read (the source reads an
    /// unopenable file as empty and then fails on its header), and everything
    /// [`read_transitions`](Self::read_transitions) and
    /// [`to_targeted_experiment`](Self::to_targeted_experiment) return.
    pub fn convert_tsv_to_targeted_experiment(
        &self,
        path: impl AsRef<Path>,
        file_type: FileType,
    ) -> Result<TargetedExperiment> {
        Ok(self
            .convert_tsv_to_targeted_experiment_with_warnings(path, file_type)?
            .0)
    }

    /// As [`convert_tsv_to_targeted_experiment`](Self::convert_tsv_to_targeted_experiment),
    /// also returning the warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`convert_tsv_to_targeted_experiment`](Self::convert_tsv_to_targeted_experiment).
    pub fn convert_tsv_to_targeted_experiment_with_warnings(
        &self,
        path: impl AsRef<Path>,
        file_type: FileType,
    ) -> Result<(TargetedExperiment, Vec<String>)> {
        let file = std::io::BufReader::new(std::fs::File::open(path.as_ref())?);
        self.read_targeted_experiment(file, file_type)
    }

    /// Read a transition list from a stream into a [`TargetedExperiment`].
    ///
    /// # Errors
    ///
    /// As [`read_transitions`](Self::read_transitions) and
    /// [`to_targeted_experiment`](Self::to_targeted_experiment).
    pub fn read_targeted_experiment(
        &self,
        input: impl BufRead,
        file_type: FileType,
    ) -> Result<(TargetedExperiment, Vec<String>)> {
        let (mut transitions, mut warnings) = self.read_transitions(input, file_type)?;
        let exp = self.to_targeted_experiment(&mut transitions, &mut warnings)?;
        Ok((exp, warnings))
    }

    /// Convert read lines into a [`TargetedExperiment`] with transitions,
    /// peptides or compounds, and proteins (`TSVToTargetedExperiment_`).
    ///
    /// Mixed peptide label groups are resolved first (see
    /// [`resolve_mixed_sequence_groups`](Self::resolve_mixed_sequence_groups)).
    /// Every line becomes a transition; the first line of a transition group
    /// creates its peptide or compound, and the first line naming a protein
    /// creates it. A protein's UniProt accession is taken from the same
    /// position of the line's `UniprotId` list when both lists have the same
    /// length.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a charge or fragment number that does not
    /// convert or does not fit, and for a peptide sequence that does not parse
    /// unless [`force_invalid_mods`](Self::force_invalid_mods) is set.
    pub fn to_targeted_experiment(
        &self,
        transitions: &mut [TsvTransition],
        warnings: &mut Vec<String>,
    ) -> Result<TargetedExperiment> {
        self.resolve_mixed_sequence_groups(transitions, warnings);
        let mut exp = TargetedExperiment::new();
        let mut compounds = Vec::new();
        let mut peptides = Vec::new();
        let mut proteins = Vec::new();
        let mut peptide_ids = std::collections::BTreeSet::new();
        let mut compound_ids = std::collections::BTreeSet::new();
        let mut protein_ids = std::collections::BTreeSet::new();
        for t in transitions.iter() {
            exp.add_transition(create_transition(t)?);
            if !peptide_ids.contains(&t.group_id) && !compound_ids.contains(&t.group_id) {
                if t.is_peptide() {
                    let peptide = self.create_peptide(t, warnings)?;
                    peptide_ids.insert(peptide.id.clone());
                    peptides.push(peptide);
                } else {
                    let compound = self.create_compound(t)?;
                    compound_ids.insert(compound.id.clone());
                    compounds.push(compound);
                }
            }
            for (i, name) in t.protein_name.iter().enumerate() {
                if t.is_peptide() && !protein_ids.contains(name) {
                    let uniprot = if t.uniprot_id.len() == t.protein_name.len() {
                        t.uniprot_id[i].as_str()
                    } else {
                        ""
                    };
                    proteins.push(create_protein(name, uniprot)?);
                    protein_ids.insert(name.clone());
                }
            }
        }
        exp.set_compounds(compounds);
        exp.set_peptides(peptides);
        exp.set_proteins(proteins);
        Ok(exp)
    }

    /// Resolve peptide label groups whose members have different sequences
    /// (`resolveMixedSequenceGroups_`).
    ///
    /// Members of one label group (`MS:1000893`) should only be isotopically
    /// modified forms of the same peptide, so different sequences in a group
    /// most likely are an error. Each member whose sequence differs from the
    /// group's first is given its transition group id as label, unless
    /// [`override_group_label_check`](Self::override_group_label_check) is
    /// set, which only warns.
    pub fn resolve_mixed_sequence_groups(
        &self,
        transitions: &mut [TsvTransition],
        warnings: &mut Vec<String>,
    ) {
        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (index, t) in transitions.iter().enumerate() {
            if !t.peptide_group_label.is_empty() {
                groups
                    .entry(t.peptide_group_label.clone())
                    .or_default()
                    .push(index);
            }
        }
        for (label, members) in groups {
            let first = transitions[members[0]].peptide_sequence.clone();
            for index in members {
                let t = &mut transitions[index];
                if !first.is_empty() && t.peptide_sequence != first {
                    if self.override_group_label_check {
                        warnings.push(format!("Warning: Found multiple peptide sequences for peptide label group {label}. Since 'override_group_label_check' is on, nothing will be changed."));
                    } else {
                        warnings.push(format!("Warning: Found multiple peptide sequences for peptide label group {label}. This is most likely an error and to fix this, a new peptide label group will be inferred - to override this decision, please use the override_group_label_check parameter."));
                        t.peptide_group_label = t.group_id.clone();
                    }
                }
            }
        }
    }

    /// The retention time of a peptide or compound (`interpretRetentionTime_`).
    pub fn interpret_retention_time(&self, rt: f64) -> RetentionTime {
        let mut retention_time = RetentionTime::default();
        retention_time.set_rt(rt);
        match self.retention_time_interpretation {
            RetentionTimeInterpretation::Irt => {
                retention_time.retention_time_type = RTType::Irt;
            }
            RetentionTimeInterpretation::Seconds => {
                retention_time.retention_time_type = RTType::Local;
                retention_time.retention_time_unit = RTUnit::Second;
            }
            RetentionTimeInterpretation::Minutes => {
                retention_time.retention_time_type = RTType::Local;
                retention_time.retention_time_unit = RTUnit::Minute;
            }
        }
        retention_time
    }

    /// A peptide from a line (`createPeptide_`).
    fn create_peptide(&self, t: &TsvTransition, warnings: &mut Vec<String>) -> Result<Peptide> {
        let mut peptide = Peptide::default();
        peptide.id = t.group_id.clone();
        peptide.sequence = t.peptide_sequence.clone();
        let meta = &mut peptide.base.cv_terms.metadata;
        meta.insert(
            "full_peptide_name".into(),
            MetaValue::from(t.full_peptide_name.as_str()),
        );
        if !t.label_type.is_empty() {
            meta.insert("LabelType".into(), MetaValue::from(t.label_type.as_str()));
        }
        if !t.gene_name.is_empty() {
            meta.insert("GeneName".into(), MetaValue::from(t.gene_name.as_str()));
        }
        if !t.sum_formula.is_empty() {
            meta.insert("SumFormula".into(), MetaValue::from(t.sum_formula.as_str()));
        }
        peptide.set_peptide_group_label(t.peptide_group_label.clone());
        if !t.precursor_charge.is_empty() && t.precursor_charge != "NA" {
            peptide.set_charge_state(to_int32(&t.precursor_charge)?);
        }
        peptide.rts = vec![self.interpret_retention_time(t.rt_calibrated)];
        if t.drift_time >= 0.0 {
            peptide.set_drift_time(t.drift_time);
        }
        let sequence_text = if t.full_peptide_name.is_empty() {
            &t.peptide_sequence
        } else {
            &t.full_peptide_name
        };
        let aa_sequence = match AASequence::parse(sequence_text) {
            Ok(sequence) => sequence,
            Err(error) => {
                if self.force_invalid_mods {
                    AASequence::parse(&t.peptide_sequence)?
                } else {
                    return Err(bad(format!(
                        "Invalid input, cannot parse: {} ({error}; use 'force_invalid_mods' to override)",
                        t.full_peptide_name
                    )));
                }
            }
        };
        peptide.protein_refs = t.protein_name.clone();
        if peptide.sequence != aa_sequence.as_str() {
            if self.force_invalid_mods {
                return Ok(peptide);
            }
            warnings.push(format!("Warning: The peptide sequence {} and the full peptide name {aa_sequence} are not equal. Please check your input.", peptide.sequence));
        }
        let mut mods = Vec::new();
        if let Some(m) = aa_sequence.n_terminal_modification() {
            mods.push(modification(-1, m)?);
        }
        if let Some(m) = aa_sequence.c_terminal_modification() {
            let location = i32::try_from(aa_sequence.len())
                .map_err(|_| bad("peptide sequence too long for a modification location"))?;
            mods.push(modification(location, m)?);
        }
        for index in 0..aa_sequence.len() {
            if let Some(m) = aa_sequence.residue_modification(index)? {
                let location = i32::try_from(index)
                    .map_err(|_| bad("peptide sequence too long for a modification location"))?;
                mods.push(modification(location, m)?);
            }
        }
        peptide.mods = mods;
        Ok(peptide)
    }

    /// A compound (metabolite) from a line (`createCompound_`).
    fn create_compound(&self, t: &TsvTransition) -> Result<Compound> {
        let mut compound = Compound {
            molecular_formula: t.sum_formula.clone(),
            smiles_string: t.smiles.clone(),
            ..Compound::default()
        };
        compound.id = t.group_id.clone();
        let meta = &mut compound.base.cv_terms.metadata;
        meta.insert(
            "CompoundName".into(),
            MetaValue::from(t.compound_name.as_str()),
        );
        if !t.adducts.is_empty() {
            meta.insert("Adducts".into(), MetaValue::from(t.adducts.as_str()));
        }
        if !t.label_type.is_empty() {
            meta.insert("LabelType".into(), MetaValue::from(t.label_type.as_str()));
        }
        if t.drift_time >= 0.0 {
            compound.set_drift_time(t.drift_time);
        }
        if !t.precursor_charge.is_empty() && t.precursor_charge != "NA" {
            compound.set_charge_state(to_int32(&t.precursor_charge)?);
        }
        compound.rts = vec![self.interpret_retention_time(t.rt_calibrated)];
        Ok(compound)
    }

    /// Validate an experiment: every id unique and every reference
    /// resolvable (`validateTargetedExperiment`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when
    /// [`TargetedExperiment::contains_invalid_references`] reports a problem
    /// (the source `Exception::IllegalArgument`).
    pub fn validate_targeted_experiment(&self, exp: &TargetedExperiment) -> Result<()> {
        if let Some(problem) = exp.first_invalid_reference() {
            return Err(bad(format!(
                "Invalid input, contains duplicate or invalid references ({problem})"
            )));
        }
        Ok(())
    }

    /// Write an experiment as a transition list file
    /// (`convertTargetedExperimentToTSV`).
    ///
    /// The file is published only when the whole list was written.
    ///
    /// # Errors
    ///
    /// As [`write_targeted_experiment`](Self::write_targeted_experiment), plus
    /// [`Error::Io`] when the file cannot be written.
    pub fn convert_targeted_experiment_to_tsv(
        &self,
        path: impl AsRef<Path>,
        exp: &TargetedExperiment,
    ) -> Result<()> {
        let mut bytes = Vec::new();
        self.write_targeted_experiment(&mut bytes, exp)?;
        std::fs::write(path.as_ref(), bytes)?;
        Ok(())
    }

    /// Write an experiment as a transition list to a stream
    /// (`writeTSVOutput_` after the reference check).
    ///
    /// One line per transition under [`HEADER_NAMES`]; numbers are written as
    /// `StringUtils::toStr(double)` (full precision, e.g. `500.0`), flags as
    /// `0`/`1`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the experiment contains invalid
    /// references (the source `Exception::IllegalArgument`, "Your input file
    /// contains invalid references, cannot process file."), when a peptide's
    /// modified sequence cannot be built, or when a collision energy term's
    /// value does not convert; [`Error::Io`] from `out`.
    pub fn write_targeted_experiment(
        &self,
        mut out: impl Write,
        exp: &TargetedExperiment,
    ) -> Result<()> {
        if exp.contains_invalid_references() {
            return Err(bad(
                "Your input file contains invalid references, cannot process file.",
            ));
        }
        let rows = exp
            .transitions()
            .iter()
            .map(|t| convert_transition(t, exp))
            .collect::<Result<Vec<_>>>()?;
        use crate::param::value::format_float;
        let f = |x: f64| format_float(x, true);
        let b = |x: bool| if x { "1" } else { "0" };
        let mut text = HEADER_NAMES.join("\t");
        text.push('\n');
        for t in &rows {
            let fields = [
                f(t.precursor),
                f(t.product),
                t.precursor_charge.clone(),
                t.fragment_charge.clone(),
                f(t.library_intensity),
                f(t.rt_calibrated),
                t.peptide_sequence.clone(),
                t.full_peptide_name.clone(),
                t.peptide_group_label.clone(),
                t.label_type.clone(),
                t.compound_name.clone(),
                t.sum_formula.clone(),
                t.smiles.clone(),
                t.adducts.clone(),
                t.protein_name.join(";"),
                t.uniprot_id.join(";"),
                t.gene_name.clone(),
                t.fragment_type.clone(),
                t.fragment_nr.to_string(),
                t.annotation.clone(),
                f(t.ce),
                f(t.drift_time),
                t.group_id.clone(),
                t.transition_name.clone(),
                b(t.decoy).into(),
                b(t.detecting_transition).into(),
                b(t.identifying_transition).into(),
                b(t.quantifying_transition).into(),
                t.peptidoforms.join("|"),
            ];
            text.push_str(&fields.join("\t"));
            text.push('\n');
        }
        out.write_all(text.as_bytes())?;
        Ok(())
    }
}

/// `cleanupTransitions_`: a full peptide name like `PEPTIDE/2` becomes the
/// name and the precursor charge.
fn cleanup_transition(t: &mut TsvTransition) {
    let parts = split_str(&t.full_peptide_name, "/");
    if parts.len() == 2 {
        t.full_peptide_name = parts[0].clone();
        t.precursor_charge = parts[1].clone();
    }
}

/// `spectrastRTExtract`: `3887.50(57.30)` gives the normalized `57.30`; a
/// plain number marks a legacy (not RT-normalized) SpectraST run.
fn spectrast_rt_extract(text: &str, value: &mut f64, legacy: &mut bool) -> Result<()> {
    if let Some(start) = text.find('(') {
        if let Some(end) = text.find(')') {
            let inner = text.get(start + 1..end).unwrap_or("");
            *value = to_double(inner)?;
        }
    } else {
        *legacy = true;
        *value = to_double(text)?;
    }
    Ok(())
}

/// `spectrastAnnotationExtract`: parse the best SpectraST fragment
/// annotation (`y13^2/0.000,b16-18^2/-0.013,...`). Returns `true` when the
/// transition is to be skipped because its first annotation is not a unique
/// plain fragment (brackets, immonium `I`, precursor `p`, isotope `i`,
/// internal `m` or `?`).
fn spectrast_annotation_extract(text: &str, t: &mut TsvTransition) -> Result<bool> {
    t.annotation = text.to_owned();
    let all = split_str(text, ",");
    let first = all.first().map(String::as_str).unwrap_or("");
    if first.contains(['[', ']', 'I', 'p', 'i', 'm', '?']) {
        return Ok(true);
    }
    let with_deviation = split_str(first, "/");
    let mut best = with_deviation.first().cloned().unwrap_or_default();
    if best.contains('^') {
        let parts = split_str(&best, "^");
        t.fragment_charge = parts.get(1).cloned().unwrap_or_default();
        best = parts[0].clone();
    } else {
        t.fragment_charge = "1".into();
    }
    let head_and_number = |part: &str| -> Result<(String, i32)> {
        let kind: String = part.chars().take(1).collect();
        let number = to_int32(part.get(kind.len()..).unwrap_or(""))?;
        Ok((kind, number))
    };
    if best.contains('-') {
        let parts = split_str(&best, "-");
        let (kind, number) = head_and_number(&parts[0])?;
        t.fragment_type = kind;
        t.fragment_nr = number;
        t.fragment_modification = to_int32(parts.get(1).map_or("", String::as_str))?
            .checked_neg()
            .ok_or_else(|| bad("fragment modification out of range"))?;
    } else if best.contains('+') {
        let parts = split_str(&best, "+");
        let (kind, number) = head_and_number(&parts[0])?;
        t.fragment_type = kind;
        t.fragment_nr = number;
        t.fragment_modification = to_int32(parts.get(1).map_or("", String::as_str))?;
    } else {
        let (kind, number) = head_and_number(&best)?;
        t.fragment_type = kind;
        t.fragment_nr = number;
        t.fragment_modification = 0;
    }
    t.fragment_mzdelta = to_double(with_deviation.get(1).map_or("", String::as_str))?;
    Ok(false)
}

fn cv_term(accession: &str, name: &str, value: Option<MetaValue>) -> Result<CVTerm> {
    let mut term = CVTerm::new(accession, name, "MS");
    if let Some(value) = value {
        term.value = value;
    }
    Ok(term)
}

/// A transition from a line (`createTransition_`).
fn create_transition(t: &TsvTransition) -> Result<ReactionMonitoringTransition> {
    let mut rm = ReactionMonitoringTransition::default();
    rm.set_native_id(t.transition_name.clone());
    rm.set_precursor_mz(t.precursor);
    rm.set_product_mz(t.product);
    if t.is_peptide() {
        rm.set_peptide_ref(t.group_id.clone());
    } else {
        rm.set_compound_ref(t.group_id.clone());
    }
    rm.set_library_intensity(t.library_intensity);
    let mut product = rm.product().clone();
    if !t.fragment_charge.is_empty() && t.fragment_charge != "NA" {
        product.set_charge_state(to_int32(&t.fragment_charge)?);
    }
    let mut interpretation = Interpretation::default();
    let interpretation_set = t.fragment_nr != -1
        || t.fragment_mzdelta != -1.0
        || t.fragment_modification < 0
        || !t.fragment_type.is_empty();
    if t.fragment_nr != -1 {
        // Only the best interpretation is stored.
        interpretation.rank = 1;
        // The source assigns the int to an `unsigned char`, wrapping modulo
        // 256; a number outside 0..=255 is refused here instead.
        interpretation.ordinal = u8::try_from(t.fragment_nr).map_err(|_| {
            bad(format!(
                "fragment number {} does not fit the stored 8-bit ordinal",
                t.fragment_nr
            ))
        })?;
    }
    if t.fragment_mzdelta != -1.0 {
        interpretation.cv_terms.add(cv_term(
            "MS:1000904",
            "product ion m/z delta",
            Some(MetaValue::try_from(t.fragment_mzdelta)?),
        )?)?;
    }
    if t.fragment_modification < 0 {
        interpretation.cv_terms.add(cv_term(
            "MS:1001524",
            "fragment neutral loss",
            Some(MetaValue::from(t.fragment_modification)),
        )?)?;
    }
    let term_only = |accession: &str, name: &str, i: &mut Interpretation| -> Result<()> {
        i.cv_terms.add(cv_term(accession, name, None)?)
    };
    match t.fragment_type.as_str() {
        "v" => term_only("MS:1001237", "frag: v ion", &mut interpretation)?,
        "w" => term_only("MS:1001238", "frag: w ion", &mut interpretation)?,
        "x" => interpretation.iontype = PeptideFragmentType::XIon,
        "y" => interpretation.iontype = PeptideFragmentType::YIon,
        "z" => interpretation.iontype = PeptideFragmentType::ZIon,
        "a" => interpretation.iontype = PeptideFragmentType::AIon,
        "b" => interpretation.iontype = PeptideFragmentType::BIon,
        "c" => interpretation.iontype = PeptideFragmentType::CIon,
        "d" => term_only("MS:1001236", "frag: d ion", &mut interpretation)?,
        "unknown" => interpretation.iontype = PeptideFragmentType::NonIdentified,
        "" => interpretation.iontype = PeptideFragmentType::Unannotated,
        _ => interpretation.iontype = PeptideFragmentType::NonIdentified,
    }
    if interpretation_set {
        product.add_interpretation(interpretation);
    }
    rm.set_product(product);
    if t.ce > 0.0 {
        rm.cv_terms.add(cv_term(
            "MS:1000045",
            "collision energy",
            Some(MetaValue::try_from(t.ce)?),
        )?)?;
    }
    rm.set_decoy_transition_type(if t.decoy {
        DecoyTransitionType::Decoy
    } else {
        DecoyTransitionType::Target
    });
    if !t.annotation.is_empty() {
        rm.cv_terms
            .metadata
            .insert("annotation".into(), MetaValue::from(t.annotation.as_str()));
    }
    rm.set_detecting_transition(t.detecting_transition);
    rm.set_identifying_transition(t.identifying_transition);
    rm.set_quantifying_transition(t.quantifying_transition);
    if !t.peptidoforms.is_empty() {
        rm.cv_terms.metadata.insert(
            "Peptidoforms".into(),
            MetaValue::from(t.peptidoforms.join("|")),
        );
    }
    Ok(rm)
}

/// A protein from a line (`createProtein_`): its id and, when given, its
/// UniProt accession as `MS:1000885`.
fn create_protein(name: &str, uniprot_id: &str) -> Result<Protein> {
    let mut protein = Protein {
        id: name.to_owned(),
        ..Protein::default()
    };
    if !uniprot_id.is_empty() {
        protein.cv_terms.add(cv_term(
            "MS:1000885",
            "protein accession",
            Some(MetaValue::from(uniprot_id)),
        )?)?;
    }
    Ok(protein)
}

/// A TraML modification record from a sequence modification
/// (`addModification_`).
///
/// # Errors
///
/// [`Error::Unsupported`] for a mass tag without a registry modification
/// (`PEPTIDEK[+8.0]` style): the source creates a user-defined modification
/// whose average mass this port does not model.
fn modification(location: i32, m: &SequenceModification) -> Result<PeptideModification> {
    let (mono, avg) = match (m.diff_mono_mass(), m.diff_average_mass()) {
        (Ok(mono), Ok(avg)) => (mono, avg),
        _ => {
            return Err(Error::Unsupported(format!(
                "modification '{}' has no registry average mass; the source's user-defined \
                 modifications are not ported",
                m.name()
            )));
        }
    };
    Ok(PeptideModification {
        location,
        mono_mass_delta: mono,
        avg_mass_delta: avg,
        unimod_id: m
            .record_id()
            .and_then(|id| i32::try_from(id).ok())
            .unwrap_or(-1),
        ..PeptideModification::default()
    })
}

/// A line for the writer from a transition (`convertTransition_`).
pub(crate) fn convert_transition(
    rm: &ReactionMonitoringTransition,
    exp: &TargetedExperiment,
) -> Result<TsvTransition> {
    let mut t = TsvTransition {
        precursor: rm.precursor_mz(),
        product: rm.product_mz(),
        ..TsvTransition::default()
    };
    if !rm.peptide_ref().is_empty() {
        let peptide = exp.peptide_by_ref(rm.peptide_ref())?;
        t.group_id = rm.peptide_ref().to_owned();
        if peptide.has_retention_time() {
            t.rt_calibrated = peptide.retention_time()?;
        }
        t.peptide_sequence = peptide.sequence.clone();
        t.gene_name = "NA".into();
        for reference in &peptide.protein_refs {
            let protein = exp.protein_by_ref(reference)?;
            t.protein_name.push(protein.id.clone());
            if let Some(terms) = protein.cv_terms.get("MS:1000885") {
                if let Some(first) = terms.first() {
                    t.uniprot_id.push(value_text(&first.value));
                }
            }
        }
        // The source skips a modification it cannot resolve (with a warning
        // on stderr); the writer does the same.
        let conversion = get_aa_sequence_with(
            peptide,
            ModificationsDB::global(),
            UnresolvedModificationPolicy::SourceSkip,
        )?;
        t.full_peptide_name = conversion.sequence.to_unimod_string()?;
        if peptide.drift_time() >= 0.0 {
            t.drift_time = peptide.drift_time();
        }
        t.precursor_charge = peptide
            .charge_state()
            .map_or_else(|| "NA".to_owned(), |c| c.to_string());
        t.peptide_group_label = if peptide.peptide_group_label().is_empty() {
            "NA".into()
        } else {
            peptide.peptide_group_label().to_owned()
        };
        let meta = &peptide.cv_terms.metadata;
        if let Some(value) = meta.get("LabelType") {
            t.label_type = value_text(value);
        }
        if let Some(value) = meta.get("GeneName") {
            t.gene_name = value_text(value);
        }
    } else if !rm.compound_ref().is_empty() {
        let compound = exp.compound_by_ref(rm.compound_ref())?;
        t.group_id = rm.compound_ref().to_owned();
        if compound.has_retention_time() {
            t.rt_calibrated = compound.retention_time()?;
        }
        if compound.drift_time() >= 0.0 {
            t.drift_time = compound.drift_time();
        }
        t.precursor_charge = compound
            .charge_state()
            .map_or_else(|| "NA".to_owned(), |c| c.to_string());
        t.sum_formula = compound.molecular_formula.clone();
        t.smiles = compound.smiles_string.clone();
        let meta = &compound.cv_terms.metadata;
        if let Some(value) = meta.get("CompoundName") {
            t.compound_name = value_text(value);
        }
        if let Some(value) = meta.get("Adducts") {
            t.adducts = value_text(value);
        }
    }
    if let Some(charge) = rm.product_charge_state() {
        t.fragment_charge = charge.to_string();
    }
    let interpretations = rm.product().interpretation_list();
    for interpretation in interpretations {
        // Only the first / best interpretation is reported.
        if interpretation.rank == 1 || interpretations.len() == 1 {
            if interpretation.ordinal != 0 {
                t.fragment_nr = i32::from(interpretation.ordinal);
            }
            let kind = match interpretation.iontype {
                PeptideFragmentType::AIon => Some("a"),
                PeptideFragmentType::BIon => Some("b"),
                PeptideFragmentType::CIon => Some("c"),
                PeptideFragmentType::XIon => Some("x"),
                PeptideFragmentType::YIon => Some("y"),
                PeptideFragmentType::ZIon => Some("z"),
                PeptideFragmentType::Zp1Ion => Some("z."),
                PeptideFragmentType::Zp2Ion => Some("z'"),
                PeptideFragmentType::Precursor => Some("prec"),
                // The source spells water with a zero.
                PeptideFragmentType::BIonMinusH2O => Some("b-H20"),
                PeptideFragmentType::YIonMinusH2O => Some("y-H20"),
                PeptideFragmentType::BIonMinusNH3 => Some("b-NH3"),
                PeptideFragmentType::YIonMinusNH3 => Some("y-NH3"),
                PeptideFragmentType::NonIdentified => Some("unknown"),
                PeptideFragmentType::Unannotated => Some(""),
                PeptideFragmentType::Full
                | PeptideFragmentType::Internal
                | PeptideFragmentType::NTerminal
                | PeptideFragmentType::CTerminal => None,
            };
            if let Some(kind) = kind {
                t.fragment_type = kind.to_owned();
            }
        }
    }
    t.transition_name = rm.native_id().to_owned();
    if let Some(terms) = rm.cv_terms.get("MS:1000045") {
        if let Some(first) = terms.first() {
            t.ce = to_double(&value_text(&first.value))?;
        }
    }
    if rm.library_intensity() > -100.0 {
        t.library_intensity = rm.library_intensity();
    }
    t.decoy = rm.decoy_transition_type() == DecoyTransitionType::Decoy;
    t.annotation = "NA".into();
    if let Some(value) = rm.cv_terms.metadata.get("annotation") {
        t.annotation = value_text(value);
    }
    if let Some(value) = rm.cv_terms.metadata.get("Peptidoforms") {
        t.peptidoforms = split(&value_text(value), '|');
    }
    t.detecting_transition = rm.is_detecting_transition();
    t.identifying_transition = rm.is_identifying_transition();
    t.quantifying_transition = rm.is_quantifying_transition();
    Ok(t)
}

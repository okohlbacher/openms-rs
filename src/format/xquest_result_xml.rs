// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! xQuest result files (`.xquest.xml`), ported from
//! `FORMAT/XQuestResultXMLFile.h` and `FORMAT/HANDLERS/XQuestResultXMLHandler.h`.
//!
//! These files store the results of chemical cross-linking coupled to mass
//! spectrometry (XL-MS). They are written by xQuest and by OpenPepXL, and read
//! by `IDFileConverter`, `XFDR` and `OpenPepXL`. A documented schema is the
//! source's `share/OpenMS/SCHEMAS/xQuest_1_0.xsd`; it is not consulted, as the
//! source does not validate against it either.
//!
//! One `<xquest_results>` element carries the search settings; each
//! `<spectrum_search>` names a light/heavy spectrum pair (the same spectrum
//! twice for label-free cross-linkers) and holds its `<search_hit>`s. A hit is
//! a cross-link (`type="xlink"`, two peptides), a loop-link (`intralink`) or a
//! mono-link (`monolink`).
//!
//! # Loading
//!
//! [`XQuestResultXmlHandler::parse`](crate::format::xquest_result_xml::XQuestResultXmlHandler::parse) turns every `<search_hit>` into one
//! [`PeptideIdentification`](crate::identification::PeptideIdentification) with an alpha hit and, for a cross-link, a beta
//! hit, and collects every protein accession into the single
//! [`ProteinIdentification`](crate::identification::ProteinIdentification) (search engine `xQuest`, protocol `MS:1002494`).
//! [`XQuestResultXmlFile::load`](crate::format::xquest_result_xml::XQuestResultXmlFile::load) then applies the five `OPXLHelper` steps the
//! source's `load` applies, in its order:
//!
//! 1. `addXLTargetDecoyMV`: the alpha hit gets `xl_target_decoy_alpha` and
//!    `xl_target_decoy_beta`; a cross-link one of whose chains has no target
//!    accession becomes a decoy;
//! 2. `addBetaAccessions`: `accessions_beta` is the beta hit's accessions,
//!    comma-separated, or `-`;
//! 3. `combineTopRanksFromPairs(ids, number_of_hits)`: the identifications are
//!    regrouped by the alpha hit's `spectrum_index` (in byte order of its
//!    text), each group sorted by descending alpha score *as the Release
//!    build's `std::sort` leaves it* (not stably: ties of more than sixteen
//!    identifications are observable), consecutive equal identifications
//!    removed, and `xl_rank` renumbered from 1;
//! 4. `removeBetaPeptideHits`: each beta hit is folded into its alpha hit as
//!    the `BetaPepEv:pre`, `:post`, `:start` and `:end` meta values, and the
//!    identifications are merged by their `spectrum_reference` (in byte order
//!    of that text), one hit per cross-link;
//! 5. `computeDeltaScores`: each identification is sorted stably by descending
//!    score and every hit gets `delta_score`, the next hit's score divided by
//!    its own (0 for the last).
//!
//! # Storing
//!
//! [`XQuestResultXmlFile::store`](crate::format::xquest_result_xml::XQuestResultXmlFile::store) writes the OpenPepXL flavour of the format,
//! byte for byte as the source's `XQuestResultXMLHandler::writeTo` does
//! (numbers in the stream's 15-significant-digit form or through
//! `StringUtils::toStr`, attribute values unescaped), except that the meta
//! values written at the end of every `<search_hit>` come in byte order of
//! their names, not in the order in which the source process first registered
//! each name. `docs/XQUEST_RESULT_XML_SUPPORT.md` lists every member, the
//! preserved source conventions and every native difference.
//!
//! The spec.xml writers (`writeXQuestXMLSpec`, both overloads) are not ported:
//! they take `OPXLDataStructs::CrossLinkSpectrumMatch` and
//! `OPXLDataStructs::PreprocessedPairSpectra`, which have no port yet.

use crate::chemistry::{AASequence, EmpiricalFormula, ModificationsDB, ProteaseDB};
use crate::chemistry::{ResidueModification, SequenceModification, TermSpecificity};
use crate::comparison::Tolerance;
use crate::concept::constants::{PROTON_MASS_U, user_param};
use crate::data_structures::DateTime;
use crate::format::file_info::text_format::{list_to_string, ostream_g, to_str};
use crate::format::file_types::{FileType, has_valid_extension};
use crate::format::pepxml_mascot::permissive;
use crate::format::sax_xml::{
    self, SaxAttributes as Attributes, SaxEncoding, SaxEvent, SaxLimits, xerces_parse_int,
};
use crate::format::spectrum_meta_data_lookup::{to_double, to_int32};
use crate::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification, SearchParameters,
};
use crate::math::source_sort::source_sort_reversed_by;
use crate::metadata::{MetaInfo, MetaValue, MetaValueData};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// The schema the source's constructor names, `XMLFile("/SCHEMAS/xQuest_1_0.xsd", "1.0")`.
pub const SCHEMA_LOCATION: &str = "/SCHEMAS/xQuest_1_0.xsd";
/// The schema version the source's constructor names.
pub const SCHEMA_VERSION: &str = "1.0";

/// The source's static `XQuestResultXMLHandler::enzymes`: xQuest's
/// `enzyme_num` and the OpenMS enzyme name it maps to.
///
/// A number not listed maps to the empty name, which no enzyme has (the
/// source's `std::map::operator[]`, which also inserts that empty entry into
/// the process-wide map).
pub const ENZYMES: [(usize, &str); 11] = [
    (1, "trypsin"),
    (2, "chymotrypsin"),
    (3, "unknown_enzyme"),
    (9, "unknown_enzyme"),
    (10, "unknown_enzyme"),
    (14, "unknown_enzyme"),
    (15, "unknown_enzyme"),
    (16, "unknown_enzyme"),
    (17, "unknown_enzyme"),
    (18, "unknown_enzyme"),
    (20, "unknown_enzyme"),
];

/// The source's static `XQuestResultXMLHandler::months`: the English month
/// abbreviation of an xQuest date and its number. Any other text maps to 0,
/// which is not a valid month.
pub const MONTHS: [(&str, u32); 12] = [
    ("Jan", 1),
    ("Feb", 2),
    ("Mar", 3),
    ("Apr", 4),
    ("May", 5),
    ("Jun", 6),
    ("Jul", 7),
    ("Aug", 8),
    ("Sep", 9),
    ("Oct", 10),
    ("Nov", 11),
    ("Dec", 12),
];

/// The mass window, in Da, within which a cross-linker mass selects a
/// modification (`searchModificationsByDiffMonoMass(..., 0.01, ...)`).
pub const MODIFICATION_MASS_TOLERANCE: f64 = 0.01;

/// Resource ceilings of one read, each checked before the work it bounds.
///
/// The defaults accept every upstream fixture and any single search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum document size in bytes.
    pub max_bytes: usize,
    /// Maximum XML events, including those of ignored elements.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum attributes on one element.
    pub max_attributes: usize,
    /// Maximum `<search_hit>` elements.
    pub max_hits: usize,
    /// Maximum characters of one `seq1` or `seq2` attribute (after the
    /// source's `X` to `M(Oxidation)` substitution).
    pub max_sequence_bytes: usize,
    /// Maximum protein accessions of one `prot1` or `prot2` attribute.
    pub max_accessions_per_hit: usize,
    /// Maximum retained log lines; later ones are counted in
    /// [`XQuestDocument::suppressed_log_lines`].
    pub max_log_lines: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1 << 30,
            max_events: 100_000_000,
            max_depth: 64,
            max_attributes: 1_024,
            max_hits: 10_000_000,
            max_sequence_bytes: 100_000,
            max_accessions_per_hit: 100_000,
            max_log_lines: 10_000,
        }
    }
}

/// What one load produced: the source's `pep_ids` and `prot_ids`
/// out-parameters and what it printed.
#[derive(Clone, Debug, PartialEq)]
pub struct XQuestDocument {
    /// The source's `prot_ids`: always exactly one run.
    pub protein_identifications: Vec<ProteinIdentification>,
    /// The source's `pep_ids`.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// Lines the source prints to standard output while loading
    /// (`Set default mono-link: <modification>`), in order.
    pub log: Vec<String>,
    /// Log lines beyond [`ReadLimits::max_log_lines`].
    pub suppressed_log_lines: usize,
}

/// The source's `XQuestResultXMLFile`: loads and stores xQuest result files
/// and remembers the hit count and score range of the last load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct XQuestResultXmlFile {
    n_hits: i32,
    min_score: f64,
    max_score: f64,
}

impl Default for XQuestResultXmlFile {
    fn default() -> Self {
        Self::new()
    }
}

impl XQuestResultXmlFile {
    /// A file adapter, the source's constructor: no hits (`-1`) until a load.
    ///
    /// The source leaves the score range uninitialized until a load; it is
    /// 0 here.
    pub fn new() -> Self {
        Self {
            n_hits: -1,
            min_score: 0.0,
            max_score: 0.0,
        }
    }

    /// `getNumberOfHits`: the `<search_hit>` count of the last successful
    /// load, or -1.
    pub fn number_of_hits(&self) -> i32 {
        self.n_hits
    }

    /// `getMinScore`: the lowest `score` of the last load, starting from 0 (a
    /// file whose scores are all positive reports 0).
    pub fn min_score(&self) -> f64 {
        self.min_score
    }

    /// `getMaxScore`: the highest `score` of the last load, starting from 0.
    pub fn max_score(&self) -> f64 {
        self.max_score
    }

    /// Load an xQuest result file, as `load(filename, pep_ids, prot_ids)`.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] when the file cannot be opened or read;
    /// - [`Error::Parse`] for a document that is not well-formed XML, a
    ///   missing required attribute, an integer attribute that is not an
    ///   integer, or a sequence that does not parse (the source's
    ///   `Exception::ParseError`);
    /// - [`Error::InvalidValue`] for a number that does not convert (the
    ///   source's `Exception::ConversionError`), an enzyme or modification
    ///   the registries do not have (`Exception::ElementNotFound`), a
    ///   `spectrum` name without an odd number of underscores
    ///   (`Exception::IllegalArgument`), an unsupported cross-link `type`
    ///   (the source's bare `std::exception`), a value the metadata model
    ///   cannot hold (a non-finite number), or an exceeded [`ReadLimits`]
    ///   ceiling;
    /// - [`Error::InvalidRange`] for a link position outside its peptide
    ///   (the source's `Exception::IndexOverflow`) and for attribute text
    ///   with fewer parts than the source indexes (undefined behaviour in the
    ///   source);
    /// - [`Error::Unsupported`] for a DTD.
    ///
    /// On error this adapter keeps the state of its previous load.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<XQuestDocument> {
        self.load_with_limits(path, &ReadLimits::default())
    }

    /// [`Self::load`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn load_with_limits(
        &mut self,
        path: impl AsRef<Path>,
        limits: &ReadLimits,
    ) -> Result<XQuestDocument> {
        let path = path.as_ref();
        let file = File::open(path)?;
        self.read_with_limits(BufReader::new(file), limits)
    }

    /// Read an xQuest result document from `reader`; see [`Self::load`].
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read(&mut self, reader: impl BufRead) -> Result<XQuestDocument> {
        self.read_with_limits(reader, &ReadLimits::default())
    }

    /// [`Self::read`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read_with_limits(
        &mut self,
        reader: impl BufRead,
        limits: &ReadLimits,
    ) -> Result<XQuestDocument> {
        let parsed = XQuestResultXmlHandler::parse(reader, limits)?;
        let n_hits = parsed.n_hits;
        let mut peptides = parsed.peptide_identifications;
        add_xl_target_decoy_mv(&mut peptides)?;
        add_beta_accessions(&mut peptides);
        // `combineTopRanksFromPairs(pep_ids, n_hits_)`: the `int` member
        // converted to `Size`.
        let mut peptides = combine_top_ranks_from_pairs(peptides, n_hits as usize)?;
        let mut peptides = remove_beta_peptide_hits(std::mem::take(&mut peptides))?;
        compute_delta_scores(&mut peptides)?;
        self.n_hits = i32::try_from(n_hits).map_err(|_| limit("hit count"))?;
        self.min_score = parsed.min_score;
        self.max_score = parsed.max_score;
        Ok(XQuestDocument {
            protein_identifications: parsed.protein_identifications,
            peptide_identifications: peptides,
            log: parsed.log,
            suppressed_log_lines: parsed.suppressed_log_lines,
        })
    }

    /// Store identifications as an xQuest result file, as
    /// `store(filename, poid, peid)`, stamped with the current local time.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when `path` has a known extension other than
    ///   `.xquest.xml` (the source's `Exception::UnableToCreateFile`);
    /// - [`Error::Io`] when the file cannot be written;
    /// - otherwise as [`XQuestResultXmlHandler::write_to`].
    pub fn store(
        &self,
        path: impl AsRef<Path>,
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
    ) -> Result<()> {
        let path = path.as_ref();
        let name = path.to_string_lossy();
        if !has_valid_extension(&name, FileType::XquestXml) {
            return Err(Error::InvalidValue(format!(
                "cannot create '{name}': invalid file extension, expected 'xquest.xml'"
            )));
        }
        // Render first so that a refused record leaves no partial file.
        let mut text = Vec::new();
        XQuestResultXmlHandler::write_to(&mut text, proteins, peptides, &DateTime::now().get())?;
        let mut file = BufWriter::new(File::create(path)?);
        file.write_all(&text)?;
        file.flush()?;
        Ok(())
    }
}

/// What the handler itself produced, before the `OPXLHelper` steps of
/// [`XQuestResultXmlFile::load`].
#[derive(Clone, Debug, PartialEq)]
pub struct ParsedXQuest {
    /// One run, as the handler's constructor creates it.
    pub protein_identifications: Vec<ProteinIdentification>,
    /// One identification per `<search_hit>`, in document order.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// `getNumberOfHits`.
    pub n_hits: u32,
    /// `getMinScore`.
    pub min_score: f64,
    /// `getMaxScore`.
    pub max_score: f64,
    /// Lines the handler prints to standard output.
    pub log: Vec<String>,
    /// Log lines beyond [`ReadLimits::max_log_lines`].
    pub suppressed_log_lines: usize,
}

/// The source's `Internal::XQuestResultXMLHandler`, both its reading and its
/// writing half, and its two string helpers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct XQuestResultXmlHandler;

impl XQuestResultXmlHandler {
    /// The reading constructor plus the parse: every callback the source's
    /// `onStartElement` and `onEndElement` handle.
    ///
    /// # Errors
    ///
    /// As [`XQuestResultXmlFile::load`].
    pub fn parse(reader: impl BufRead, limits: &ReadLimits) -> Result<ParsedXQuest> {
        let sax = SaxLimits {
            max_input_bytes: u64::try_from(limits.max_bytes).unwrap_or(u64::MAX),
            max_events: u64::try_from(limits.max_events).unwrap_or(u64::MAX),
            max_depth: limits.max_depth,
            max_attributes: limits.max_attributes,
            // Character data is not read; a run cannot exceed the document.
            max_text_bytes: limits.max_bytes.max(1),
        };
        let mut handler = Reader::new(*limits);
        sax_xml::parse_with(reader, SaxEncoding::Declared, &sax, |event| match event {
            SaxEvent::Start(name, attributes) => handler.start(name, attributes),
            SaxEvent::End(name) => handler.end(name),
            SaxEvent::Characters(_) => Ok(()),
        })?;
        Ok(ParsedXQuest {
            protein_identifications: vec![handler.protein],
            peptide_identifications: handler.peptides,
            n_hits: handler.n_hits,
            min_score: handler.min_score,
            max_score: handler.max_score,
            log: handler.log,
            suppressed_log_lines: handler.suppressed_log_lines,
        })
    }

    /// `splitByNth`: `input` split at the `n`th occurrence of `separator`,
    /// that separator dropped.
    ///
    /// The source's scan starts at the *second* character (it advances before
    /// it looks), so a separator in first position is never counted. With
    /// `n` = 0 the first part is empty and the second is `input` without its
    /// first character.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `input` has fewer than `n` counted
    /// separators (the source's `std::out_of_range` from `std::string::at`),
    /// or when `n` = 0 and `input` is empty (its `Exception::IndexOverflow`).
    pub fn split_by_nth(input: &str, separator: char, n: usize) -> Result<[String; 2]> {
        let bytes = input.as_bytes();
        let mut current = 0usize;
        let mut count = 0usize;
        let mut sep = [0u8; 4];
        let sep = separator.encode_utf8(&mut sep).as_bytes();
        // The source compares `char`s; a multi-byte separator never equals one.
        let single = (sep.len() == 1).then(|| sep[0]);
        while current < bytes.len() && count < n {
            current += 1;
            let Some(byte) = bytes.get(current) else {
                return Err(Error::InvalidRange(format!(
                    "basic_string::at: __n (which is {current}) >= this->size() (which is {})",
                    bytes.len()
                )));
            };
            if Some(*byte) == single {
                count += 1;
            }
        }
        let rest = bytes.len().checked_sub(current + 1).ok_or_else(|| {
            Error::InvalidRange(format!(
                "suffix length exceeds the {}-byte input",
                bytes.len()
            ))
        })?;
        let first = &bytes[..current];
        let second = &bytes[bytes.len() - rest..];
        Ok([lossy(first), lossy(second)])
    }

    /// `splitByMiddle`: `input` split at the middle one of an odd number of
    /// `separator`s (the `n / 2 + 1`th), through [`Self::split_by_nth`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `separator` does not occur an odd number
    /// of times (the source's `Exception::IllegalArgument`); otherwise as
    /// [`Self::split_by_nth`].
    pub fn split_by_middle(input: &str, separator: char) -> Result<[String; 2]> {
        let n = input.chars().filter(|c| *c == separator).count();
        if n != 0 && n % 2 == 1 {
            Self::split_by_nth(input, separator, n / 2 + 1)
        } else {
            Err(Error::InvalidValue(
                "The separator has to occur in the input string an uneven number of times (and at least once).".into(),
            ))
        }
    }

    /// The writing constructor plus `writeTo(os)`: the OpenPepXL flavour of
    /// the format, with `date` as the `date` attribute (the source writes
    /// `DateTime::now()`, `yyyy-MM-dd hh:mm:ss`).
    ///
    /// Every hit of every identification is written as its own
    /// `<search_hit>`; a hit starts a new `<spectrum_search>` unless its
    /// `spectrum_reference` (and `spectrum_reference_heavy`, if present) is
    /// a *string* equal to the previous one, as the source's `DataValue !=
    /// std::string` comparison decides. A loaded file stores these references
    /// as integers, so it is written with one `<spectrum_search>` per hit.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] where the source throws: no protein
    /// identification (undefined behaviour in the source); a meta value it
    /// converts that is missing or of another type
    /// (`Exception::ConversionError`: `fragment:mass_tolerance_xlinks` and
    /// `cross_link:mass` on the run, `spectrum_index`, `xl_pos1`, `xl_pos2`,
    /// `xl_mass` and the four `matched_*` counts on a hit); a hit without a
    /// peptide evidence, or a string where the source converts to a number
    /// (both undefined behaviour); a sequence without a monoisotopic mass;
    /// and [`Error::Io`] when `writer` fails.
    pub fn write_to(
        writer: &mut impl Write,
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        date: &str,
    ) -> Result<()> {
        let mut os = String::new();
        write_document(&mut os, proteins, peptides, date)?;
        writer.write_all(os.as_bytes())?;
        Ok(())
    }
}

fn lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("xQuest result {what} limit exceeded"))
}

fn short(what: &str, text: &str) -> Error {
    Error::InvalidRange(format!(
        "xQuest {what} '{text}' has fewer parts than the source reads"
    ))
}

/// `StringUtils::split(s, c)` (also the `std::string` splitter overload with
/// a one-character separator): an empty subject yields nothing.
fn split<'a>(text: &'a str, separator: &str) -> Vec<&'a str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(separator).collect()
    }
}

/// `StringUtils::toInt32` for the source's `ListUtils::create<Int>`, which
/// trims first.
fn list_ints(text: &str) -> Result<Vec<i32>> {
    split(text, ",")
        .into_iter()
        .map(|part| {
            to_int32(part.trim_matches([' ', '\t', '\n', '\r']))
                .map_err(|_| Error::InvalidValue(format!("Could not convert string '{part}'")))
        })
        .collect()
}

/// `std::istringstream(text) >> bool` without `boolalpha`: a leading
/// integer, 0 false and anything else true; no integer at all is false.
fn stream_bool(text: &str) -> bool {
    let token = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let body = token.strip_prefix(['+', '-']).unwrap_or(token);
    let digits = body.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return false;
    }
    // Any digit other than 0 makes a nonzero value; an overflowing value is
    // also true.
    body.as_bytes()[..digits].iter().any(|digit| *digit != b'0')
}

fn meta_set(meta: &mut MetaInfo, key: &str, value: MetaValueData) -> Result<()> {
    meta.insert(key.to_owned(), MetaValue::new(value)?);
    Ok(())
}

fn meta_str(meta: &mut MetaInfo, key: &str, value: impl Into<String>) {
    meta.insert(key.to_owned(), MetaValue::from(value.into()));
}

fn meta_int(meta: &mut MetaInfo, key: &str, value: i64) {
    meta.insert(key.to_owned(), MetaValue::from(value));
}

fn meta_f64(meta: &mut MetaInfo, key: &str, value: f64) -> Result<()> {
    meta_set(meta, key, MetaValueData::Float(value))
}

/// Whether a meta value is the string `text`, as the source's
/// `DataValue == std::string` decides: a value of any other type is unequal.
fn is_string(value: Option<&MetaValue>, text: &str) -> bool {
    matches!(value.map(MetaValue::data), Some(MetaValueData::String(s)) if s == text)
}

/// The handler state of the reading half: the source's members.
struct Reader {
    limits: ReadLimits,
    /// `decoy_string_`, default `decoy_`.
    decoy_string: String,
    spectrum_index_light: i32,
    spectrum_index_heavy: i32,
    cross_linker_name: String,
    protein: ProteinIdentification,
    peptides: Vec<PeptideIdentification>,
    n_hits: u32,
    min_score: f64,
    max_score: f64,
    /// `is_openpepxl_`: uninitialized in the source until `<xquest_results>`.
    is_openpepxl: bool,
    accessions: BTreeSet<String>,
    charges: BTreeSet<u32>,
    /// Uninitialized in the source; the Release build was observed to read 0.
    min_precursor_charge: u32,
    max_precursor_charge: u32,
    rt_light: f64,
    rt_heavy: f64,
    mz_light: f64,
    mz_heavy: f64,
    spectrum_input_file: String,
    log: Vec<String>,
    suppressed_log_lines: usize,
}

impl Reader {
    fn new(limits: ReadLimits) -> Self {
        let mut protein = ProteinIdentification {
            search_engine: "xQuest".into(),
            search_engine_version: crate::CORE_SDK_VERSION.into(),
            ..Default::default()
        };
        // crosslinking search = MS:1002494
        meta_str(
            &mut protein.metadata,
            "SpectrumIdentificationProtocol",
            "MS:1002494",
        );
        Self {
            limits,
            decoy_string: "decoy_".into(),
            spectrum_index_light: 0,
            spectrum_index_heavy: 0,
            cross_linker_name: String::new(),
            protein,
            peptides: Vec::new(),
            n_hits: 0,
            min_score: 0.0,
            max_score: 0.0,
            is_openpepxl: false,
            accessions: BTreeSet::new(),
            charges: BTreeSet::new(),
            min_precursor_charge: 0,
            max_precursor_charge: 0,
            rt_light: 0.0,
            rt_heavy: 0.0,
            mz_light: 0.0,
            mz_heavy: 0.0,
            spectrum_input_file: String::new(),
            log: Vec::new(),
            suppressed_log_lines: 0,
        }
    }

    fn print(&mut self, line: String) {
        if self.log.len() < self.limits.max_log_lines {
            self.log.push(line);
        } else {
            self.suppressed_log_lines += 1;
        }
    }

    fn end(&mut self, tag: &str) -> Result<()> {
        if tag == "xquest_results" && !self.is_openpepxl {
            let params = &mut self.protein.search_parameters;
            params.charges = self
                .charges
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(",");
            // Min and max searched precursor charge are not written by
            // xQuest; the charges of the results are the best estimate.
            meta_int(
                &mut params.metadata,
                "precursor:min_charge",
                i64::from(self.min_precursor_charge),
            );
            meta_int(
                &mut params.metadata,
                "precursor:max_charge",
                i64::from(self.max_precursor_charge),
            );
        }
        Ok(())
    }

    fn start(&mut self, tag: &str, attributes: &Attributes) -> Result<()> {
        match tag {
            "xquest_results" => self.results(attributes),
            "spectrum_search" => self.spectrum_search(attributes),
            "search_hit" => self.search_hit(attributes),
            _ => Ok(()),
        }
    }

    /// `extractDateTime_`.
    fn date_time(&self, text: &str) -> Result<DateTime> {
        let parts = split(text, " ");
        let mut date_time = DateTime::default();
        if self.is_openpepxl {
            // Example: 2017-03-17 23:04:50
            let (Some(date), Some(time)) = (parts.first(), parts.get(1)) else {
                return Err(short("date", text));
            };
            date_time.set_date(date)?;
            date_time.set_time(time)?;
        } else {
            // Example: Fri Dec 18 12:28:42 2015, or Fri Dec  8 12:28:42 2015
            // with an empty third field for a single-digit day.
            let correction = usize::from(parts.get(2).is_some_and(|day| day.is_empty()));
            let field = |index: usize| parts.get(index).copied().ok_or_else(|| short("date", text));
            let day = to_int32(field(2 + correction)?)? as u32;
            let year = to_int32(field(4 + correction)?)? as u32;
            let month_name = field(1)?;
            let month = MONTHS
                .iter()
                .find(|(name, _)| *name == month_name)
                .map_or(0, |(_, number)| *number);
            date_time.set_date_components(month, day, year)?;
            date_time.set_time(field(3 + correction)?)?;
        }
        Ok(date_time)
    }

    fn results(&mut self, attributes: &Attributes) -> Result<()> {
        let version = attributes.required("xquest_version")?;
        self.is_openpepxl = version.contains("OpenPepXL");

        let date_time = self.date_time(attributes.required("date")?)?;
        self.protein.date_time = Some(date_time.iso_string());

        let mut params = SearchParameters::default();
        let enzyme_name = if self.is_openpepxl {
            attributes.required("enzyme_name")?.to_owned()
        } else {
            let number = xerces_parse_int(attributes.required("enzyme_num")?)?;
            ENZYMES
                .iter()
                .find(|(key, _)| usize::try_from(number).ok() == Some(*key))
                .map(|(_, name)| (*name).to_owned())
                .unwrap_or_default()
        };
        params.digestion_enzyme = enzyme(&enzyme_name)?;

        params.missed_cleavages =
            xerces_parse_int(attributes.required("missed_cleavages")?)? as u32;
        params.database = attributes.required("database")?.to_owned();
        let ms1 = attribute_double(attributes, "ms1tolerance")?;
        let ms1_unit = attributes.required(if self.is_openpepxl {
            "tolerancemeasure_ms1"
        } else {
            "tolerancemeasure"
        })?;
        params.precursor_tolerance = if ms1_unit == "ppm" {
            Tolerance::Ppm(ms1)
        } else {
            Tolerance::Absolute(ms1)
        };
        let ms2 = attribute_double(attributes, "ms2tolerance")?;
        params.fragment_tolerance = if attributes.required("tolerancemeasure_ms2")? != "Da" {
            Tolerance::Ppm(ms2)
        } else {
            Tolerance::Absolute(ms2)
        };

        // Variable modifications. Only the xQuest style `<letter>,<mass>` is
        // looked up; the source reads that mass through
        // `double(DataValue(std::string))`, which returns the bits of the
        // string's storage pointer rather than the number, so no modification
        // is ever found (observed in the Release build). The list is set,
        // empty, whenever the attribute is present and not empty.
        if let Some(text) = attributes.get("variable_mod") {
            if !text.is_empty() {
                params.variable_modifications = Vec::new();
            }
        }
        // Fixed modifications, taken as written.
        if let Some(text) = attributes.get("fixed_mod") {
            if !text.is_empty() {
                params.fixed_modifications =
                    split(text, ",").into_iter().map(str::to_owned).collect();
            }
        }

        // Without the attribute the decoy string is taken to be a prefix.
        let decoy_prefix = attributes.get("decoy_prefix").unwrap_or("1");
        if let Some(text) = attributes.get("decoy_string") {
            if !text.is_empty() {
                self.decoy_string = text.to_owned();
            }
        }
        let decoy_prefix = stream_bool(decoy_prefix);

        if let Some(text) = attributes.get("database_dc") {
            if !text.is_empty() {
                meta_str(&mut params.metadata, "input_decoys", text);
            }
        }
        meta_int(
            &mut params.metadata,
            "decoy_prefix",
            i64::from(decoy_prefix),
        );
        meta_str(
            &mut params.metadata,
            "decoy_string",
            self.decoy_string.clone(),
        );
        let xlink_tolerance = attribute_double(attributes, "xlink_ms2tolerance")?;
        meta_f64(
            &mut params.metadata,
            "fragment:mass_tolerance_xlinks",
            xlink_tolerance,
        )?;

        if let Some(text) = attributes.get("monolinkmw") {
            if !text.is_empty() {
                let masses = split(text, ",")
                    .into_iter()
                    .map(|mass| to_double(mass.trim_matches([' ', '\t', '\n', '\r'])))
                    .collect::<Result<Vec<f64>>>()?;
                if !masses.is_empty() {
                    meta_set(
                        &mut params.metadata,
                        "cross_link:mass_monolink",
                        MetaValueData::FloatList(masses),
                    )?;
                }
            }
        }

        // xQuest writes U+2212 for the minus of a zero-length cross-linker.
        let xlinker_mass = attributes.required("xlinkermw")?.replace('\u{2212}', "-");
        meta_f64(
            &mut params.metadata,
            "cross_link:mass_mass",
            to_double(&xlinker_mass)?,
        )?;
        self.cross_linker_name = attributes.required("crosslinkername")?.to_owned();
        meta_str(
            &mut params.metadata,
            "cross_link:name",
            self.cross_linker_name.clone(),
        );
        let isotope_shift = attributes.required("cp_isotopediff")?;
        if !isotope_shift.is_empty() {
            meta_f64(
                &mut params.metadata,
                "cross_link:mass_isoshift",
                to_double(isotope_shift)?,
            )?;
        }

        let nterm_linkable = stream_bool(attributes.required("ntermxlinkable")?);
        let residues = |text: &str| -> Vec<String> {
            split(text, ",").into_iter().map(str::to_owned).collect()
        };
        match attributes.get("AArequired").filter(|text| !text.is_empty()) {
            // Older xQuest versions only allowed homobifunctional linkers.
            Some(required) => {
                let mut required = required.to_owned();
                if nterm_linkable {
                    required.push_str(",N-term");
                }
                meta_set(
                    &mut params.metadata,
                    "cross_link:residue1",
                    MetaValueData::StringList(residues(&required)),
                )?;
                meta_set(
                    &mut params.metadata,
                    "cross_link:residue2",
                    MetaValueData::StringList(residues(&required)),
                )?;
            }
            None => {
                let mut first = attributes.required("AArequired1")?.to_owned();
                let mut second = attributes.required("AArequired2")?.to_owned();
                if nterm_linkable && !(first.contains("N-term") || second.contains("N-term")) {
                    first.push_str(",N-term");
                    second.push_str(",N-term");
                }
                meta_set(
                    &mut params.metadata,
                    "cross_link:residue1",
                    MetaValueData::StringList(residues(&first)),
                )?;
                meta_set(
                    &mut params.metadata,
                    "cross_link:residue2",
                    MetaValueData::StringList(residues(&second)),
                )?;
            }
        }

        if self.is_openpepxl {
            let searched = attributes.required("charges")?;
            params.charges = searched.to_owned();
            let mut charges = list_ints(searched)?;
            charges.sort_unstable();
            let (Some(min), Some(max)) = (charges.first(), charges.last()) else {
                return Err(short("charges", searched));
            };
            meta_int(
                &mut params.metadata,
                "precursor:min_charge",
                i64::from(*min),
            );
            meta_int(
                &mut params.metadata,
                "precursor:max_charge",
                i64::from(*max),
            );
            // Required, read and discarded by the source.
            attributes.required("run_path")?;
        }
        self.protein.search_parameters = params;
        Ok(())
    }

    fn spectrum_search(&mut self, attributes: &Attributes) -> Result<()> {
        let rt_text = attributes.required("rtsecscans")?;
        let rt = split(rt_text, ":");
        let (Some(light), Some(heavy)) = (rt.first(), rt.get(1)) else {
            return Err(short("rtsecscans", rt_text));
        };
        self.rt_light = to_double(light)?;
        self.rt_heavy = to_double(heavy)?;

        let mz_text = attributes.required("mzscans")?;
        if mz_text.is_empty() {
            let precursor = to_double(attributes.required("mz_precursor")?)?;
            self.mz_light = precursor;
            self.mz_heavy = precursor;
        } else {
            let mz = split(mz_text, ":");
            let (Some(light), Some(heavy)) = (mz.first(), mz.get(1)) else {
                return Err(short("mzscans", mz_text));
            };
            self.mz_light = to_double(light)?;
            self.mz_heavy = to_double(heavy)?;
        }

        // `UInt charge_precursor = attributeAsInt_(...)`.
        let charge = xerces_parse_int(attributes.required("charge_precursor")?)? as u32;
        if self.is_openpepxl {
            self.spectrum_index_light = xerces_parse_int(attributes.required("scan_index_light")?)?;
            self.spectrum_index_heavy = xerces_parse_int(attributes.required("scan_index_heavy")?)?;
            let params = &mut self.protein.search_parameters;
            if !params.metadata.contains_key("input_mzML") {
                let spectrum = attributes.required("spectrum")?;
                let file = split(spectrum, ".");
                let Some(file) = file.first() else {
                    return Err(short("spectrum", spectrum));
                };
                meta_str(&mut params.metadata, "input_mzML", format!("{file}.mzML"));
            }
        } else {
            if charge < self.min_precursor_charge {
                self.min_precursor_charge = charge;
            }
            if charge > self.max_precursor_charge {
                self.max_precursor_charge = charge;
            }
            self.charges.insert(charge);

            // For example "C_Lee_141014_CRM_dialysis_NCE20_1.06904.06904.3_C_Lee_...":
            // light and heavy (or twice light) spectrum, separated in the
            // middle, each with its file name and indices.
            let spectrum = attributes.required("spectrum")?;
            let [light, heavy] = XQuestResultXmlHandler::split_by_middle(spectrum, '_')?;
            let light = split(&light, ".");
            let heavy = split(&heavy, ".");
            let Some(file) = light.first() else {
                return Err(short("spectrum", spectrum));
            };
            self.spectrum_input_file = (*file).to_owned();
            let from_end = |parts: &[&str], back: usize| -> Result<i32> {
                let index = parts
                    .len()
                    .checked_sub(back)
                    .ok_or_else(|| short("spectrum", spectrum))?;
                to_int32(parts[index])
            };
            let last_light = light.last().copied().unwrap_or_default();
            let back = if last_light.len() > 1 { 1 } else { 2 };
            self.spectrum_index_light = from_end(&light, back)?;
            self.spectrum_index_heavy = from_end(&heavy, back)?;
        }
        Ok(())
    }

    /// `setPeptideEvidence_`: one evidence per comma-separated accession, and
    /// one protein hit per accession not seen before.
    fn evidences(&mut self, proteins: &str, hit: &mut PeptideHit) -> Result<()> {
        let accessions = split(proteins, ",");
        if accessions.len() > self.limits.max_accessions_per_hit {
            return Err(limit("accessions per hit"));
        }
        let mut evidences = Vec::with_capacity(accessions.len());
        for accession in accessions {
            if self.accessions.insert(accession.to_owned()) {
                let mut protein = ProteinHit {
                    accession: accession.to_owned(),
                    ..Default::default()
                };
                let state = if accession.contains(&self.decoy_string) {
                    "decoy"
                } else {
                    "target"
                };
                meta_str(&mut protein.metadata, "target_decoy", state);
                self.protein.hits.push(protein);
            }
            // Positions and flanking residues are not in the file.
            evidences.push(PeptideEvidence {
                protein_accession: accession.to_owned(),
                start: None,
                end: None,
                aa_before: FlankingResidue::Unknown,
                aa_after: FlankingResidue::Unknown,
            });
        }
        hit.evidences = evidences;
        Ok(())
    }

    fn sequence_text(&self, attributes: &Attributes, name: &str) -> Result<String> {
        let mut text = attributes.required(name)?.to_owned();
        if !self.is_openpepxl {
            // xQuest writes X for an oxidised methionine.
            text = text.replace('X', "M(Oxidation)");
        }
        if text.len() > self.limits.max_sequence_bytes {
            return Err(limit("sequence length"));
        }
        Ok(text)
    }

    fn search_hit(&mut self, attributes: &Attributes) -> Result<()> {
        if self.peptides.len() >= self.limits.max_hits {
            return Err(limit("hit count"));
        }
        let mut identification = PeptideIdentification::default();
        let mut alpha = PeptideHit::default();
        let mut beta = PeptideHit::default();
        let mut beta_used = false;

        let seq1 = self.sequence_text(attributes, "seq1")?;
        let link_type = attributes.required("type")?.to_owned();
        let mut alpha_sequence = from_string(&seq1)?;
        let positions = link_position(attributes)?;
        let xl_pos = positions.0 - 1;

        // U+2212 minus, as for `xlinkermw`.
        let xl_mass_text = attributes.required("xlinkermass")?.replace('\u{2212}', "-");
        let xl_mass = to_double(&xl_mass_text)?;

        // `alpha_seq[xl_pos]` checks the index (`Exception::IndexOverflow`).
        let residue = usize::try_from(xl_pos)
            .ok()
            .and_then(|index| alpha_sequence.as_str().as_bytes().get(index).copied())
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "link position {} is outside the {}-residue peptide {seq1}",
                    positions.0,
                    alpha_sequence.len()
                ))
            })?;
        let modifications = search_by_diff_mono_mass(
            ModificationsDB::global(),
            xl_mass,
            MODIFICATION_MASS_TOLERANCE,
            char::from(residue),
            Some(TermSpecificity::Anywhere),
        );

        if link_type == "monolink" {
            if let Some(first) = modifications.first() {
                let index = xl_pos as usize;
                let mut set = false;
                // Every match is set in turn, so the last one containing the
                // cross-linker name wins: the XLMOD entries come after the
                // Unimod ones and are meant to take precedence.
                for modification in &modifications {
                    if modification.full_id().contains(&self.cross_linker_name) {
                        alpha_sequence.set_modification(index, modification.full_id())?;
                        set = true;
                    }
                }
                if !set {
                    let name = first.full_id().to_owned();
                    self.print(format!("Set default mono-link: {name}"));
                    alpha_sequence.set_modification(index, &name)?;
                }
            }
        } else {
            meta_str(
                &mut alpha.metadata,
                "xl_mod",
                self.cross_linker_name.clone(),
            );
        }
        alpha.sequence = alpha_sequence;

        // `UInt charge` passed to `setCharge(Int)`.
        let charge = xerces_parse_int(attributes.required("charge")?)?;
        alpha.charge = charge;

        let light = i64::from(self.spectrum_index_light);
        let heavy = i64::from(self.spectrum_index_heavy);
        meta_int(&mut alpha.metadata, user_param::SPECTRUM_REFERENCE, light);
        meta_int(&mut alpha.metadata, "spectrum_index", light);
        meta_str(
            &mut alpha.metadata,
            "spectrum_input_file",
            self.spectrum_input_file.clone(),
        );
        let references = if light != heavy {
            meta_f64(
                &mut alpha.metadata,
                user_param::OPENPEPXL_HEAVY_SPEC_RT,
                self.rt_heavy,
            )?;
            meta_f64(
                &mut alpha.metadata,
                user_param::OPENPEPXL_HEAVY_SPEC_MZ,
                self.mz_heavy,
            )?;
            meta_int(
                &mut alpha.metadata,
                user_param::OPENPEPXL_HEAVY_SPEC_REF,
                heavy,
            );
            meta_int(&mut alpha.metadata, "spectrum_index_heavy", heavy);
            format!("{light},{heavy}")
        } else {
            light.to_string()
        };
        identification.set_spectrum_reference(references);
        meta_str(&mut alpha.metadata, "xl_chain", "MS:1002509");

        identification.mz = Some(self.mz_light);
        identification.rt = Some(self.rt_light);
        // Hard-coded in the source's MzIdentMLHandler.
        identification.score_type = "OpenPepXL:score".into();

        let prot1 = attributes.required("prot1")?.to_owned();
        let state = if prot1.contains(&self.decoy_string) {
            "decoy"
        } else {
            "target"
        };
        meta_str(&mut alpha.metadata, "target_decoy", state);

        let score = attribute_double(attributes, "score")?;
        if score < self.min_score {
            self.min_score = score;
        }
        if score > self.max_score {
            self.max_score = score;
        }
        alpha.score = score;
        let error_rel = attribute_double(attributes, "error_rel")?;
        meta_f64(
            &mut alpha.metadata,
            user_param::PRECURSOR_ERROR_PPM_USERPARAM,
            error_rel,
        )?;

        // The attributes every hit of this record shares.
        let mut shared = MetaInfo::new();
        meta_str(&mut shared, "OpenPepXL:id", attributes.required("id")?);
        meta_f64(&mut shared, "OpenPepXL:xlinkermass", xl_mass)?;
        let rank = xerces_parse_int(attributes.required("search_hit_rank")?)?;
        meta_int(&mut shared, user_param::OPENPEPXL_XL_RANK, i64::from(rank));
        meta_f64(&mut shared, "OpenPepXL:score", score)?;
        meta_str(
            &mut shared,
            "OpenPepXL:structure",
            attributes.required("structure")?,
        );
        for (attribute, key) in [
            ("wTIC", "OpenPepXL:wTIC"),
            ("TIC", "OpenPepXL:percTIC"),
            ("intsum", "OpenPepXL:intsum"),
            ("match_odds", "OpenPepXL:match-odds"),
            ("fdr", user_param::XFDR_FDR),
        ] {
            if let Some(text) = attributes.get(attribute).filter(|text| !text.is_empty()) {
                meta_f64(&mut shared, key, to_double(text)?)?;
            }
        }
        if let Some(text) = attributes.get("xprophet_f").filter(|text| !text.is_empty()) {
            meta_int(
                &mut alpha.metadata,
                "XFDR:used_for_FDR",
                i64::from(to_int32(text)?),
            );
        }
        if let Some(text) = attributes.get("fdr_type").filter(|text| !text.is_empty()) {
            meta_str(&mut alpha.metadata, "XFDR:fdr_type", text);
        }
        alpha.metadata.extend(shared.clone());

        let matched = [
            ("matched_linear_alpha", "num_of_matched_common_ions_alpha"),
            ("matched_xlink_alpha", "num_of_matched_xlink_ions_alpha"),
            ("matched_linear_beta", "num_of_matched_common_ions_beta"),
            ("matched_xlink_beta", "num_of_matched_xlink_ions_beta"),
        ];
        for (key, attribute) in matched {
            let value = xerces_parse_int(attributes.required(attribute)?)?;
            meta_int(&mut alpha.metadata, key, i64::from(value));
        }

        meta_str(&mut alpha.metadata, "prot1", prot1.clone());
        meta_str(&mut alpha.metadata, "prot2", "-");
        meta_str(
            &mut alpha.metadata,
            user_param::OPENPEPXL_BETA_ACCESSIONS,
            "-",
        );
        meta_f64(&mut alpha.metadata, user_param::OPENPEPXL_XL_MASS, xl_mass)?;
        self.evidences(&prot1, &mut alpha)?;

        match link_type.as_str() {
            "xlink" => {
                {
                    let params = &mut self.protein.search_parameters;
                    if !params.metadata.contains_key("cross_link:mass") {
                        meta_f64(&mut params.metadata, "cross_link:mass", xl_mass)?;
                    }
                }
                beta.score = score;
                meta_f64(
                    &mut beta.metadata,
                    user_param::PRECURSOR_ERROR_PPM_USERPARAM,
                    error_rel,
                )?;
                let seq2 = self.sequence_text(attributes, "seq2")?;
                beta.sequence = from_string(&seq2)?;
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_BETA_SEQUENCE,
                    seq2.clone(),
                );
                meta_str(
                    &mut beta.metadata,
                    user_param::OPENPEPXL_BETA_SEQUENCE,
                    seq2.clone(),
                );
                beta.charge = charge;
                meta_int(&mut beta.metadata, user_param::SPECTRUM_REFERENCE, light);
                if light != heavy {
                    meta_f64(
                        &mut beta.metadata,
                        user_param::OPENPEPXL_HEAVY_SPEC_RT,
                        self.rt_heavy,
                    )?;
                    meta_f64(
                        &mut beta.metadata,
                        user_param::OPENPEPXL_HEAVY_SPEC_MZ,
                        self.mz_heavy,
                    )?;
                    meta_int(
                        &mut beta.metadata,
                        user_param::OPENPEPXL_HEAVY_SPEC_REF,
                        heavy,
                    );
                    meta_int(&mut beta.metadata, "spectrum_index_heavy", heavy);
                }
                beta.metadata.extend(shared);
                for hit in [&mut alpha, &mut beta] {
                    meta_str(
                        &mut hit.metadata,
                        user_param::OPENPEPXL_XL_TYPE,
                        "cross-link",
                    );
                    meta_int(
                        &mut hit.metadata,
                        user_param::OPENPEPXL_XL_POS1,
                        positions.0 - 1,
                    );
                    meta_int(
                        &mut hit.metadata,
                        user_param::OPENPEPXL_XL_POS2,
                        positions.1 - 1,
                    );
                }

                // A link at a terminus whose residue the cross-linker does
                // not accept there is a terminal link. The source indexes the
                // attribute text, modifications included, not the residues.
                let params = &self.protein.search_parameters;
                let residue1 = string_list(params.metadata.get("cross_link:residue1"))?;
                let residue2 = string_list(params.metadata.get("cross_link:residue2"))?;
                let at = |text: &str, index: i64| -> String {
                    // `std::string::operator[]`: position size() reads '\0'.
                    let byte = usize::try_from(index)
                        .ok()
                        .and_then(|i| text.as_bytes().get(i).copied())
                        .unwrap_or(0);
                    lossy(&[byte])
                };
                let accepts =
                    |list: &[String], letter: &str| list.iter().any(|r| r.ends_with(letter));
                let mut alpha_spec = "ANYWHERE";
                let mut beta_spec = "ANYWHERE";
                if positions.0 == 1 && !accepts(&residue1, &at(&seq1, 0)) {
                    alpha_spec = "N_TERM";
                }
                if positions.1 == 1 && !accepts(&residue2, &at(&seq2, 0)) {
                    beta_spec = "N_TERM";
                }
                if positions.0 == seq1.len() as i64
                    && !accepts(&residue1, &at(&seq1, positions.0 - 1))
                {
                    alpha_spec = "C_TERM";
                }
                if positions.1 == seq2.len() as i64
                    && !accepts(&residue2, &at(&seq2, positions.1 - 1))
                {
                    beta_spec = "C_TERM";
                }
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_TERM_SPEC_ALPHA,
                    alpha_spec,
                );
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_TERM_SPEC_BETA,
                    beta_spec,
                );

                let prot2 = attributes.required("prot2")?.to_owned();
                let state = if prot2.contains(&self.decoy_string) {
                    "decoy"
                } else {
                    "target"
                };
                meta_str(&mut beta.metadata, user_param::TARGET_DECOY, state);
                meta_str(&mut beta.metadata, "xl_chain", "MS:1002510");
                for (key, attribute) in matched {
                    let value = xerces_parse_int(attributes.required(attribute)?)?;
                    meta_int(&mut beta.metadata, key, i64::from(value));
                }
                meta_str(&mut alpha.metadata, "prot2", prot2.clone());
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_BETA_ACCESSIONS,
                    prot2.clone(),
                );
                meta_str(&mut beta.metadata, "prot1", prot1);
                meta_str(
                    &mut beta.metadata,
                    user_param::OPENPEPXL_BETA_ACCESSIONS,
                    prot2.clone(),
                );
                meta_f64(&mut beta.metadata, user_param::OPENPEPXL_XL_MASS, xl_mass)?;
                meta_str(
                    &mut beta.metadata,
                    user_param::OPENPEPXL_XL_MOD,
                    self.cross_linker_name.clone(),
                );
                self.evidences(&prot2, &mut beta)?;
                beta_used = true;
            }
            "intralink" => {
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_TYPE,
                    "loop-link",
                );
                meta_int(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_POS1,
                    positions.0 - 1,
                );
                meta_int(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_POS2,
                    positions.1 - 1,
                );
                meta_f64(&mut alpha.metadata, user_param::OPENPEPXL_XL_MASS, xl_mass)?;
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_MOD,
                    self.cross_linker_name.clone(),
                );
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_BETA_SEQUENCE,
                    "-",
                );
            }
            "monolink" => {
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_TYPE,
                    "mono-link",
                );
                meta_int(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_XL_POS1,
                    positions.0 - 1,
                );
                meta_str(&mut alpha.metadata, user_param::OPENPEPXL_XL_POS2, "-");
                meta_str(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_BETA_SEQUENCE,
                    "-",
                );
            }
            other => {
                // The source logs this and throws a bare `std::exception`.
                return Err(Error::InvalidValue(format!(
                    "ERROR: Unsupported Cross-Link type: {other}"
                )));
            }
        }

        identification.hits.push(alpha);
        if beta_used {
            identification.hits.push(beta);
        }
        self.peptides.push(identification);
        self.n_hits += 1;
        Ok(())
    }
}

/// `attributeAsDouble_`: a required attribute through `StringUtils::toDouble`.
fn attribute_double(attributes: &Attributes, name: &str) -> Result<f64> {
    to_double(attributes.required(name)?)
}

/// `getLinkPosition_`: `xlinkposition`, one or two comma-separated 1-based
/// positions; the second is 0 when there is not exactly a second one.
fn link_position(attributes: &Attributes) -> Result<(i64, i64)> {
    let text = attributes.required("xlinkposition")?;
    let parts = split(text, ",");
    let first = parts.first().ok_or_else(|| short("xlinkposition", text))?;
    let first = i64::from(to_int32(first)?);
    let second = if parts.len() == 2 {
        i64::from(to_int32(parts[1])?)
    } else {
        0
    };
    Ok((first, second))
}

/// `AASequence::fromString` with its permissive default.
fn from_string(text: &str) -> Result<AASequence> {
    AASequence::parse(&permissive(text)).map_err(|error| Error::Parse {
        line: 0,
        message: format!("Cannot convert string to amino acid sequence: {error} in: {text}"),
    })
}

/// `DataValue::toStringList`, which throws for any other type.
fn string_list(value: Option<&MetaValue>) -> Result<Vec<String>> {
    match value.map(MetaValue::data) {
        Some(MetaValueData::StringList(list)) => Ok(list.clone()),
        _ => Err(Error::InvalidValue(
            "Could not convert non-StringList DataValue to StringList".into(),
        )),
    }
}

/// `ProteaseDB::getEnzyme(name)`: the source registers every enzyme under
/// its name, its lower-case name and each synonym, a later registration
/// replacing an earlier one under the same key.
fn enzyme(name: &str) -> Result<String> {
    ProteaseDB::global()
        .enzymes()
        .iter()
        .rev()
        .find(|enzyme| {
            enzyme.name() == name
                || enzyme.name().to_ascii_lowercase() == name
                || enzyme.synonyms().contains(&name)
        })
        .map(|enzyme| enzyme.name().to_owned())
        .ok_or_else(|| Error::InvalidValue(format!("the element '{name}' could not be found")))
}

/// `ModificationsDB::searchModificationsByDiffMonoMass(mods, mass,
/// max_error, residue, term_spec)` over the registry in its order.
fn search_by_diff_mono_mass(
    db: &ModificationsDB,
    mass: f64,
    max_error: f64,
    residue: char,
    term: Option<TermSpecificity>,
) -> Vec<&ResidueModification> {
    db.entries()
        .iter()
        .map(AsRef::as_ref)
        .filter(|modification: &&ResidueModification| {
            let diff = modification.diff_mono_mass();
            (diff - mass).abs() <= max_error
                && (diff != 0.0 || mass == 0.0)
                && residues_match(residue, modification)
                && term.is_none_or(|term| term == modification.term_specificity())
        })
        .collect()
}

/// `ModificationsDB::residuesMatch_`. The global registry holds no
/// user-defined modifications, so origin `X` matches every residue.
fn residues_match(residue: char, modification: &ResidueModification) -> bool {
    match modification.origin() {
        Some(origin) if origin != 'X' => origin == residue || matches!(residue, 'X' | '.' | '?'),
        _ => true,
    }
}

// ---------------------------------------------------------------------------
// `OPXLHelper` steps of `XQuestResultXMLFile::load`, as the source runs them.
// `ANALYSIS/XLMS/OPXLHelper.h` has no port of its own yet; these follow it
// statement by statement (the mzIdentML reader's copies diverge on purpose).
// ---------------------------------------------------------------------------

/// `OPXLHelper::addXLTargetDecoyMV`.
fn add_xl_target_decoy_mv(peptides: &mut [PeptideIdentification]) -> Result<()> {
    for identification in peptides {
        if identification.hits.is_empty() {
            continue;
        }
        let beta_state = identification
            .hits
            .get(1)
            .filter(|_| identification.hits.len() == 2)
            .map(|hit| hit.metadata.get(user_param::TARGET_DECOY).cloned());
        let alpha = &mut identification.hits[0];
        let alpha_state = alpha.metadata.get(user_param::TARGET_DECOY).cloned();
        set_or_empty(
            &mut alpha.metadata,
            user_param::OPENPEPXL_TARGET_DECOY_ALPHA,
            alpha_state.clone(),
        );
        match beta_state {
            Some(beta_state) => {
                set_or_empty(
                    &mut alpha.metadata,
                    user_param::OPENPEPXL_TARGET_DECOY_BETA,
                    beta_state.clone(),
                );
                // A cross-link one of whose accession lists holds only decoys
                // is a decoy.
                if !lenient_text(alpha_state.as_ref())?.contains("target")
                    || !lenient_text(beta_state.as_ref())?.contains("target")
                {
                    meta_str(&mut alpha.metadata, user_param::TARGET_DECOY, "decoy");
                }
            }
            None => meta_str(
                &mut alpha.metadata,
                user_param::OPENPEPXL_TARGET_DECOY_BETA,
                "-",
            ),
        }
    }
    Ok(())
}

/// `setMetaValue(key, getMetaValue(other))`: a missing value is the empty
/// `DataValue`, which `setMetaValue` stores as such; the port's metadata
/// model has an empty value too.
fn set_or_empty(meta: &mut MetaInfo, key: &str, value: Option<MetaValue>) {
    meta.insert(key.to_owned(), value.unwrap_or_default());
}

/// `StringUtils::toStr(const DataValue&)`: `DataValue::toString(true)`.
fn lenient_text(value: Option<&MetaValue>) -> Result<String> {
    let Some(value) = value else {
        return Ok(String::new());
    };
    Ok(match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => to_str(*number),
        MetaValueData::StringList(list) => list_to_string(list)?,
        MetaValueData::IntegerList(list) => list_to_string(list)?,
        MetaValueData::FloatList(list) => list_to_string(list)?,
    })
}

/// `OPXLHelper::addBetaAccessions`.
fn add_beta_accessions(peptides: &mut [PeptideIdentification]) {
    for identification in peptides {
        if identification.hits.is_empty() {
            continue;
        }
        if identification.hits.len() == 2 {
            let accessions = identification.hits[1]
                .evidences
                .iter()
                .map(|evidence| evidence.protein_accession.as_str())
                .collect::<Vec<_>>()
                .join(",");
            for hit in identification.hits.iter_mut() {
                meta_str(
                    &mut hit.metadata,
                    user_param::OPENPEPXL_BETA_ACCESSIONS,
                    accessions.clone(),
                );
            }
        } else {
            meta_str(
                &mut identification.hits[0].metadata,
                user_param::OPENPEPXL_BETA_ACCESSIONS,
                "-",
            );
        }
    }
}

/// `OPXLHelper::combineTopRanksFromPairs(peptide_ids, number_top_hits)`.
///
/// The source scans every identification once per distinct index; grouping
/// once yields the same groups in the same order.
fn combine_top_ranks_from_pairs(
    peptides: Vec<PeptideIdentification>,
    number_top_hits: usize,
) -> Result<Vec<PeptideIdentification>> {
    let mut groups: BTreeMap<String, Vec<PeptideIdentification>> = BTreeMap::new();
    for identification in peptides {
        let Some(first) = identification.hits.first() else {
            continue;
        };
        let index = lenient_text(first.metadata.get("spectrum_index"))?;
        groups.entry(index).or_default().push(identification);
    }
    let mut combined = Vec::new();
    for (_, mut group) in groups {
        // `std::sort(rbegin, rend, PeptideIDScoreComparator())`.
        source_sort_reversed_by(&mut group, |a, b| match (a.hits.first(), b.hits.first()) {
            (Some(a), Some(b)) => a.score < b.score,
            _ => false,
        })?;
        group.dedup_by(|later, earlier| source_equal(earlier, later));
        group.truncate(number_top_hits);
        for (rank, identification) in group.iter_mut().enumerate() {
            let rank = i64::try_from(rank + 1).map_err(|_| limit("rank"))?;
            let two = identification.hits.len() == 2;
            for (index, hit) in identification.hits.iter_mut().enumerate() {
                if index == 0 || (index == 1 && two) {
                    meta_int(&mut hit.metadata, user_param::OPENPEPXL_XL_RANK, rank);
                }
            }
        }
        combined.extend(group);
    }
    Ok(combined)
}

/// `OPXLHelper::removeBetaPeptideHits`.
fn remove_beta_peptide_hits(
    mut peptides: Vec<PeptideIdentification>,
) -> Result<Vec<PeptideIdentification>> {
    for identification in &mut peptides {
        let folds = identification.hits.len() == 2
            && is_string(
                identification.hits[1].metadata.get("xl_chain"),
                "MS:1002510",
            );
        if !folds {
            continue;
        }
        let beta = &identification.hits[1];
        let mut pre = String::new();
        let mut post = String::new();
        let mut start = String::new();
        let mut end = String::new();
        for evidence in &beta.evidences {
            if !pre.is_empty() {
                pre.push(',');
                post.push(',');
                start.push(',');
                end.push(',');
            }
            pre.push(evidence.aa_before.code());
            post.push(evidence.aa_after.code());
            // The source's numeric `std::string += Int` overload writes the
            // decimal text; an unknown position is -1.
            start.push_str(&position_text(evidence.start));
            end.push_str(&position_text(evidence.end));
        }
        let mut hit = identification.hits[0].clone();
        hit.metadata.remove("xl_chain");
        meta_str(&mut hit.metadata, user_param::OPENPEPXL_BETA_PEPEV_PRE, pre);
        meta_str(
            &mut hit.metadata,
            user_param::OPENPEPXL_BETA_PEPEV_POST,
            post,
        );
        meta_str(
            &mut hit.metadata,
            user_param::OPENPEPXL_BETA_PEPEV_START,
            start,
        );
        meta_str(&mut hit.metadata, user_param::OPENPEPXL_BETA_PEPEV_END, end);
        identification.hits = vec![hit];
    }

    // Collect the hits of one spectrum under one identification; only the
    // first hit of each identification is kept, as in the source.
    let mut merged: BTreeMap<String, PeptideIdentification> = BTreeMap::new();
    for identification in peptides {
        let Some(mut hit) = identification.hits.first().cloned() else {
            continue;
        };
        let reference = lenient_text(identification.metadata.get(user_param::SPECTRUM_REFERENCE))?;
        let entry = merged.entry(reference.clone()).or_insert_with(|| {
            let mut fresh = PeptideIdentification {
                rt: identification.rt,
                mz: identification.mz,
                score_type: user_param::OPENPEPXL_SCORE.into(),
                ..Default::default()
            };
            fresh.set_spectrum_reference(reference);
            fresh
        });
        hit.metadata.remove("xl_chain");
        entry.hits.push(hit);
    }
    Ok(merged.into_values().collect())
}

fn position_text(position: Option<usize>) -> String {
    position.map_or_else(|| "-1".to_owned(), |value| value.to_string())
}

/// `OPXLHelper::computeDeltaScores`.
///
/// A ratio that is not finite (a score of 0 followed by a nonzero one) is
/// left unset: the port's metadata holds finite numbers only; the source
/// stores the infinity or NaN.
fn compute_delta_scores(peptides: &mut [PeptideIdentification]) -> Result<()> {
    for identification in peptides {
        // `PeptideIdentification::sort`, a `std::stable_sort` by score.
        identification.sort()?;
        let scores: Vec<f64> = identification.hits.iter().map(|hit| hit.score).collect();
        let last = scores.len().saturating_sub(1);
        for (rank, hit) in identification.hits.iter_mut().enumerate() {
            let delta = if rank == last {
                0.0
            } else {
                scores[rank + 1] / scores[rank]
            };
            if delta.is_finite() {
                meta_f64(&mut hit.metadata, user_param::DELTA_SCORE, delta)?;
            }
        }
    }
    Ok(())
}

/// `PeptideIdentification::operator==`, with `DataValue`'s equality for
/// meta values: equal types, and doubles within 1e-6 of each other.
fn source_equal(a: &PeptideIdentification, b: &PeptideIdentification) -> bool {
    meta_equal(&a.metadata, &b.metadata)
        && a.identifier == b.identifier
        && a.hits.len() == b.hits.len()
        && a.hits.iter().zip(&b.hits).all(|(a, b)| hit_equal(a, b))
        && a.significance_threshold == b.significance_threshold
        && a.score_type == b.score_type
        && a.higher_score_better == b.higher_score_better
        && a.mz == b.mz
        && a.rt == b.rt
}

/// `PeptideHit::operator==`.
fn hit_equal(a: &PeptideHit, b: &PeptideHit) -> bool {
    meta_equal(&a.metadata, &b.metadata)
        && a.sequence == b.sequence
        && a.score == b.score
        && a.rank == b.rank
        && a.charge == b.charge
        && a.evidences == b.evidences
        && a.peak_annotations == b.peak_annotations
        && a.analysis_results == b.analysis_results
}

fn meta_equal(a: &MetaInfo, b: &MetaInfo) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|((key_a, a), (key_b, b))| {
            key_a == key_b
                && a.unit() == b.unit()
                && match (a.data(), b.data()) {
                    (MetaValueData::Float(x), MetaValueData::Float(y)) => (x - y).abs() < 1e-6,
                    (x, y) => x == y,
                }
        })
}

// ---------------------------------------------------------------------------
// Writing: `XQuestResultXMLHandler::writeTo`.
// ---------------------------------------------------------------------------

/// `os << double` after the source's `os.precision(writtenDigits<double>())`.
fn stream(value: f64) -> String {
    ostream_g(value, 15)
}

fn conversion(key: &str, target: &str) -> Error {
    Error::InvalidValue(format!("Could not convert DataValue '{key}' to {target}"))
}

/// `double(DataValue)`: an integer converts; an empty value throws; a string
/// would read the storage bits (undefined), refused here.
fn value_double(meta: &MetaInfo, key: &str) -> Result<f64> {
    match meta.get(key).map(MetaValue::data) {
        Some(MetaValueData::Float(value)) => Ok(*value),
        Some(MetaValueData::Integer(value)) => Ok(*value as f64),
        _ => Err(conversion(key, "double")),
    }
}

/// `int(DataValue)`: integers only, truncated to 32 bits.
fn value_int(meta: &MetaInfo, key: &str) -> Result<i32> {
    match meta.get(key).map(MetaValue::data) {
        Some(MetaValueData::Integer(value)) => Ok(*value as i32),
        _ => Err(conversion(key, "int")),
    }
}

/// `DataValue::toString()` of a meta value, empty when it is missing.
fn meta_text(meta: &MetaInfo, key: &str) -> Result<String> {
    lenient_text(meta.get(key))
}

/// `StringUtils::substr(text, pos, n)`, clamping the start.
fn substr(text: &str, pos: usize, n: usize) -> String {
    let bytes = text.as_bytes();
    let begin = pos.min(bytes.len());
    let end = begin.saturating_add(n).min(bytes.len());
    lossy(&bytes[begin..end])
}

/// `substr(1)` then `chop(1)`: the brackets of a list's text removed.
fn strip_brackets(text: &str) -> String {
    let inner = substr(text, 1, usize::MAX);
    let keep = inner.len().saturating_sub(1);
    lossy(&inner.as_bytes()[..keep])
}

/// The keys `writeTo` writes under their own attribute names or does not
/// write, removed before the remaining meta values are written.
const CONSUMED: [&str; 32] = [
    user_param::OPENPEPXL_XL_MASS,
    user_param::OPENPEPXL_XL_RANK,
    user_param::OPENPEPXL_XL_POS1,
    user_param::OPENPEPXL_XL_TYPE,
    user_param::OPENPEPXL_XL_POS2,
    "matched_xlink_alpha",
    "matched_linear_alpha",
    "matched_xlink_beta",
    "matched_linear_beta",
    "OpenPepXL:xcorr xlink",
    "OpenPepXL:xcorr common",
    "OpenPepXL:match-odds",
    "OpenPepXL:prescore",
    "OpenPepXL:TIC",
    "OpenPepXL:wTIC",
    "OpenPepXL:intsum",
    "spectrum_reference",
    user_param::OPENPEPXL_HEAVY_SPEC_REF,
    "spectrum_index",
    "spectrum_index_heavy",
    user_param::OPENPEPXL_HEAVY_SPEC_RT,
    user_param::OPENPEPXL_HEAVY_SPEC_MZ,
    user_param::PRECURSOR_ERROR_PPM_USERPARAM,
    user_param::XFDR_FDR,
    user_param::OPENPEPXL_XL_MOD,
    "prot1",
    "prot2",
    "OpenPepXL:id",
    "OpenPepXL:percTIC",
    "OpenPepXL:score",
    "OpenPepXL:structure",
    "OpenPepXL:xlinkermass",
];

fn write_document(
    os: &mut String,
    proteins: &[ProteinIdentification],
    peptides: &[PeptideIdentification],
    date: &str,
) -> Result<()> {
    use std::fmt::Write as _;
    let run = proteins.first().ok_or_else(|| {
        Error::InvalidValue("an xQuest result file needs a protein identification".into())
    })?;
    let params = &run.search_parameters;
    let meta = &params.metadata;
    let input_filename = if meta.contains_key("input_mzML") {
        meta_text(meta, "input_mzML")?
    } else {
        String::new()
    };
    let spec_xml_name = meta_text(meta, "out_xquest_specxml")?;

    os.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    os.push_str("<?xml-stylesheet type=\"text/xsl\" href=\"\"?>\n");

    let mono_masses = strip_brackets(&meta_text(meta, "cross_link:mass_monolink")?);
    let (precursor_unit, precursor_tolerance) = unit(params.precursor_tolerance);
    let (fragment_unit, fragment_tolerance) = unit(params.fragment_tolerance);
    let xlink_tolerance = value_double(meta, "fragment:mass_tolerance_xlinks")?;
    let cross_link_name = meta_text(meta, "cross_link:name")?;
    let cross_link_mass = value_double(meta, "cross_link:mass")?;
    let isotope_shift = if meta.contains_key("cross_link:mass_isoshift") {
        value_double(meta, "cross_link:mass_isoshift")?
    } else {
        0.0
    };
    let required1 = strip_brackets(&meta_text(meta, "cross_link:residue1")?);
    let required2 = strip_brackets(&meta_text(meta, "cross_link:residue2")?);
    let nterm_linkable = required1.contains("N-term") || required2.contains("N-term");
    let decoy_fasta = meta_text(meta, "input_decoys")?;
    let variable = params.variable_modifications.join(",");
    let fixed = params.fixed_modifications.join(",");
    let decoy_prefix = meta_text(meta, "decoy_prefix")?;
    let decoy_string = meta_text(meta, "decoy_string")?;
    let runs = run.primary_ms_run_paths.join(",");

    let _ = writeln!(
        os,
        "<xquest_results xquest_version=\"OpenPepXL 1.0\" date=\"{date}\" author=\"Eugen Netz\" \
         tolerancemeasure_ms1=\"{precursor_unit}\" tolerancemeasure_ms2=\"{fragment_unit}\" \
         ms1tolerance=\"{}\" ms2tolerance=\"{}\" xlink_ms2tolerance=\"{}\" crosslinkername=\"{cross_link_name}\" \
         xlinkermw=\"{}\" monolinkmw=\"{mono_masses}\" database=\"{}\" database_dc=\"{decoy_fasta}\" \
         xlinktypes=\"1111\" AArequired1=\"{required1}\" AArequired2=\"{required2}\" cp_isotopediff=\"{}\" \
         enzyme_name=\"{}\" outputpath=\"{spec_xml_name}\" missed_cleavages=\"{}\" ntermxlinkable=\"{}\" \
         CID_match2ndisotope=\"1\" variable_mod=\"{variable}\" fixed_mod=\"{fixed}\" decoy_prefix=\"{decoy_prefix}\" \
         decoy_string=\"{decoy_string}\" charges=\"{}\" run_path=\"{runs}\" nocutatxlink=\"1\">",
        stream(precursor_tolerance),
        stream(fragment_tolerance),
        stream(xlink_tolerance),
        stream(cross_link_mass),
        params.database,
        stream(isotope_shift),
        // The source holds the enzyme object and writes its name.
        enzyme(&params.digestion_enzyme).unwrap_or_else(|_| params.digestion_enzyme.clone()),
        // `int missed_cleavages = search_params.missed_cleavages`.
        params.missed_cleavages as i32,
        u8::from(nterm_linkable),
        params.charges,
    );

    let mut current_light = String::new();
    let mut current_heavy = String::new();
    for identification in peptides {
        for hit in &identification.hits {
            write_hit(
                os,
                identification,
                hit,
                &input_filename,
                cross_link_mass,
                &mut current_light,
                &mut current_heavy,
            )?;
        }
    }
    os.push_str("</spectrum_search>\n");
    os.push_str("</xquest_results>\n");
    Ok(())
}

/// `AASequence::getMonoWeight()` (`Residue::Full`, charge 0) as the source
/// computes it: terminal delta masses, then each residue's internal weight
/// (its `Residue` mono weight minus water), then water. A modified residue's
/// weight is the formula mass of the free residue plus the delta formula, or
/// the declared absolute mass, or the free residue plus the delta mass, as
/// `Residue::setModification` sets it. The crate's
/// [`AASequence::mono_mass_for`] adds the delta formula's mass to the internal
/// residue's instead, which rounds differently in the last bits, and the
/// written `error` magnifies that.
fn source_mono_weight(sequence: &AASequence) -> Result<f64> {
    if sequence.is_empty() {
        return Ok(0.0);
    }
    let water = EmpiricalFormula::parse("H2O")?;
    let water_mass = water.mono_mass();
    let mut weight = 0.0;
    for modification in [
        sequence.n_terminal_modification(),
        sequence.c_terminal_modification(),
    ]
    .into_iter()
    .flatten()
    {
        weight += modification.diff_mono_mass()?;
    }
    let mut free: BTreeMap<char, EmpiricalFormula> = BTreeMap::new();
    for (index, residue) in sequence.as_str().chars().enumerate() {
        if residue == 'X' {
            return Err(Error::InvalidValue(
                "Cannot get weight of sequence with unknown AA 'X' with unknown mass.".into(),
            ));
        }
        let formula = match free.entry(residue) {
            std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
            // The free amino acid's formula: a one-residue peptide's.
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(AASequence::parse(&residue.to_string())?.formula()?)
            }
        };
        let mut mono = formula.mono_mass();
        match sequence.residue_modification(index)? {
            None => {}
            Some(SequenceModification::Known(modification)) => {
                let changes = modification.diff_mono_mass() != 0.0
                    || modification.diff_average_mass() != 0.0
                    || !modification.diff_formula().is_empty();
                if changes {
                    if modification.mono_mass() != 0.0 {
                        mono = modification.mono_mass();
                    } else if modification.diff_mono_mass() != 0.0 {
                        mono += modification.diff_mono_mass();
                    }
                    if !modification.diff_formula().is_empty() {
                        mono = formula
                            .checked_add(modification.diff_formula())?
                            .mono_mass();
                    } else if let Some(absolute) = modification.absolute_formula() {
                        mono = absolute.mono_mass();
                    }
                }
            }
            Some(_) => {
                return Err(Error::Unsupported(
                    "the xQuest writer needs named residue modifications".into(),
                ));
            }
        }
        weight += mono - water_mass;
    }
    Ok(weight + water_mass)
}

fn unit(tolerance: Tolerance) -> (&'static str, f64) {
    match tolerance {
        Tolerance::Ppm(value) => ("ppm", value),
        Tolerance::Absolute(value) => ("Da", value),
    }
}

#[allow(clippy::too_many_arguments)]
fn write_hit(
    os: &mut String,
    identification: &PeptideIdentification,
    hit: &PeptideHit,
    input_filename: &str,
    cross_link_mass: f64,
    current_light: &mut String,
    current_heavy: &mut String,
) -> Result<()> {
    use std::fmt::Write as _;
    let meta = &hit.metadata;
    // An unset m/z is the source's NaN.
    let precursor_mz = identification.mz.unwrap_or(f64::NAN);
    let charge = hit.charge;
    let z = f64::from(charge);
    let precursor_mass = precursor_mz * z - z * PROTON_MASS_U;

    let new_spectrum = !is_string(meta.get(user_param::SPECTRUM_REFERENCE), current_light)
        || (meta.contains_key(user_param::OPENPEPXL_HEAVY_SPEC_REF)
            && !is_string(
                meta.get(user_param::OPENPEPXL_HEAVY_SPEC_REF),
                current_heavy,
            ));
    if new_spectrum {
        if !current_light.is_empty() {
            os.push_str("</spectrum_search>\n");
        }
        *current_light = meta_text(meta, user_param::SPECTRUM_REFERENCE)?;
        current_heavy.clear();
        if meta.contains_key(user_param::OPENPEPXL_HEAVY_SPEC_REF) {
            *current_heavy = meta_text(meta, user_param::OPENPEPXL_HEAVY_SPEC_REF)?;
        }

        let base_name = if !input_filename.is_empty() {
            let file = input_filename.split('/').next_back().unwrap_or_default();
            let parts = split(file, ".");
            parts
                .first()
                .map(|part| (*part).to_owned())
                .ok_or_else(|| short("input_mzML", input_filename))?
        } else if meta.contains_key("spectrum_input_file") {
            meta_text(meta, "spectrum_input_file")?
        } else {
            String::new()
        };

        // `(Size)getMetaValue("spectrum_index")`: a non-negative integer.
        let light = match meta.get("spectrum_index").map(MetaValue::data) {
            Some(MetaValueData::Integer(value)) if *value >= 0 => *value as u64,
            _ => return Err(conversion("spectrum_index", "unsigned long")),
        };
        // `(Size)(Int)getMetaValue("spectrum_index_heavy")`.
        let heavy = if meta.contains_key("spectrum_index_heavy") {
            i64::from(value_int(meta, "spectrum_index_heavy")?) as u64
        } else {
            light
        };
        let light_name = format!("{base_name}.light.{light}");
        let heavy_name = format!("{base_name}.heavy.{heavy}");
        let rt = identification.rt.unwrap_or(f64::NAN);
        let mut rt_scans = format!("{}:", to_str(rt));
        let mut mz_scans = format!("{}:", to_str(precursor_mz));
        let scan_type = if light == heavy {
            rt_scans.push_str(&to_str(rt));
            mz_scans.push_str(&to_str(precursor_mz));
            "light"
        } else {
            rt_scans.push_str(&meta_text(meta, user_param::OPENPEPXL_HEAVY_SPEC_RT)?);
            mz_scans.push_str(&meta_text(meta, user_param::OPENPEPXL_HEAVY_SPEC_MZ)?);
            "light_heavy"
        };
        let _ = writeln!(
            os,
            "<spectrum_search spectrum=\"{light_name}_{heavy_name}\" mz_precursor=\"{}\" scantype=\"{scan_type}\" \
             charge_precursor=\"{charge}\" Mr_precursor=\"{}\" rtsecscans=\"{rt_scans}\" mzscans=\"{mz_scans}\" \
             scan_index_light=\"{light}\" scan_index_heavy=\"{heavy}\" >",
            stream(precursor_mz),
            stream(precursor_mass),
        );
    }

    let link_type = meta_text(meta, user_param::OPENPEPXL_XL_TYPE)?;
    let mut structure = hit.sequence.as_str().to_owned();
    let pos1 = value_int(meta, user_param::OPENPEPXL_XL_POS1)?;
    // `substr(structure, Int, 1)`: a negative position converts to a huge
    // `size_t`, which is clamped to the end.
    let letter_first = substr(&structure, usize::try_from(pos1).unwrap_or(usize::MAX), 1);
    let mut weight = source_mono_weight(&hit.sequence)?;
    let alpha_pos = i64::from(pos1) + 1;
    let mut beta_pos = 0i64;
    let mut topology = format!("a{alpha_pos}");
    let id;
    let mut seq_beta = String::new();
    let xltype = match link_type.as_str() {
        "cross-link" => {
            beta_pos = i64::from(value_int(meta, user_param::OPENPEPXL_XL_POS2)?) + 1;
            let beta_text = match meta
                .get(user_param::OPENPEPXL_BETA_SEQUENCE)
                .map(MetaValue::data)
            {
                Some(MetaValueData::String(text)) => text.clone(),
                _ => return Err(conversion(user_param::OPENPEPXL_BETA_SEQUENCE, "string")),
            };
            let beta = from_string(&beta_text)?;
            structure.push('-');
            structure.push_str(beta.as_str());
            topology.push_str(&format!("-b{beta_pos}"));
            weight +=
                source_mono_weight(&beta)? + value_double(meta, user_param::OPENPEPXL_XL_MASS)?;
            id = format!("{structure}-{topology}");
            seq_beta = meta_text(meta, user_param::OPENPEPXL_BETA_SEQUENCE)?;
            "xlink"
        }
        "loop-link" => {
            beta_pos = i64::from(value_int(meta, user_param::OPENPEPXL_XL_POS2)?) + 1;
            topology.push_str(&format!("-b{beta_pos}"));
            let letter_second = substr(
                &structure,
                usize::try_from(beta_pos - 1).unwrap_or(usize::MAX),
                1,
            );
            id = format!("{structure}-{letter_first}{alpha_pos}-{letter_second}{beta_pos}");
            weight += cross_link_mass;
            "intralink"
        }
        _ => {
            if meta_text(meta, user_param::OPENPEPXL_XL_MOD)?.starts_with("unknown") {
                weight += value_double(meta, user_param::OPENPEPXL_XL_MASS)?;
            }
            // `Int(double(...))`: truncation toward zero.
            let mass = value_double(meta, user_param::OPENPEPXL_XL_MASS)? as i32;
            id = format!("{structure}-{letter_first}{alpha_pos}-{mass}");
            "monolink"
        }
    };

    // The relative error is taken from the meta value, for consistency; the
    // formula applies without it.
    let theo_mz = (weight + z * PROTON_MASS_U) / z;
    let error = precursor_mz - theo_mz;
    let rel_error = if meta.contains_key(user_param::PRECURSOR_ERROR_PPM_USERPARAM) {
        value_double(meta, user_param::PRECURSOR_ERROR_PPM_USERPARAM)?
    } else {
        (error / theo_mz) / 1e-6
    };

    let prot_alpha = if hit.evidences.is_empty() {
        return Err(Error::InvalidValue(
            "an xQuest search hit needs a peptide evidence".into(),
        ));
    } else {
        hit.evidences
            .iter()
            .map(|evidence| evidence.protein_accession.as_str())
            .collect::<Vec<_>>()
            .join(",")
    };
    let prot_beta = if meta.contains_key(user_param::OPENPEPXL_BETA_ACCESSIONS)
        && !is_string(meta.get(user_param::OPENPEPXL_BETA_ACCESSIONS), "-")
    {
        meta_text(meta, user_param::OPENPEPXL_BETA_ACCESSIONS)?
    } else {
        String::new()
    };
    let mut position = alpha_pos.to_string();
    if beta_pos > 0 {
        position.push_str(&format!(",{beta_pos}"));
    }

    let xlink_alpha = value_int(meta, "matched_xlink_alpha")?;
    let xlink_beta = value_int(meta, "matched_xlink_beta")?;
    let linear_alpha = value_int(meta, "matched_linear_alpha")?;
    let linear_beta = value_int(meta, "matched_linear_beta")?;
    let _ = write!(
        os,
        "<search_hit search_hit_rank=\"{}\" id=\"{id}\" type=\"{xltype}\" structure=\"{structure}\" seq1=\"{}\" \
         seq2=\"{seq_beta}\" prot1=\"{prot_alpha}\" prot2=\"{prot_beta}\" topology=\"{topology}\" \
         xlinkposition=\"{position}\" Mr=\"{}\" mz=\"{}\" charge=\"{charge}\" xlinkermass=\"{}\" \
         measured_mass=\"{}\" error=\"{}\" error_rel=\"{}\" xlinkions_matched=\"{}\" backboneions_matched=\"{}\" \
         xcorrx=\"{}\" xcorrb=\"{}\" match_odds=\"{}\" prescore=\"{}\" num_of_matched_ions_alpha=\"{}\" \
         num_of_matched_ions_beta=\"{}\" num_of_matched_common_ions_alpha=\"{}\" num_of_matched_common_ions_beta=\"{}\" \
         num_of_matched_xlink_ions_alpha=\"{}\" num_of_matched_xlink_ions_beta=\"{}\" TIC=\"{}\" wTIC=\"{}\" intsum=\"{}",
        meta_text(meta, user_param::OPENPEPXL_XL_RANK)?,
        hit.sequence,
        stream(weight),
        stream(theo_mz),
        meta_text(meta, user_param::OPENPEPXL_XL_MASS)?,
        stream(precursor_mass),
        stream(error),
        stream(rel_error),
        xlink_alpha.wrapping_add(xlink_beta),
        linear_alpha.wrapping_add(linear_beta),
        meta_text(meta, "OpenPepXL:xcorr xlink")?,
        meta_text(meta, "OpenPepXL:xcorr common")?,
        meta_text(meta, "OpenPepXL:match-odds")?,
        meta_text(meta, "OpenPepXL:prescore")?,
        linear_alpha.wrapping_add(xlink_alpha),
        xlink_beta.wrapping_add(linear_beta),
        meta_text(meta, "matched_linear_alpha")?,
        meta_text(meta, "matched_linear_beta")?,
        meta_text(meta, "matched_xlink_alpha")?,
        meta_text(meta, "matched_xlink_beta")?,
        meta_text(meta, "OpenPepXL:TIC")?,
        meta_text(meta, "OpenPepXL:wTIC")?,
        meta_text(meta, "OpenPepXL:intsum")?,
    );
    for (key, attribute) in [
        ("XFDR:used_for_FDR", "xprophet_f"),
        ("XFDR:fdr_type", "fdr_type"),
        (user_param::XFDR_FDR, "fdr"),
    ] {
        if meta.contains_key(key) {
            let _ = write!(os, "\" {attribute}=\"{}", meta_text(meta, key)?);
        }
    }
    // Every other meta value, under its own name.
    for (key, value) in meta {
        if CONSUMED.contains(&key.as_str()) {
            continue;
        }
        let _ = write!(os, "\" {key}=\"{}", lenient_text(Some(value))?);
    }
    let _ = writeln!(
        os,
        "\" annotated_spec=\"\" score=\"{}\" >\n</search_hit>",
        stream(hit.score)
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_bool_follows_num_get() {
        assert!(stream_bool("1"));
        assert!(!stream_bool("0"));
        assert!(stream_bool(" 2"));
        assert!(stream_bool("-1"));
        assert!(!stream_bool(""));
        assert!(!stream_bool("abc"));
        assert!(!stream_bool("00"));
        assert!(stream_bool("1abc"));
    }

    #[test]
    fn strip_brackets_matches_substr_and_chop() {
        assert_eq!(strip_brackets("[K, N-term]"), "K, N-term");
        assert_eq!(strip_brackets(""), "");
        assert_eq!(strip_brackets("["), "");
    }
}

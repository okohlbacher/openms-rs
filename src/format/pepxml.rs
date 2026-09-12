// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded pepXML (Trans-Proteomic Pipeline search results) reader and writer.
//!
//! Ports `FORMAT/PepXMLFile.h` and its implementation. A pepXML document holds
//! one `<msms_pipeline_analysis>` root with one `<msms_run_summary>` per searched
//! spectra file; each run declares its enzyme (`<sample_enzyme>`), one or more
//! `<search_summary>` blocks with the search database and the modification
//! declarations, and then one `<spectrum_query>` per spectrum holding ranked
//! `<search_hit>` entries with scores, protein references and modification
//! positions.
//!
//! A documented schema for the format ships with the TPP and is also mirrored in
//! the OpenMS `share/OpenMS/SCHEMAS` directory. This port neither reads nor
//! validates against that schema; see `docs/PEPXML_SUPPORT.md` for the API
//! mapping, the preserved source conventions, the native differences and the
//! checked boundaries.
//!
//! Reading returns
//! [`PepXmlDocument`](crate::format::pepxml::PepXmlDocument), which carries the
//! run-level [`ProteinIdentification`](crate::identification::ProteinIdentification)
//! records, the spectrum-level
//! [`PeptideIdentification`](crate::identification::PeptideIdentification)
//! records and the non-fatal diagnostics the source writes to its log stream.
//!
//! ```no_run
//! use openms::format::pepxml;
//! let mut options = pepxml::ReadOptions::default();
//! options.experiment_name = "PepXMLFile_test".into();
//! let document = pepxml::load_with_options("run.pepxml", &options)?;
//! assert!(!document.peptide_identifications.is_empty());
//! # Ok::<(), openms::Error>(())
//! ```

use crate::chemistry::{
    AASequence, EmpiricalFormula, ModificationsDB, ProteaseDB, ResidueModification,
    SequenceModification, TermSpecificity,
};
use crate::comparison::Tolerance;
use crate::constants::C13C12_MASSDIFF_U;
use crate::data_structures::DateTime;
use crate::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeakMassType, PeptideEvidence, PeptideHit,
    PeptideIdentification, ProteinHit, ProteinIdentification, SearchParameters, TargetDecoyType,
};
use crate::kernel::MSSpectrum;
use crate::metadata::MetaValue;
use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::path::Path;
use std::sync::Arc;

/// Mass window used when a declared modification mass is matched against the
/// registry or against the run header, from `PepXMLFile::mod_tol_`.
pub const MODIFICATION_TOLERANCE: f64 = 0.002;

/// Terminal modifications whose absolute `massdiff` is below this are dropped.
///
/// Source `PepXMLFile::xtandem_artificial_mod_tol_`: some X!Tandem versions
/// annotate a spurious very small fixed terminal modification (possibly the
/// electron mass) which interferes with the real ones.
pub const XTANDEM_ARTIFICIAL_MODIFICATION_TOLERANCE: f64 = 0.0005;

fn bad(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}

// ---------------------------------------------------------------------------
// Options and document
// ---------------------------------------------------------------------------

/// Resource ceilings, run selection and modification preferences for reading.
///
/// Every ceiling is checked before the corresponding allocation, so exceeding
/// one leaves no partially built document. The defaults accept every upstream
/// fixture and any pepXML a search engine writes for a single LC-MS/MS run.
#[derive(Clone, Debug)]
pub struct ReadOptions {
    /// Maximum decoded input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum number of XML elements in the document.
    pub max_elements: usize,
    /// Maximum XML nesting depth.
    pub max_depth: usize,
    /// Maximum `<spectrum_query>`/`<search_result>` records.
    pub max_identifications: usize,
    /// Maximum `<search_hit>` records per `<search_result>`.
    pub max_hits: usize,
    /// Maximum modification entries per header or per hit.
    pub max_modifications: usize,
    /// Maximum `<search_hit>`/`<alternative_protein>` protein references per run.
    pub max_protein_hits: usize,
    /// Maximum retained diagnostics; further ones are counted, not stored.
    pub max_warnings: usize,
    /// Name of the spectra file whose results are wanted, extension optional.
    ///
    /// Empty accepts every run, as the source's defaulted `experiment_name`.
    /// A nonempty name that no run's `base_name` ends with is an error, matching
    /// the source's `Exception::ParseError` "Found no experiment with name".
    pub experiment_name: String,
    /// Retain the pepXML `spectrum` attribute as `pepxml_spectrum_name`.
    pub keep_native_spectrum_name: bool,
    /// Also record unrecognised `<search_score>` entries as meta values.
    pub parse_unknown_scores: bool,
    /// Fixed modifications preferred over a registry mass search.
    pub preferred_fixed_modifications: Vec<Arc<ResidueModification>>,
    /// Variable modifications preferred over a registry mass search.
    pub preferred_variable_modifications: Vec<Arc<ResidueModification>>,
    /// Spectrum metadata used when a query carries no `retention_time_sec`.
    pub lookup: SpectrumIndex,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 512 * 1024 * 1024,
            max_elements: 20_000_000,
            max_depth: 64,
            max_identifications: 2_000_000,
            max_hits: 10_000,
            max_modifications: 100_000,
            max_protein_hits: 5_000_000,
            max_warnings: 1_000,
            experiment_name: String::new(),
            keep_native_spectrum_name: false,
            parse_unknown_scores: false,
            preferred_fixed_modifications: Vec::new(),
            preferred_variable_modifications: Vec::new(),
            lookup: SpectrumIndex::new(),
        }
    }
}

/// Output ceilings and the run-level attributes the source `store` takes as
/// arguments.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    /// Spectra file: its stem becomes `base_name`, its type `raw_data`.
    ///
    /// Source `store` additionally *loads* this file to build its retention-time
    /// lookup; this port never reads spectra, so supply [`WriteOptions::lookup`]
    /// when scan numbers must come from the spectra rather than from the
    /// identification order.
    pub mz_file: String,
    /// Overrides the `base_name` attribute, as the source `mz_name` argument.
    pub mz_name: String,
    /// Fallback base name when `mz_file` and `mz_name` are both empty.
    ///
    /// [`store`] fills this from the output path, as the source does; [`write()`]
    /// has no path and uses this value verbatim.
    pub output_name: String,
    /// Emit `<analysis_timestamp>` and a PeptideProphet analysis result.
    pub peptideprophet_analyzed: bool,
    /// Retain a `pepxml_spectrum_name` meta value as the `spectrum` attribute.
    pub keep_native_spectrum_name: bool,
    /// Spectrum metadata used to resolve scan numbers and indices.
    pub lookup: SpectrumIndex,
    /// Maximum serialized output size in bytes.
    pub max_output_bytes: usize,
    /// Maximum identifications written.
    pub max_identifications: usize,
}
impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            mz_file: String::new(),
            mz_name: String::new(),
            output_name: String::new(),
            peptideprophet_analyzed: false,
            keep_native_spectrum_name: false,
            lookup: SpectrumIndex::new(),
            max_output_bytes: 512 * 1024 * 1024,
            max_identifications: 2_000_000,
        }
    }
}

/// Run-level and spectrum-level identifications read from one pepXML document.
///
/// `protein_identifications` and `peptide_identifications` are linked by
/// [`ProteinIdentification::identifier`], exactly as the source's two output
/// containers are.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PepXmlDocument {
    /// One entry per identification run, in document order.
    pub protein_identifications: Vec<ProteinIdentification>,
    /// One entry per accepted `<search_result>`, in document order.
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// Non-fatal diagnostics; the source writes these to its log stream.
    pub warnings: Vec<String>,
    /// Diagnostics suppressed once `max_warnings` was reached.
    pub suppressed_warnings: usize,
}

// ---------------------------------------------------------------------------
// Spectrum lookup
// ---------------------------------------------------------------------------

/// One spectrum's identifying metadata, as `SpectrumMetaDataLookup` records it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumMetaData {
    /// mzML `id` attribute of the spectrum.
    pub native_id: String,
    /// Retention time in seconds.
    pub rt: f64,
    /// MS level; the reader accepts a retention time only from level 2.
    pub ms_level: u32,
    /// One-based scan number, when one can be extracted from the native ID.
    pub scan_number: Option<u64>,
}

/// Native stand-in for the parts of `METADATA/SpectrumMetaDataLookup.h` and its
/// `SpectrumLookup` base that pepXML needs; both are separate unported headers.
///
/// The name is deliberately not `SpectrumLookup`: this is not a port of that
/// class, only the four queries pepXML performs - by native ID, by scan number,
/// by retention time and by index. Entry order is the spectrum order of the
/// source experiment, so an entry's position is its spectrum index.
///
/// Scan numbers are extracted with the source's `default_scan_regexp`,
/// `=(?<SCAN>\d+)$`: the digits of a trailing `=`-separated token. The source's
/// caller-registered reference formats (`addReferenceFormat`, `findByReference`)
/// and its index-based lookup are not provided.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumIndex {
    entries: Vec<SpectrumMetaData>,
    rt_tolerance: f64,
}
impl SpectrumIndex {
    /// Empty lookup with the source's 0.01 s retention-time tolerance.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            rt_tolerance: 0.01,
        }
    }
    /// Build a lookup from spectra, as `SpectrumMetaDataLookup::readSpectra`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a retention time is not finite.
    pub fn read_spectra(spectra: &[MSSpectrum]) -> Result<Self> {
        let mut result = Self::new();
        for spectrum in spectra {
            if !spectrum.rt.is_finite() {
                return Err(Error::InvalidValue(
                    "spectrum retention time must be finite".into(),
                ));
            }
            result.entries.push(SpectrumMetaData {
                native_id: spectrum.native_id.clone(),
                rt: spectrum.rt,
                ms_level: spectrum.ms_level,
                scan_number: extract_scan_number(&spectrum.native_id),
            });
        }
        Ok(result)
    }
    /// Build a lookup from already-assembled entries, in spectrum order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a retention time is not finite.
    pub fn from_entries(entries: Vec<SpectrumMetaData>) -> Result<Self> {
        if entries.iter().any(|e| !e.rt.is_finite()) {
            return Err(Error::InvalidValue(
                "spectrum retention time must be finite".into(),
            ));
        }
        Ok(Self {
            entries,
            rt_tolerance: 0.01,
        })
    }
    /// Retention-time window used by [`SpectrumIndex::find_by_rt`], in seconds.
    pub fn rt_tolerance(&self) -> f64 {
        self.rt_tolerance
    }
    /// Set the retention-time window, as the source's public `rt_tolerance`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] unless the value is finite and
    /// nonnegative. The source performs no such check.
    pub fn set_rt_tolerance(&mut self, value: f64) -> Result<()> {
        if !value.is_finite() || value < 0.0 {
            return Err(Error::InvalidValue(
                "retention-time tolerance must be finite and nonnegative".into(),
            ));
        }
        self.rt_tolerance = value;
        Ok(())
    }
    /// Number of known spectra.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether no spectra are known, as `SpectrumLookup::empty`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Metadata of one spectrum index, as `getSpectrumMetaData`.
    pub fn get(&self, index: usize) -> Option<&SpectrumMetaData> {
        self.entries.get(index)
    }
    /// Index of the spectrum with this native ID.
    ///
    /// Source `findByNativeID` throws `Exception::ElementNotFound`; this returns
    /// `None` so absence is not an error path.
    pub fn find_by_native_id(&self, native_id: &str) -> Option<usize> {
        self.entries.iter().position(|e| e.native_id == native_id)
    }
    /// Index of the spectrum with this one-based scan number.
    pub fn find_by_scan_number(&self, scan_number: u64) -> Option<usize> {
        self.entries
            .iter()
            .position(|e| e.scan_number == Some(scan_number))
    }
    /// Index of the spectrum nearest `rt` within the tolerance.
    ///
    /// Ties between an equally distant earlier and later spectrum choose the
    /// later one, because the source compares the lower difference with `<`.
    pub fn find_by_rt(&self, rt: f64) -> Option<usize> {
        if !rt.is_finite() {
            return None;
        }
        let mut best: Option<(f64, usize)> = None;
        for (index, entry) in self.entries.iter().enumerate() {
            let difference = (entry.rt - rt).abs();
            if difference > self.rt_tolerance {
                continue;
            }
            // A non-strict comparison keeps the LAST equally distant entry,
            // which is what the source produces: it keys a std::map on the
            // retention time, so a repeated time overwrites the earlier index,
            // and it returns the upper neighbour when both are equally far.
            let better = best.is_none_or(|(current, _)| difference <= current);
            if better {
                best = Some((difference, index));
            }
        }
        best.map(|(_, index)| index)
    }
    /// Whether an identifier looks like an mzML native ID.
    ///
    /// Reproduces `SpectrumNativeIDParser::isNativeID`, a prefix test against
    /// `scan=`, `scanId=`, `scanID=`, `controllerType=`, `function=`, `sample=`,
    /// `index=`, `spectrum=`, `file=` and `frame=`.
    pub fn is_native_id(id: &str) -> bool {
        [
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
        ]
        .iter()
        .any(|prefix| id.starts_with(prefix))
    }
}

/// Scan number from a native ID, per the source `default_scan_regexp`.
///
/// `=(?<SCAN>\d+)$` is anchored at the end, so it reads the digits of the LAST
/// `=`-separated token, not the `scan=` one. For a Bruker TDF identifier such
/// as `frame=5 scan=7 precursor=3` that yields 3; the source's
/// `getRegExFromNativeID` would target `scan=`, but `readSpectra` uses the
/// anchored default. The port keeps the source behaviour.
fn extract_scan_number(native_id: &str) -> Option<u64> {
    let digits = native_id.len()
        - native_id
            .trim_end_matches(|c: char| c.is_ascii_digit())
            .len();
    if digits == 0 {
        return None;
    }
    let head = &native_id[..native_id.len() - digits];
    if !head.ends_with('=') {
        return None;
    }
    native_id[native_id.len() - digits..].parse().ok()
}

// ---------------------------------------------------------------------------
// Bounded work accounting
// ---------------------------------------------------------------------------

struct Meter {
    work: usize,
    bytes: usize,
}
impl Meter {
    fn new(options: &ReadOptions) -> Self {
        Self {
            work: options.max_input_bytes.saturating_mul(32).max(1 << 24),
            bytes: options.max_input_bytes.saturating_mul(16).max(1 << 24),
        }
    }
    fn limit(&self) -> Error {
        bad("pepXML resource limit exceeded")
    }
    fn spend(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.work = self.work.checked_sub(work).ok_or_else(|| self.limit())?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(|| self.limit())?;
        Ok(())
    }
    fn cap(&self, value: usize, maximum: usize) -> Result<()> {
        if value > maximum {
            Err(self.limit())
        } else {
            Ok(())
        }
    }
}

// ---------------------------------------------------------------------------
// Header modifications (source PepXMLFile::AminoAcidModification)
// ---------------------------------------------------------------------------

/// A modification resolved from a pepXML header declaration.
///
/// Source `AminoAcidModification` always ends up holding a
/// `const ResidueModification*`, registering a `MASS_ONLY` record in the global
/// `ModificationsDB` when no known modification explains the declared mass.
/// This port never mutates a registry: an unexplained mass becomes an anonymous
/// annotation whose identifier is spelled exactly as the source's
/// `createUnknownFromMassString` spells it, so `M[+1.0]`, `.n[+2.5]` and
/// `.c[+3.4]` are unchanged.
#[derive(Clone, Debug, PartialEq)]
enum ResolvedModification {
    /// A registry entry, shared and immutable.
    Known(Arc<ResidueModification>),
    /// An anonymous monoisotopic mass annotation.
    Mass {
        /// Signed decimal spelling of the mass difference, e.g. `+2.5`.
        text: String,
        /// Identifier, e.g. `M[+1.0]` or `.n[+2.5]`.
        full_id: String,
        /// Positional specificity the declaration implies.
        term: TermSpecificity,
    },
}
impl ResolvedModification {
    /// Identifier the source reports through `AminoAcidModification::getDescription`.
    fn full_id(&self) -> &str {
        match self {
            Self::Known(value) => value.full_id(),
            Self::Mass { full_id, .. } => full_id,
        }
    }
    /// Positional specificity.
    fn term_specificity(&self) -> TermSpecificity {
        match self {
            Self::Known(value) => value.term_specificity(),
            Self::Mass { term, .. } => *term,
        }
    }
    fn is_n_terminal(&self) -> bool {
        matches!(
            self.term_specificity(),
            TermSpecificity::NTerm | TermSpecificity::ProteinNTerm
        )
    }
    fn is_c_terminal(&self) -> bool {
        matches!(
            self.term_specificity(),
            TermSpecificity::CTerm | TermSpecificity::ProteinCTerm
        )
    }
}

/// One `<aminoacid_modification>` or `<terminal_modification>` declaration.
///
/// The source uses both elements ambiguously through one class, so this type
/// does too. `mass` is the absolute mass of the modified residue or terminus and
/// `mass_diff` the difference it adds; `<mod_aminoacid_mass>` positions in the
/// hits reference the header by `mass`, not by name.
#[derive(Clone, Debug, PartialEq)]
struct HeaderModification {
    /// `aminoacid` attribute; empty for an unrestricted terminal modification.
    amino_acid: String,
    /// `massdiff` attribute.
    mass_diff: f64,
    /// `mass` attribute, corrected when the file sets it equal to `massdiff`.
    mass: f64,
    /// Whether `variable` was `Y`, case-insensitively.
    variable: bool,
    /// `description` attribute, when present.
    description: String,
    /// Lower-cased `terminus`/`peptide_terminus` attribute: `""`, `n`, `c`, `nc`.
    terminus: String,
    /// Whether the declaration names a protein rather than a peptide terminus.
    protein_terminus: bool,
    /// Positional specificity, `None` when the declaration left it open.
    term: Option<TermSpecificity>,
    /// The modification this declaration resolved to, if any.
    resolved: Option<ResolvedModification>,
    /// Diagnostics the source collects in `AminoAcidModification::errors_`.
    warnings: Vec<String>,
}

/// Signed decimal spelling of a mass difference.
///
/// Reproduces `ResidueModification::getDiffMonoMassString`: a `+`/`-` sign and
/// then the magnitude through the source's full-precision double formatter.
fn diff_mono_mass_string(value: f64) -> String {
    let sign = if value < 0.0 { "-" } else { "+" };
    format!(
        "{sign}{}",
        crate::param::value::format_float(value.abs(), true)
    )
}

/// Internal-residue monoisotopic mass of one residue code.
///
/// Source `ResidueDB::getResidue(code)->getMonoWeight(Residue::Internal)`.
/// Built from public chemistry only: the single-residue peptide's neutral mass
/// less water. Ambiguous codes (B/Z/X) have no monoisotopic composition here,
/// where the source supplies an averaged placeholder.
fn internal_residue_mass(code: char) -> Result<f64> {
    let water = EmpiricalFormula::parse("H2O")?.mono_mass();
    Ok(AASequence::parse(&code.to_string())?.mono_mass()? - water)
}

/// Monoisotopic mass added to an internal residue to make an N- or C-terminus.
///
/// Source `Residue::getInternalToNTerm()` (a hydrogen) and
/// `Residue::getInternalToCTerm()` (a hydroxyl).
fn terminal_gain(n_terminal: bool) -> Result<f64> {
    Ok(EmpiricalFormula::parse(if n_terminal { "H" } else { "OH" })?.mono_mass())
}

impl HeaderModification {
    /// Resolve one header declaration, as `AminoAcidModification`'s constructor.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when neither `aminoacid` nor a
    /// terminus is given, exactly as the source constructor throws, and
    /// [`Error::Parse`] when `massdiff` or `mass` is not a finite number.
    ///
    /// Diagnostics that the source records without throwing are collected in
    /// [`HeaderModification::warnings`].
    #[allow(clippy::too_many_arguments)]
    fn new(
        amino_acid: &str,
        mass_diff: &str,
        mass: &str,
        variable: &str,
        description: &str,
        terminus: &str,
        protein_terminus: &str,
        options: &ReadOptions,
        registry: &ModificationsDB,
    ) -> Result<Self> {
        if amino_acid.is_empty() && terminus.is_empty() {
            return Err(Error::MissingInformation(
                "either terminus or amino acid origin or both needs to be set".into(),
            ));
        }
        let mut value = Self {
            amino_acid: amino_acid.to_owned(),
            mass_diff: finite(mass_diff, "modification massdiff")?,
            mass: finite(mass, "modification mass")?,
            variable: variable.eq_ignore_ascii_case("y"),
            description: description.to_owned(),
            terminus: terminus.to_ascii_lowercase(),
            protein_terminus: false,
            term: None,
            resolved: None,
            warnings: Vec::new(),
        };
        if value.terminus == "nc" {
            value.warnings.push(
                "value 'nc' for aminoacid terminus not supported. The modification will be parsed \
                 as an unrestricted modification."
                    .into(),
            );
        }
        if value.amino_acid.chars().count() > 1 {
            value.warnings.push(
                "single modification specified for multiple amino acids, which is not supported; \
                 proceeding with the first"
                    .into(),
            );
        }
        // Source note: the schema allows only "", "n" or "c" for protein_terminus,
        // but many tools write "Y"/"N" for yes/no, which collides with "n". Upper
        // and lower case are therefore both interpreted, case-sensitively.
        let lowered = protein_terminus.to_ascii_lowercase();
        if lowered == "y" {
            value.protein_terminus = true;
        } else if lowered == "c" {
            value.protein_terminus = true;
            value.terminus = lowered;
        } else if protein_terminus == "n" {
            value.protein_terminus = true;
            value.terminus = protein_terminus.to_owned();
        } else if protein_terminus == "N" {
            value.protein_terminus = false;
        }
        value.term = match value.terminus.as_str() {
            "n" if value.protein_terminus => Some(TermSpecificity::ProteinNTerm),
            "n" => Some(TermSpecificity::NTerm),
            "c" if value.protein_terminus => Some(TermSpecificity::ProteinCTerm),
            "c" => Some(TermSpecificity::CTerm),
            _ => None,
        };
        if value.mass == value.mass_diff {
            value.warnings.push(format!(
                "for modification with mass {mass}, mass == massdiff, which is wrong; recomputing \
                 the absolute mass from massdiff"
            ));
            let recomputed = match value.term {
                Some(TermSpecificity::NTerm | TermSpecificity::ProteinNTerm) => {
                    terminal_gain(true).map(|gain| value.mass_diff + gain)
                }
                Some(TermSpecificity::CTerm | TermSpecificity::ProteinCTerm) => {
                    terminal_gain(false).map(|gain| value.mass_diff + gain)
                }
                _ => value
                    .amino_acid
                    .chars()
                    .next()
                    .ok_or_else(|| bad("modification without terminus needs an amino acid"))
                    .and_then(internal_residue_mass)
                    .map(|residue| value.mass_diff + residue),
            };
            match recomputed {
                Ok(recomputed) if recomputed.is_finite() => value.mass = recomputed,
                // Source dereferences the ResidueDB lookup unconditionally and
                // therefore crashes on an unknown residue code; this keeps the
                // declared mass and records the failure instead.
                _ => value.warnings.push(
                    "the absolute mass could not be recomputed because the residue has no known \
                     monoisotopic mass; keeping the declared value"
                        .into(),
                ),
            }
        }
        value.resolve(options, registry)?;
        Ok(value)
    }

    fn resolve(&mut self, options: &ReadOptions, registry: &ModificationsDB) -> Result<()> {
        let preferred = if self.variable {
            &options.preferred_variable_modifications
        } else {
            &options.preferred_fixed_modifications
        };
        self.resolved = lookup_preferred(
            preferred,
            &self.amino_acid,
            self.mass_diff,
            &self.description,
            self.term,
        )
        .map(ResolvedModification::Known);
        if self.resolved.is_none() && !self.description.is_empty() {
            match registry.get_modification_handle(
                &self.description,
                self.amino_acid.chars().next(),
                self.term,
            ) {
                Ok(handle) => self.resolved = Some(ResolvedModification::Known(handle)),
                Err(_) => self.warnings.push(format!(
                    "modification '{}' of residue '{}' could not be matched; trying by \
                     modification mass",
                    self.description, self.amino_acid
                )),
            }
        } else if self.resolved.is_none() {
            self.warnings.push(
                "no modification description given; trying to define by modification mass".into(),
            );
        }
        if self.resolved.is_some() {
            return Ok(());
        }
        let origin = self.amino_acid.chars().next();
        let mut matches = Vec::new();
        // Source tries the least specific search first, but only when the
        // declaration left the terminus open.
        if self.term.is_none() {
            matches = sorted_by_mass(
                registry,
                self.mass_diff,
                origin,
                Some(TermSpecificity::Anywhere),
            )?;
        }
        if matches.is_empty() {
            matches = sorted_by_mass(registry, self.mass_diff, origin, self.term)?;
        }
        if let Some(first) = matches.first() {
            if matches.len() > 1 {
                let names: Vec<&str> = matches.iter().map(|m| m.full_id()).collect();
                self.warnings.push(format!(
                    "modification '{}' is not uniquely defined by the given data; using '{}' to \
                     represent any of '{}'",
                    crate::param::value::format_float(self.mass, true),
                    first.full_id(),
                    names.join(", ")
                ));
            }
            self.resolved = Some(ResolvedModification::Known(Arc::clone(first)));
            return Ok(());
        }
        if self.mass_diff == 0.0 {
            // Source leaves registered_mod_ null here, and the caller then drops
            // the declaration entirely.
            return Ok(());
        }
        let text = diff_mono_mass_string(self.mass_diff);
        let term = self.term.unwrap_or(TermSpecificity::Anywhere);
        let full_id = match term {
            TermSpecificity::NTerm | TermSpecificity::ProteinNTerm => format!(".n[{text}]"),
            TermSpecificity::CTerm | TermSpecificity::ProteinCTerm => format!(".c[{text}]"),
            TermSpecificity::Anywhere => {
                let origin = origin.ok_or_else(|| {
                    Error::InvalidValue(
                        "cannot create non-terminal mod without origin AA residue".into(),
                    )
                })?;
                format!("{origin}[{text}]")
            }
        };
        self.warnings.push(format!(
            "modification '{}/delta {}' is unknown; resuming with '{full_id}', which could lead to \
             failures using the data downstream",
            crate::param::value::format_float(self.mass, true),
            crate::param::value::format_float(self.mass_diff, true)
        ));
        self.resolved = Some(ResolvedModification::Mass {
            text,
            full_id,
            term,
        });
        Ok(())
    }
}

/// Preferred-modification lookup, as `lookupModInPreferredMods_`.
///
/// Matches the declared description against a full identifier first, then falls
/// back to residue, terminal specificity and a mass window.
fn lookup_preferred(
    preferred: &[Arc<ResidueModification>],
    amino_acid: &str,
    mass_diff: f64,
    description: &str,
    term: Option<TermSpecificity>,
) -> Option<Arc<ResidueModification>> {
    for candidate in preferred {
        if description == candidate.full_id() {
            return Some(Arc::clone(candidate));
        }
    }
    let origin = amino_acid.chars().next();
    for candidate in preferred {
        let residue_matches = match origin {
            None => true,
            Some(residue) => candidate.origin() == Some(residue),
        };
        let term_matches = term.is_none_or(|value| value == candidate.term_specificity());
        if residue_matches
            && term_matches
            && (mass_diff - candidate.diff_mono_mass()).abs() < MODIFICATION_TOLERANCE
        {
            return Some(Arc::clone(candidate));
        }
    }
    None
}

/// Registry mass search ordered as `searchModificationsByDiffMonoMassSorted`.
///
/// The source keys a `std::map` on `(|massdiff error|, insertion counter)`, so
/// equal errors keep registry order; [`slice::sort_by`] is stable and the
/// registry search already yields registry order, which reproduces that.
fn sorted_by_mass(
    registry: &ModificationsDB,
    mass: f64,
    origin: Option<char>,
    term: Option<TermSpecificity>,
) -> Result<Vec<Arc<ResidueModification>>> {
    let mut matches: Vec<Arc<ResidueModification>> =
        registry.search_by_mass_handles(mass, MODIFICATION_TOLERANCE, origin, term)?;
    matches.sort_by(|a, b| {
        (a.diff_mono_mass() - mass)
            .abs()
            .total_cmp(&(b.diff_mono_mass() - mass).abs())
    });
    Ok(matches)
}

fn finite(value: &str, label: &str) -> Result<f64> {
    value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| bad(format!("invalid finite {label}: {value:?}")))
}

// ---------------------------------------------------------------------------
// Bounded XML input
// ---------------------------------------------------------------------------

fn xml_char(c: char) -> bool {
    matches!(c, '\t' | '\n' | '\r')
        || ('\u{20}'..='\u{d7ff}').contains(&c)
        || ('\u{e000}'..='\u{fffd}').contains(&c)
        || c >= '\u{10000}'
}
fn space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\r' | b'\n')
}

/// Read the whole input, bounded, and normalise line endings.
fn document(mut input: impl BufRead, options: &ReadOptions, meter: &mut Meter) -> Result<String> {
    let mut bytes = Vec::new();
    loop {
        let chunk = input.fill_buf()?;
        if chunk.is_empty() {
            break;
        }
        let remaining = options.max_input_bytes.saturating_sub(bytes.len());
        let take = chunk.len().min(1 << 16).min(remaining.saturating_add(1));
        meter.spend(take, 0)?;
        if take > remaining {
            return Err(bad("pepXML input byte limit exceeded"));
        }
        let wanted = bytes.len().checked_add(take).ok_or_else(|| meter.limit())?;
        if wanted > bytes.capacity() {
            let capacity = bytes
                .capacity()
                .saturating_mul(2)
                .max(wanted)
                .min(options.max_input_bytes);
            meter.spend(bytes.len(), capacity)?;
            bytes.reserve_exact(capacity.saturating_sub(bytes.len()));
        }
        bytes.extend_from_slice(&chunk[..take]);
        input.consume(take);
    }
    meter.spend(bytes.len().saturating_mul(4), bytes.len().saturating_mul(3))?;
    let content = bytes.strip_prefix(&[239, 187, 191]).unwrap_or(&bytes);
    let text = std::str::from_utf8(content)
        .map_err(|_| unsupported("pepXML input requires UTF-8 or ASCII-compatible bytes"))?;
    if !text.chars().all(xml_char) {
        return Err(bad("invalid XML 1.0 character"));
    }
    meter.spend(text.len().saturating_mul(2), text.len().saturating_mul(2))?;
    let text = text.replace("\r\n", "\n").replace('\r', "\n");
    meter.cap(text.len(), options.max_input_bytes)?;
    Ok(text)
}

/// Attribute list of one element, entity-expanded and whitespace-normalised.
type Attributes = Vec<(String, String)>;

fn attribute_spacing(element: &BytesStart<'_>) -> Result<()> {
    let tail = element.attributes_raw();
    if tail.first().is_some_and(|b| !space(*b)) {
        return Err(bad("missing XML attribute separator"));
    }
    let mut quote = None;
    let mut closed = false;
    for &b in tail {
        match quote {
            Some(q) if b == q => {
                quote = None;
                closed = true;
            }
            Some(_) => {}
            None => {
                if closed && !space(b) {
                    return Err(bad("missing XML attribute separator"));
                }
                closed = false;
                if b == b'\'' || b == b'"' {
                    quote = Some(b);
                }
            }
        }
    }
    Ok(())
}

fn attributes(element: &BytesStart<'_>, meter: &mut Meter) -> Result<Attributes> {
    attribute_spacing(element)?;
    let mut out: Attributes = Vec::new();
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|e| bad(e.to_string()))?;
        let key = std::str::from_utf8(attribute.key.as_ref())
            .map_err(|_| bad("invalid XML attribute name"))?;
        for (previous, _) in &out {
            meter.spend(previous.len().saturating_add(key.len()), 0)?;
            if previous == key {
                return Err(bad("duplicate XML attribute"));
            }
        }
        if attribute.value.contains(&b'<') {
            return Err(bad("raw '<' in XML attribute"));
        }
        let raw =
            std::str::from_utf8(&attribute.value).map_err(|_| bad("invalid XML attribute text"))?;
        meter.spend(raw.len().saturating_mul(4), raw.len().saturating_mul(4))?;
        let mut normalized = String::with_capacity(raw.len());
        for c in raw.chars() {
            normalized.push(if matches!(c, '\t' | '\n' | '\r') {
                ' '
            } else {
                c
            });
        }
        let value = quick_xml::escape::unescape(&normalized).map_err(|e| bad(e.to_string()))?;
        if !value.chars().all(xml_char) {
            return Err(bad("invalid XML 1.0 character in attribute"));
        }
        out.push((key.to_owned(), value.into_owned()));
    }
    Ok(out)
}

fn attribute<'a>(attributes: &'a Attributes, name: &str) -> Option<&'a str> {
    attributes
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value.as_str())
}
fn required<'a>(attributes: &'a Attributes, element: &str, name: &str) -> Result<&'a str> {
    attribute(attributes, name).ok_or_else(|| {
        bad(format!(
            "'{element}' element requires the '{name}' attribute"
        ))
    })
}
fn required_f64(attributes: &Attributes, element: &str, name: &str) -> Result<f64> {
    finite(required(attributes, element, name)?, name)
}
fn optional_f64(attributes: &Attributes, name: &str) -> Result<Option<f64>> {
    attribute(attributes, name)
        .map(|value| finite(value, name))
        .transpose()
}
fn integer<T: std::str::FromStr>(value: &str, name: &str) -> Result<T> {
    value
        .trim()
        .parse()
        .map_err(|_| bad(format!("invalid integer {name}: {value:?}")))
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

struct RunState {
    /// Indices into the document's protein identifications for this run.
    proteins: Vec<usize>,
    base_name: String,
    ms_run_path: String,
    enzyme: String,
    parameters: SearchParameters,
    precursor_tolerance: Option<f64>,
    precursor_tolerance_ppm: bool,
    fixed: Vec<HeaderModification>,
    variable: Vec<HeaderModification>,
    protein_hits: usize,
}
impl RunState {
    fn new() -> Self {
        Self {
            proteins: Vec::new(),
            base_name: String::new(),
            ms_run_path: String::new(),
            enzyme: "unknown_enzyme".into(),
            parameters: SearchParameters::default(),
            precursor_tolerance: None,
            precursor_tolerance_ppm: false,
            fixed: Vec::new(),
            variable: Vec::new(),
            protein_hits: 0,
        }
    }
}

struct ReaderState<'a> {
    options: &'a ReadOptions,
    registry: &'a ModificationsDB,
    document: PepXmlDocument,
    run: RunState,
    experiment_name: String,
    analysis_summary: bool,
    search_score_summary: bool,
    search_summary: bool,
    wrong_experiment: bool,
    seen_experiment: bool,
    checked_base_name: bool,
    has_decoys: bool,
    decoy_prefix: String,
    search_engine: String,
    date: DateTime,
    hydrogen_mass: f64,
    hydrogen_mono: f64,
    hydrogen_average: f64,
    native_spectrum_name: String,
    experiment_label: String,
    swath_assay: String,
    status: String,
    rt: Option<f64>,
    mz: Option<f64>,
    scan_number: u64,
    charge: i32,
    search_id: usize,
    current_peptide: PeptideIdentification,
    current_hit: PeptideHit,
    current_sequence: String,
    current_analysis: crate::identification::AnalysisResult,
    current_modifications: Vec<(ResolvedModification, usize)>,
}

impl<'a> ReaderState<'a> {
    fn new(options: &'a ReadOptions, registry: &'a ModificationsDB) -> Result<Self> {
        let hydrogen = EmpiricalFormula::parse("H")?;
        let experiment_name = strip_extension(&options.experiment_name).to_owned();
        Ok(Self {
            options,
            registry,
            document: PepXmlDocument::default(),
            run: RunState::new(),
            analysis_summary: false,
            search_score_summary: false,
            search_summary: false,
            wrong_experiment: false,
            seen_experiment: experiment_name.is_empty(),
            checked_base_name: experiment_name.is_empty(),
            experiment_name,
            has_decoys: false,
            decoy_prefix: String::new(),
            search_engine: String::new(),
            date: DateTime::default(),
            // Source assumes "average" until a search_summary says otherwise.
            hydrogen_mass: hydrogen.average_mass(),
            hydrogen_mono: hydrogen.mono_mass(),
            hydrogen_average: hydrogen.average_mass(),
            native_spectrum_name: String::new(),
            experiment_label: String::new(),
            swath_assay: String::new(),
            status: String::new(),
            rt: None,
            mz: None,
            scan_number: 0,
            charge: 0,
            search_id: 1,
            current_peptide: PeptideIdentification::default(),
            current_hit: PeptideHit::default(),
            current_sequence: String::new(),
            current_analysis: crate::identification::AnalysisResult::default(),
            current_modifications: Vec::new(),
        })
    }

    fn warn(&mut self, message: impl Into<String>) {
        if self.document.warnings.len() < self.options.max_warnings {
            self.document.warnings.push(message.into());
        } else {
            self.document.suppressed_warnings += 1;
        }
    }

    /// Index of the protein identification the current `search_id` refers to.
    ///
    /// pepXML numbers searches either per run or sequentially across runs, so
    /// `search_id` may exceed the number of runs recorded for this
    /// `msms_run_summary`; the source clamps with
    /// `min(current_proteins_.size(), search_id_) - 1`, which underflows for an
    /// empty run or for `search_id="0"`. This returns an explicit error instead.
    fn current_protein(&self) -> Result<usize> {
        let position = self
            .run
            .proteins
            .len()
            .min(self.search_id)
            .checked_sub(1)
            .ok_or_else(|| bad("'search_id' must be at least 1 and a run must be open"))?;
        self.run
            .proteins
            .get(position)
            .copied()
            .ok_or_else(|| bad("no identification run is open"))
    }

    fn read_rt_mz_charge(&mut self, attributes: &Attributes) -> Result<()> {
        let mass = required_f64(attributes, "spectrum_query", "precursor_neutral_mass")?;
        self.charge = integer(
            required(attributes, "spectrum_query", "assumed_charge")?,
            "assumed_charge",
        )?;
        if self.charge == 0 {
            // Source divides by the charge unconditionally, yielding infinity.
            return Err(bad("'assumed_charge' must not be zero"));
        }
        let mz = (mass + self.hydrogen_mass * f64::from(self.charge)) / f64::from(self.charge);
        if !mz.is_finite() {
            return Err(bad("recomputed precursor m/z is not finite"));
        }
        self.mz = Some(mz);
        self.rt = None;
        // Source ignores "end_scan" and assumes a single scan; it also compares
        // start_scan with itself, so its "endscan not equal to startscan" error
        // is unreachable. This port compares the two attributes as intended.
        self.scan_number = integer(
            required(attributes, "spectrum_query", "start_scan")?,
            "start_scan",
        )?;
        if let Some(end) = attribute(attributes, "end_scan") {
            let end: u64 = integer(end, "end_scan")?;
            if end != self.scan_number {
                self.warn(
                    "endscan not equal to startscan. Merged spectrum queries not supported. \
                     Parsing start scan nr. only.",
                );
            }
        }
        if let Some(rt) = optional_f64(attributes, "retention_time_sec")? {
            self.rt = Some(rt);
            return Ok(());
        }
        if self.options.lookup.is_empty() {
            self.warn("Cannot get RT information - no spectra given");
            return Ok(());
        }
        let index = if self.scan_number != 0 {
            self.options.lookup.find_by_scan_number(self.scan_number)
        } else {
            attribute(attributes, "spectrum")
                .and_then(|value| self.options.lookup.find_by_native_id(value))
        };
        match index.and_then(|index| self.options.lookup.get(index)) {
            Some(meta) if meta.ms_level == 2 => self.rt = Some(meta.rt),
            _ => self.warn("Cannot get RT information - scan mapping is incorrect"),
        }
        Ok(())
    }

    fn start(&mut self, element: &str, attributes: &Attributes, meter: &mut Meter) -> Result<()> {
        if element == "msms_run_summary" {
            return self.start_run(attributes);
        }
        if element == "analysis_summary" {
            // This element can nest "search_summary" elements, which are only
            // expected under "msms_run_summary", so the whole subtree is skipped.
            self.analysis_summary = true;
            return Ok(());
        }
        if self.wrong_experiment || self.analysis_summary {
            return Ok(());
        }
        match element {
            "search_score" => self.start_search_score(attributes),
            "search_hit" => self.start_search_hit(attributes, meter),
            "search_result" => self.start_search_result(attributes),
            "spectrum_query" => self.start_spectrum_query(attributes),
            "analysis_result" => {
                self.current_analysis = crate::identification::AnalysisResult::default();
                self.current_analysis.score_type =
                    required(attributes, "analysis_result", "analysis")?.to_owned();
                Ok(())
            }
            "search_score_summary" => {
                self.search_score_summary = true;
                Ok(())
            }
            "parameter" => self.start_parameter(attributes),
            "peptideprophet_result" => {
                let value = required_f64(attributes, "peptideprophet_result", "probability")?;
                if self.current_peptide.score_type != "InterProphet probability" {
                    self.current_hit.score = value;
                    self.current_peptide.score_type = "PeptideProphet probability".into();
                    self.current_peptide.higher_score_better = true;
                }
                self.current_analysis.main_score = value;
                self.current_analysis.higher_is_better = true;
                Ok(())
            }
            "interprophet_result" => {
                let value = required_f64(attributes, "interprophet_result", "probability")?;
                self.current_hit.score = value;
                self.current_peptide.score_type = "InterProphet probability".into();
                self.current_peptide.higher_score_better = true;
                self.current_analysis.main_score = value;
                self.current_analysis.higher_is_better = true;
                Ok(())
            }
            "modification_info" => self.start_modification_info(attributes, meter),
            "alternative_protein" => self.start_alternative_protein(attributes, meter),
            "mod_aminoacid_mass" => self.start_mod_aminoacid_mass(attributes, meter),
            "aminoacid_modification" | "terminal_modification" => {
                self.start_modification_declaration(element, attributes, meter)
            }
            "search_summary" => self.start_search_summary(attributes),
            "sample_enzyme" => {
                // Special case: a search parameter that occurs *before*
                // "search_summary".
                let mut name = required(attributes, "sample_enzyme", "name")?.to_owned();
                if name == "stricttrypsin" {
                    name = "Trypsin/P".into(); // MSFragger synonym
                }
                self.run.enzyme = name;
                if let Some(enzyme) = find_enzyme(&self.run.enzyme) {
                    self.run.parameters.digestion_enzyme = enzyme.name().into();
                }
                Ok(())
            }
            "specificity" if self.run.parameters.digestion_enzyme == "unknown_enzyme" => {
                let cut = required(attributes, "specificity", "cut")?;
                let no_cut = attribute(attributes, "no_cut").unwrap_or("");
                let sense = required(attributes, "specificity", "sense")?;
                self.run.parameters.digestion_enzyme =
                    format!("user-defined,{},{cut},{no_cut},{sense}", self.run.enzyme);
                self.run.parameters.digestion_regex = cut.to_owned();
                Ok(())
            }
            "enzymatic_search_constraint" => self.start_enzymatic_constraint(attributes),
            "search_database" => {
                let path = required(attributes, "search_database", "local_path")?;
                self.run.parameters.database = if path.is_empty() {
                    attribute(attributes, "database_name")
                        .unwrap_or("")
                        .to_owned()
                } else {
                    path.to_owned()
                };
                Ok(())
            }
            "msms_pipeline_analysis" => {
                let raw = required(attributes, "msms_pipeline_analysis", "date")?;
                let mut date = raw.to_owned();
                // Source repairs a corrupted xs:dateTime of the shape
                // "yyyy:MM:dd:hh:mm:ss" by indexing bytes 4, 7 and 10 without a
                // length check; a shorter value reads out of bounds there.
                let bytes = date.as_bytes();
                if bytes.len() > 10 && bytes[4] == b':' && bytes[7] == b':' && bytes[10] == b':' {
                    self.warn(
                        "Format of attribute 'date' in tag 'msms_pipeline_analysis' does not \
                         comply with standard 'xs:dateTime'",
                    );
                    date.replace_range(4..5, "-");
                    date.replace_range(7..8, "-");
                    date.replace_range(10..11, "T");
                }
                // Source asDateTime_ throws on an unparsable value; DateTime::set
                // clears itself instead, which the empty-date identifiers below
                // then reflect.
                self.date = DateTime::parse(&date).unwrap_or_default();
                Ok(())
            }
            _ => Ok(()),
        }
    }

    fn start_run(&mut self, attributes: &Attributes) -> Result<()> {
        let base_name = attribute(attributes, "base_name").unwrap_or("").to_owned();
        self.run.ms_run_path.clear();
        if !self.experiment_name.is_empty() {
            if base_name.is_empty() {
                // Really should not happen, but does for Mascot pepXML exports.
                self.warn("'base_name' attribute of 'msms_run_summary' element is empty");
                self.wrong_experiment = false;
                self.checked_base_name = false;
            } else {
                self.wrong_experiment = !base_name.ends_with(&self.experiment_name);
                self.seen_experiment = self.seen_experiment || !self.wrong_experiment;
                self.checked_base_name = true;
            }
        }
        if self.wrong_experiment {
            return Ok(());
        }
        // Compose this run's primary MS run path from base_name and raw_data,
        // e.g. base_name="run1" + raw_data=".mzML" -> "run1.mzML". raw_data
        // already carries a leading dot in TPP and OpenMS exports, so one is
        // only added when missing.
        if !base_name.is_empty() {
            let raw = attribute(attributes, "raw_data").unwrap_or("");
            self.run.ms_run_path = if raw.is_empty() || raw.starts_with('.') {
                format!("{base_name}{raw}")
            } else {
                format!("{base_name}.{raw}")
            };
        }
        // Create a ProteinIdentification in case "search_summary" is missing.
        let mut protein = ProteinIdentification::new();
        protein.date_time = date_time(&self.date);
        protein.identifier = format!("unknown_{}", self.date.date_string());
        self.run.enzyme = "unknown_enzyme".into();
        if !self.run.ms_run_path.is_empty() {
            protein.primary_ms_run_paths = vec![self.run.ms_run_path.clone()];
        }
        self.document.protein_identifications.push(protein);
        self.run.proteins.clear();
        self.run
            .proteins
            .push(self.document.protein_identifications.len() - 1);
        self.run.protein_hits = 0;
        Ok(())
    }

    fn start_search_score(&mut self, attributes: &Attributes) -> Result<()> {
        let name = required(attributes, "search_score", "name")?.to_owned();
        let number = |attributes: &Attributes| required_f64(attributes, "search_score", "value");
        match name.as_str() {
            // X!Tandem, Mascot or MSFragger E-value.
            "expect" => {
                let value = number(attributes)?;
                self.current_hit.score = value;
                self.current_peptide.score_type = name.clone();
                self.current_peptide.higher_score_better = false;
                let accession = match self.search_engine.as_str() {
                    "Comet" => Some("MS:1002257"),
                    // No separate CV terms are known for these two.
                    "X! Tandem" | "MSFragger" => Some("MS:1001330"),
                    "Mascot" => Some("MS:1001172"),
                    // There is no generic umbrella term for an expectation value.
                    _ => None,
                };
                if let Some(accession) = accession {
                    self.set_hit_meta(accession, value)?;
                }
            }
            "mvh" => {
                // MyriMatch score.
                self.current_hit.score = number(attributes)?;
                self.current_peptide.score_type = name.clone();
                self.current_peptide.higher_score_better = true;
            }
            "xcorr" => {
                let value = number(attributes)?;
                // MyriMatch also reports an xcorr, which the source ignores.
                if self.search_engine != "MyriMatch" {
                    self.current_hit.score = value;
                    self.current_peptide.score_type = name.clone();
                    self.current_peptide.higher_score_better = true;
                }
                let accession = if self.search_engine == "Comet" {
                    "MS:1002252"
                } else {
                    // No other or generic xcorr term exists; SEQUEST's is used.
                    "MS:1001155"
                };
                self.set_hit_meta(accession, value)?;
            }
            "fval" => {
                // SpectraST score.
                let value = number(attributes)?;
                self.current_hit.score = value;
                self.current_peptide.score_type = name.clone();
                self.current_peptide.higher_score_better = true;
                self.set_hit_meta("MS:1001419", value)?;
            }
            "hyperscore" | "nextscore" => {
                let value = number(attributes)?;
                self.set_hit_meta(&name, value)?;
            }
            _ if self.search_engine == "Comet" => {
                let accessions: &[(&str, &[&str])] = &[
                    ("deltacn", &["MS:1002253", "COMET:deltaCn"]),
                    ("spscore", &["MS:1002255"]),
                    ("sprank", &["MS:1002256"]),
                    ("deltacnstar", &["MS:1002254"]),
                    ("lnrSp", &["COMET:lnRankSP"]),
                    ("deltLCn", &["COMET:deltaLCn"]),
                    ("lnExpect", &["COMET:lnExpect"]),
                    // matched_ions / total_ions
                    ("IonFrac", &["COMET:IonFrac"]),
                    ("lnNumSP", &["COMET:lnNumSP"]),
                ];
                if let Some((_, keys)) = accessions.iter().find(|(key, _)| *key == name) {
                    let value = number(attributes)?;
                    for key in *keys {
                        self.set_hit_meta(key, value)?;
                    }
                }
            }
            _ if self.options.parse_unknown_scores => {
                let raw = required(attributes, "search_score", "value")?.to_owned();
                if name.ends_with("_ions") {
                    let value: i64 = integer(&raw, "search_score value")?;
                    self.current_hit.metadata.insert(name, value.into());
                } else {
                    // Source falls back to the raw string when the value does not
                    // convert; the numeric attempt is not reported.
                    match finite(&raw, "search_score value") {
                        Ok(value) => self.set_hit_meta(&name, value)?,
                        Err(_) => {
                            self.current_hit.metadata.insert(name, raw.into());
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn set_hit_meta(&mut self, key: &str, value: f64) -> Result<()> {
        self.current_hit
            .metadata
            .insert(key.to_owned(), MetaValue::try_from(value)?);
        Ok(())
    }

    fn start_search_hit(&mut self, attributes: &Attributes, meter: &mut Meter) -> Result<()> {
        meter.cap(self.current_peptide.hits.len() + 1, self.options.max_hits)?;
        meter.spend(256, 512)?;
        self.current_sequence = required(attributes, "search_hit", "peptide")?.to_owned();
        self.current_modifications.clear();
        self.current_hit = PeptideHit::default();
        let rank: u32 = integer(required(attributes, "search_hit", "hit_rank")?, "hit_rank")?;
        // Ranks are 1-based in pepXML and 0-based in OpenMS. The source stores
        // rank - 1 in an unsigned field, so hit_rank="0" wraps.
        self.current_hit.rank = rank
            .checked_sub(1)
            .ok_or_else(|| bad("'hit_rank' is 1-based and must not be zero"))?;
        self.current_hit.charge = self.charge;
        let mut aa_before = FlankingResidue::default();
        let mut aa_after = FlankingResidue::default();
        if let Some(value) = attribute(attributes, "peptide_prev_aa") {
            aa_before = flanking(value, true)?;
        }
        if let Some(value) = attribute(attributes, "peptide_next_aa") {
            aa_after = flanking(value, false)?;
        }
        if self.search_engine == "Comet" {
            for (attribute_name, key) in [
                ("num_matched_ions", "MS:1002258"),
                ("tot_num_ions", "MS:1002259"),
                ("num_matched_peptides", "num_matched_peptides"),
            ] {
                if let Some(value) = attribute(attributes, attribute_name) {
                    self.current_hit
                        .metadata
                        .insert(key.to_owned(), value.into());
                }
            }
        }
        let mass_difference = required_f64(attributes, "search_hit", "massdiff")?;
        let isotope_error = (mass_difference / C13C12_MASSDIFF_U).round();
        if !isotope_error.is_finite() || isotope_error.abs() > 9.0e15 {
            return Err(bad("'massdiff' implies an unrepresentable isotope error"));
        }
        self.current_hit.metadata.insert(
            "isotope_error".into(),
            MetaValue::from(isotope_error as i64),
        );
        let protein = required(attributes, "search_hit", "protein")?
            .trim()
            .to_owned();
        let evidence = PeptideEvidence {
            protein_accession: protein.clone(),
            aa_before,
            aa_after,
            ..Default::default()
        };
        let mut hit = ProteinHit {
            accession: protein.clone(),
            ..Default::default()
        };
        if self.has_decoys {
            self.annotate_decoy(&protein, &mut hit)?;
        }
        self.current_hit.evidences.push(evidence);
        let index = self.current_protein()?;
        self.insert_protein_hit(index, hit, meter)
    }

    fn annotate_decoy(&mut self, protein: &str, hit: &mut ProteinHit) -> Result<()> {
        let decoy = protein.starts_with(&self.decoy_prefix);
        let current = self.current_hit.target_decoy_type()?;
        let updated = match current {
            TargetDecoyType::Unknown if decoy => Some(TargetDecoyType::Decoy),
            TargetDecoyType::Unknown => Some(TargetDecoyType::Target),
            TargetDecoyType::Target if decoy => Some(TargetDecoyType::TargetAndDecoy),
            TargetDecoyType::Decoy if !decoy => Some(TargetDecoyType::TargetAndDecoy),
            _ => None,
        };
        if let Some(updated) = updated {
            self.current_hit.set_target_decoy_type(updated);
        }
        hit.set_target_decoy_type(if decoy {
            TargetDecoyType::Decoy
        } else {
            TargetDecoyType::Target
        })
    }

    fn insert_protein_hit(
        &mut self,
        index: usize,
        hit: ProteinHit,
        meter: &mut Meter,
    ) -> Result<()> {
        self.run.protein_hits += 1;
        meter.cap(self.run.protein_hits, self.options.max_protein_hits)?;
        meter.spend(
            hit.accession.len().saturating_add(64),
            hit.accession.len().saturating_mul(2).saturating_add(256),
        )?;
        self.document
            .protein_identifications
            .get_mut(index)
            .ok_or_else(|| bad("no identification run is open"))?
            .hits
            .push(hit);
        Ok(())
    }

    fn start_search_result(&mut self, attributes: &Attributes) -> Result<()> {
        self.current_peptide = PeptideIdentification::default();
        self.current_peptide.rt = self.rt;
        self.current_peptide.mz = self.mz;
        // The source MS run is recorded on this run's ProteinIdentification; a
        // PeptideIdentification is linked to it through the shared identifier.
        self.search_id = 1;
        if let Some(value) = attribute(attributes, "search_id") {
            self.search_id = integer(value, "search_id")?;
        }
        let index = self.current_protein()?;
        self.current_peptide.identifier = self
            .document
            .protein_identifications
            .get(index)
            .ok_or_else(|| bad("no identification run is open"))?
            .identifier
            .clone();
        if !self.native_spectrum_name.is_empty() && self.options.keep_native_spectrum_name {
            self.current_peptide.metadata.insert(
                "pepxml_spectrum_name".into(),
                self.native_spectrum_name.clone().into(),
            );
        }
        if SpectrumIndex::is_native_id(&self.native_spectrum_name) {
            self.current_peptide
                .set_spectrum_reference(self.native_spectrum_name.clone());
        } else if self.scan_number != 0 {
            self.current_peptide
                .set_spectrum_reference(format!("scan={}", self.scan_number));
        }
        if !self.experiment_label.is_empty() {
            self.current_peptide
                .set_experiment_label(self.experiment_label.clone());
        }
        if !self.swath_assay.is_empty() {
            self.current_peptide
                .metadata
                .insert("swath_assay".into(), self.swath_assay.clone().into());
        }
        if !self.status.is_empty() {
            self.current_peptide
                .metadata
                .insert("status".into(), self.status.clone().into());
        }
        Ok(())
    }

    fn start_spectrum_query(&mut self, attributes: &Attributes) -> Result<()> {
        self.read_rt_mz_charge(attributes)?;
        self.native_spectrum_name.clear();
        self.experiment_label.clear();
        self.swath_assay.clear();
        self.status.clear();
        // Later assignments win: "spectrumNativeID" is preferred over
        // "spectrum", and MSFragger's "native_id" over both.
        for name in ["spectrum", "spectrumNativeID", "native_id"] {
            if let Some(value) = attribute(attributes, name) {
                self.native_spectrum_name = value.to_owned();
            }
        }
        if let Some(value) = attribute(attributes, "experiment_label") {
            self.experiment_label = value.to_owned();
        }
        if let Some(value) = attribute(attributes, "swath_assay") {
            self.swath_assay = value.to_owned();
        }
        if let Some(value) = attribute(attributes, "status") {
            self.status = value.to_owned();
        }
        Ok(())
    }

    fn start_parameter(&mut self, attributes: &Attributes) -> Result<()> {
        let name = required(attributes, "parameter", "name")?;
        if self.search_score_summary {
            let value = required_f64(attributes, "parameter", "value")?;
            self.current_analysis
                .sub_scores
                .insert(name.to_owned(), value);
            return Ok(());
        }
        if !self.search_summary {
            // Every other parameter block is currently not handled.
            return Ok(());
        }
        match name {
            "fragment_bin_tol" => {
                let value = required_f64(attributes, "parameter", "value")?;
                if value < 0.0 {
                    return Err(bad("'fragment_bin_tol' must not be negative"));
                }
                self.run.parameters.fragment_tolerance = Tolerance::Absolute(value / 2.0);
            }
            "peptide_mass_tolerance" => {
                let value = required_f64(attributes, "parameter", "value")?;
                if value < 0.0 {
                    return Err(bad("'peptide_mass_tolerance' must not be negative"));
                }
                self.run.precursor_tolerance = Some(value);
            }
            // Comet-specific, as is the parameter name itself.
            "peptide_mass_units" => {
                let value: i32 = integer(required(attributes, "parameter", "value")?, "value")?;
                match value {
                    // 0 = amu, 1 = mmu, 2 = ppm.
                    0 | 1 => self.run.precursor_tolerance_ppm = false,
                    2 => self.run.precursor_tolerance_ppm = true,
                    _ => {}
                }
            }
            "decoy_search" => {
                let value: i32 = integer(required(attributes, "parameter", "value")?, "value")?;
                self.has_decoys = value != 0;
            }
            "decoy_prefix" => {
                self.decoy_prefix = required(attributes, "parameter", "value")?.to_owned();
            }
            _ => {}
        }
        Ok(())
    }

    fn start_modification_info(
        &mut self,
        attributes: &Attributes,
        meter: &mut Meter,
    ) -> Result<()> {
        for n_terminal in [true, false] {
            let name = if n_terminal {
                "mod_nterm_mass"
            } else {
                "mod_cterm_mass"
            };
            let Some(mass) = optional_f64(attributes, name)? else {
                continue;
            };
            let wanted = if n_terminal { "n" } else { "c" };
            let mut found = None;
            for source in [&self.run.variable, &self.run.fixed] {
                for candidate in source {
                    if (mass - candidate.mass).abs() < MODIFICATION_TOLERANCE
                        && candidate.terminus == wanted
                    {
                        found = candidate.resolved.clone();
                        break;
                    }
                }
                if found.is_some() {
                    break;
                }
            }
            if found.is_none() {
                // Not declared in the pepXML header: search the registry by the
                // difference against the bare terminus mass. The source tries the
                // peptide terminus and, if that throws, the protein terminus.
                let difference = mass - terminal_gain(n_terminal)?;
                let terms = if n_terminal {
                    [TermSpecificity::NTerm, TermSpecificity::ProteinNTerm]
                } else {
                    [TermSpecificity::CTerm, TermSpecificity::ProteinCTerm]
                };
                // The residue carrying the terminus, used to prefer a candidate
                // that can actually sit there. The source searches with an empty
                // residue and takes the first match, which may name a
                // specificity for a different amino acid entirely.
                let residue = if n_terminal {
                    self.current_sequence.chars().next()
                } else {
                    self.current_sequence.chars().next_back()
                };
                for term in terms {
                    let matches = self.registry.search_by_mass_handles(
                        difference,
                        MODIFICATION_TOLERANCE,
                        None,
                        Some(term),
                    )?;
                    let compatible =
                        matches
                            .iter()
                            .find(|candidate| match (candidate.origin(), residue) {
                                (None, _) => true,
                                (Some(origin), Some(residue)) => origin == residue || origin == 'X',
                                (Some(_), None) => false,
                            });
                    if let Some(first) = compatible.or_else(|| matches.first()) {
                        found = Some(ResolvedModification::Known(Arc::clone(first)));
                        break;
                    }
                }
            }
            match found {
                // The position is irrelevant for a terminal modification; the
                // source stores Size(-1) and never reads it.
                Some(resolved) => {
                    meter.cap(
                        self.current_modifications.len() + 1,
                        self.options.max_modifications,
                    )?;
                    meter.spend(64, 256)?;
                    self.current_modifications.push((resolved, usize::MAX));
                }
                None => {
                    let terminus = if n_terminal { "N" } else { "C" };
                    let text = crate::param::value::format_float(mass, true);
                    self.warn(format!(
                        "Cannot find {terminus}-terminal modification with mass {text}."
                    ));
                }
            }
        }
        Ok(())
    }

    fn start_alternative_protein(
        &mut self,
        attributes: &Attributes,
        meter: &mut Meter,
    ) -> Result<()> {
        let protein = required(attributes, "alternative_protein", "protein")?.to_owned();
        let evidence = PeptideEvidence {
            protein_accession: protein.clone(),
            ..Default::default()
        };
        let mut hit = ProteinHit {
            accession: protein.clone(),
            ..Default::default()
        };
        if self.has_decoys {
            self.annotate_decoy(&protein, &mut hit)?;
        }
        self.current_hit.evidences.push(evidence);
        let index = self.current_protein()?;
        self.insert_protein_hit(index, hit, meter)
    }

    fn start_mod_aminoacid_mass(
        &mut self,
        attributes: &Attributes,
        meter: &mut Meter,
    ) -> Result<()> {
        // This element should only carry internal residue modifications, or a
        // terminal modification at a specific amino acid (a pepXML limitation).
        let mass = required_f64(attributes, "mod_aminoacid_mass", "mass")?;
        let position: usize = integer(
            required(attributes, "mod_aminoacid_mass", "position")?,
            "position",
        )?;
        // Positions are 1-based. The source subtracts one from an unsigned value
        // and then indexes the sequence, so position="0" or a position past the
        // peptide reads out of bounds.
        let index = position
            .checked_sub(1)
            .ok_or_else(|| bad("'position' is 1-based and must not be zero"))?;
        let origin = self
            .current_sequence
            .chars()
            .nth(index)
            .ok_or_else(|| bad("'position' is outside the peptide sequence"))?;
        // The source cannot infer fixed vs variable from pepXML reliably, so it
        // tries the fixed declarations first and then the variable ones.
        let mut found = lookup_from_header(&self.run.fixed, mass, origin);
        if found.is_none() {
            found = lookup_from_header(&self.run.variable, mass, origin);
        }
        if found.is_none() {
            if let Some(psi_mod) = attribute(attributes, "id").filter(|v| !v.is_empty()) {
                if let Ok(handle) =
                    self.registry
                        .get_modification_handle(psi_mod, Some(origin), None)
                {
                    found = Some(ResolvedModification::Known(handle));
                }
            }
        }
        if found.is_none() {
            let difference = mass - internal_residue_mass(origin)?;
            // Try the least specific search first.
            let mut matches = self.registry.search_by_mass_handles(
                difference,
                MODIFICATION_TOLERANCE,
                Some(origin),
                Some(TermSpecificity::Anywhere),
            )?;
            if matches.is_empty() {
                let terms: &[TermSpecificity] = if position == 1 {
                    &[TermSpecificity::NTerm, TermSpecificity::ProteinNTerm]
                } else if position == self.current_sequence.chars().count() {
                    &[TermSpecificity::CTerm, TermSpecificity::ProteinCTerm]
                } else {
                    &[]
                };
                for term in terms {
                    matches = self.registry.search_by_mass_handles(
                        difference,
                        MODIFICATION_TOLERANCE,
                        Some(origin),
                        Some(*term),
                    )?;
                    if !matches.is_empty() {
                        break;
                    }
                }
            }
            match matches.first() {
                Some(first) => {
                    if matches.len() > 1 {
                        self.warn(format!(
                            "Modification '{}' of residue {origin} at position {position} in \
                             '{}' not registered in pepXML header nor uniquely defined in DB. \
                             Using {}",
                            crate::param::value::format_float(mass, true),
                            self.current_sequence,
                            first.full_id()
                        ));
                    }
                    found = Some(ResolvedModification::Known(Arc::clone(first)));
                }
                None => {
                    // Nothing matched: keep the declared mass as an anonymous
                    // annotation. Being attached to an amino acid it is probably
                    // not a terminal modification.
                    let text = crate::param::value::format_float(mass, true);
                    found = Some(ResolvedModification::Mass {
                        full_id: format!("{origin}[{text}]"),
                        text,
                        term: TermSpecificity::Anywhere,
                    });
                }
            }
        }
        if let Some(resolved) = found {
            meter.cap(
                self.current_modifications.len() + 1,
                self.options.max_modifications,
            )?;
            meter.spend(64, 256)?;
            self.current_modifications.push((resolved, index));
        }
        Ok(())
    }

    fn start_modification_declaration(
        &mut self,
        element: &str,
        attributes: &Attributes,
        meter: &mut Meter,
    ) -> Result<()> {
        let description = attribute(attributes, "description")
            .unwrap_or("")
            .to_owned();
        let mass_diff = required(attributes, element, "massdiff")?.to_owned();
        // The source reformats mass through its double parser and printer before
        // handing it on, which normalises the spelling.
        let mass =
            crate::param::value::format_float(required_f64(attributes, element, "mass")?, true);
        let variable = required(attributes, element, "variable")?.to_owned();
        let (amino_acid, terminus, protein_terminus) = if element == "aminoacid_modification" {
            // The specificity cannot be forced to ANYWHERE, because terminal
            // modifications may not be registered as such (notably by Comet).
            (
                required(attributes, element, "aminoacid")?.to_owned(),
                attribute(attributes, "peptide_terminus")
                    .unwrap_or("")
                    .to_owned(),
                attribute(attributes, "protein_terminus")
                    .unwrap_or("")
                    .to_owned(),
            )
        } else {
            // Very small fixed modifications - possibly the electron mass - are
            // annotated by some X!Tandem versions; they interfere with the real
            // modifications and are dropped.
            if finite(&mass_diff, "massdiff")?.abs() < XTANDEM_ARTIFICIAL_MODIFICATION_TOLERANCE {
                return Ok(());
            }
            (
                attribute(attributes, "aminoacid").unwrap_or("").to_owned(),
                required(attributes, element, "terminus")?.to_owned(),
                // Many search engines write "Y"/"N" here although the schema
                // allows only "", "n" or "c"; the constructor handles both.
                required(attributes, element, "protein_terminus")?.to_owned(),
            )
        };
        meter.cap(
            self.run.fixed.len() + self.run.variable.len() + 1,
            self.options.max_modifications,
        )?;
        meter.spend(1024, 2048)?;
        let modification = HeaderModification::new(
            &amino_acid,
            &mass_diff,
            &mass,
            &variable,
            &description,
            &terminus,
            &protein_terminus,
            self.options,
            self.registry,
        )?;
        if !modification.warnings.is_empty() {
            self.warn("Errors during parsing of aminoacid/terminal modification element:");
            for warning in modification.warnings.clone() {
                self.warn(warning);
            }
        }
        let Some(resolved) = modification.resolved.clone() else {
            return Ok(());
        };
        let description = resolved.full_id().to_owned();
        if modification.variable {
            self.run.variable.push(modification);
            self.run.parameters.variable_modifications.push(description);
        } else {
            self.run.fixed.push(modification);
            self.run.parameters.fixed_modifications.push(description);
        }
        Ok(())
    }

    fn start_search_summary(&mut self, attributes: &Attributes) -> Result<()> {
        self.search_summary = true;
        self.run.base_name = attribute(attributes, "base_name").unwrap_or("").to_owned();
        if !self.checked_base_name {
            // Work-around for files exported by Mascot.
            if self.run.base_name.ends_with(&self.experiment_name) {
                self.seen_experiment = true;
            } else {
                // Wrong experiment after all: roll back the run.
                self.document.protein_identifications.pop();
                self.run.proteins.clear();
                self.wrong_experiment = true;
                return Ok(());
            }
        }
        self.run.fixed.clear();
        self.run.variable.clear();
        let enzyme = find_enzyme(&self.run.enzyme);
        self.run.parameters = SearchParameters::default();
        self.run.precursor_tolerance = None;
        self.run.precursor_tolerance_ppm = false;
        if let Some(enzyme) = enzyme {
            self.run.parameters.digestion_enzyme = enzyme.name().into();
        }
        let mass_type = required(attributes, "search_summary", "precursor_mass_type")?;
        if mass_type == "monoisotopic" {
            self.hydrogen_mass = self.hydrogen_mono;
        } else {
            self.hydrogen_mass = self.hydrogen_average;
            if mass_type != "average" {
                self.warn(format!(
                    "'precursor_mass_type' attribute of 'search_summary' tag should be \
                     'monoisotopic' or 'average', not '{mass_type}' (assuming 'average')"
                ));
            }
        }
        // SearchParameters::mass_type is taken to refer to the fragment mass.
        let mass_type = required(attributes, "search_summary", "fragment_mass_type")?;
        match mass_type {
            "monoisotopic" => self.run.parameters.mass_type = PeakMassType::Monoisotopic,
            "average" => self.run.parameters.mass_type = PeakMassType::Average,
            other => self.warn(format!(
                "'fragment_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' \
                 or 'average', not '{other}'"
            )),
        }
        self.search_engine = required(attributes, "search_summary", "search_engine")?.to_owned();
        let version = attribute(attributes, "search_engine_version")
            .unwrap_or("")
            .to_owned();
        if self.search_engine == "X! Tandem" && version.starts_with("MSFragger") {
            self.search_engine = "MSFragger".into();
        }
        // Generate a unique identifier for every search engine run.
        let mut identifier = format!(
            "{}_{}_{}",
            self.search_engine,
            self.date.date_string(),
            self.date.time_string()
        );
        self.search_id = 1;
        if let Some(value) = attribute(attributes, "search_id") {
            self.search_id = integer(value, "search_id")?;
        }
        let index = if self.search_id <= self.document.protein_identifications.len() {
            // A ProteinIdentification was already created for "msms_run_summary".
            *self
                .run
                .proteins
                .last()
                .ok_or_else(|| bad("no identification run is open"))?
        } else {
            let mut protein = ProteinIdentification::new();
            protein.date_time = date_time(&self.date);
            self.document.protein_identifications.push(protein);
            identifier = format!("{identifier}_{}", self.search_id);
            let index = self.document.protein_identifications.len() - 1;
            self.run.proteins.push(index);
            index
        };
        let run = self
            .document
            .protein_identifications
            .get_mut(index)
            .ok_or_else(|| bad("no identification run is open"))?;
        run.search_engine = self.search_engine.clone();
        run.search_engine_version = version;
        run.identifier = identifier;
        // Record the source MS run - the spectra file - on the
        // ProteinIdentification. The path composed in "msms_run_summary" from
        // base_name and raw_data is preferred, because this element's base_name
        // can reference the identification result file rather than the spectra;
        // it is only used when the run summary had none, as in Mascot exports.
        let path = if self.run.ms_run_path.is_empty() {
            self.run.base_name.clone()
        } else {
            self.run.ms_run_path.clone()
        };
        if !path.is_empty() && run.primary_ms_run_paths.is_empty() {
            run.primary_ms_run_paths = vec![path];
        }
        Ok(())
    }

    fn start_enzymatic_constraint(&mut self, attributes: &Attributes) -> Result<()> {
        // Source note: the enzyme should not be overwritten here, but in most
        // files it is the same as in sample_enzyme or something useless such as
        // "default".
        let mut name = required(attributes, "enzymatic_search_constraint", "enzyme")?.to_owned();
        if name == "stricttrypsin" {
            name = "Trypsin/P".into(); // MSFragger synonym
        }
        self.run.enzyme = name;
        if let Some(enzyme) = find_enzyme(&self.run.enzyme) {
            if self.run.parameters.digestion_enzyme != "unknown_enzyme"
                && self.run.parameters.digestion_enzyme != enzyme.name()
            {
                self.warn(
                    "More than one enzyme found. This is currently not supported. Proceeding with \
                     last encountered only.",
                );
            }
            self.run.parameters.digestion_enzyme = enzyme.name().into();
        }
        // The additional information is always recorded, even when the enzymes
        // disagree, both for backwards compatibility and because it is the only
        // information available (MSFragger writes enzyme="default" here).
        let cleavages: i64 = integer(
            required(
                attributes,
                "enzymatic_search_constraint",
                "max_num_internal_cleavages",
            )?,
            "max_num_internal_cleavages",
        )?;
        self.run.parameters.missed_cleavages = u32::try_from(cleavages)
            .map_err(|_| bad("'max_num_internal_cleavages' must fit an unsigned 32-bit value"))?;
        let termini: i32 = integer(
            required(
                attributes,
                "enzymatic_search_constraint",
                "min_number_termini",
            )?,
            "min_number_termini",
        )?;
        // EnzymaticDigestion::Specificity is numbered so that engines can report
        // the number of required termini directly: 0 none, 1 semi, 2 full.
        self.run.parameters.enzyme_specificity = match termini {
            0 => EnzymeTermSpecificity::None,
            1 => EnzymeTermSpecificity::Semi,
            2 => EnzymeTermSpecificity::Full,
            _ => self.run.parameters.enzyme_specificity,
        };
        Ok(())
    }

    fn end(&mut self, element: &str, meter: &mut Meter) -> Result<()> {
        match element {
            "analysis_summary" => {
                self.analysis_summary = false;
                return Ok(());
            }
            "search_score_summary" => {
                self.search_score_summary = false;
                return Ok(());
            }
            "analysis_result" => {
                let result = std::mem::take(&mut self.current_analysis);
                self.current_hit.analysis_results.push(result);
                return Ok(());
            }
            _ => {}
        }
        if self.wrong_experiment || self.analysis_summary {
            return Ok(());
        }
        match element {
            "spectrum_query" => {
                self.native_spectrum_name.clear();
                self.experiment_label.clear();
                self.swath_assay.clear();
                self.status.clear();
            }
            "search_hit" => self.end_search_hit(meter)?,
            "search_result" => {
                meter.cap(
                    self.document.peptide_identifications.len() + 1,
                    self.options.max_identifications,
                )?;
                meter.spend(512, 1024)?;
                let peptide = std::mem::take(&mut self.current_peptide);
                self.document.peptide_identifications.push(peptide);
            }
            "search_summary" => {
                // idXML stores only the search engine and the date as a run
                // identifier, and two runs must differ. As a work-around for
                // multiple runs the date is advanced by one second per run.
                let (hour, minute, second) = self.date.time_components();
                let total = i64::from(hour) * 3600 + i64::from(minute) * 60 + i64::from(second) + 1;
                let total = total.rem_euclid(86_400);
                self.date.set_time_components(
                    u32::try_from(total / 3600).unwrap_or(0),
                    u32::try_from(total / 60 % 60).unwrap_or(0),
                    u32::try_from(total % 60).unwrap_or(0),
                )?;
                if let Some(tolerance) = self.run.precursor_tolerance {
                    self.run.parameters.precursor_tolerance = if self.run.precursor_tolerance_ppm {
                        Tolerance::Ppm(tolerance)
                    } else {
                        Tolerance::Absolute(tolerance)
                    };
                }
                let index = *self
                    .run
                    .proteins
                    .last()
                    .ok_or_else(|| bad("no identification run is open"))?;
                self.document
                    .protein_identifications
                    .get_mut(index)
                    .ok_or_else(|| bad("no identification run is open"))?
                    .search_parameters = self.run.parameters.clone();
                self.search_summary = false;
            }
            _ => {}
        }
        Ok(())
    }

    fn end_search_hit(&mut self, meter: &mut Meter) -> Result<()> {
        meter.spend(
            self.current_sequence.len().saturating_mul(64).max(1024),
            self.current_sequence.len().saturating_mul(64).max(4096),
        )?;
        let mut sequence = AASequence::parse_with_registry(&self.current_sequence, self.registry)?;
        // Walked once rather than per residue, so annotating stays linear in the
        // peptide length even for a sequence of non-ASCII placeholders.
        let residues: Vec<char> = sequence.as_str().chars().collect();
        meter.spend(residues.len(), residues.len().saturating_mul(4))?;
        // Applying AASequence parsing to the modification_info "modified_peptide"
        // attribute is not possible in general, because modifications may carry
        // symbols that would have to be stored and looked up separately.
        //
        // Modifications annotated on this search_hit take precedence over the
        // implicit fixed modifications applied afterwards.
        let modifications = std::mem::take(&mut self.current_modifications);
        let mut diagnostics = Vec::new();
        for (modification, position) in &modifications {
            if modification.is_n_terminal() {
                if sequence.n_terminal_modification().is_none() {
                    diagnostics.push(set_terminal(
                        &mut sequence,
                        modification,
                        true,
                        self.registry,
                    )?);
                } else {
                    self.warn(format!(
                        "Multiple N-term mods specified for search_hit with sequence {} \
                         proceeding with first.",
                        self.current_sequence
                    ));
                }
            } else if modification.is_c_terminal() {
                if sequence.c_terminal_modification().is_none() {
                    diagnostics.push(set_terminal(
                        &mut sequence,
                        modification,
                        false,
                        self.registry,
                    )?);
                } else {
                    self.warn(format!(
                        "Multiple C-term mods specified for search_hit with sequence {} \
                         proceeding with first.",
                        self.current_sequence
                    ));
                }
            } else if sequence.residue_modification(*position)?.is_none() {
                diagnostics.push(set_residue(
                    &mut sequence,
                    *position,
                    modification,
                    self.registry,
                )?);
            } else {
                self.warn(format!(
                    "Multiple mods for position {position} specified for search_hit with \
                     sequence {} proceeding with first.",
                    self.current_sequence
                ));
            }
        }
        // Apply the implicit fixed modifications wherever nothing is annotated.
        for declaration in self.run.fixed.clone() {
            let Some(modification) = declaration.resolved.as_ref() else {
                continue;
            };
            if modification.is_n_terminal() {
                if sequence.n_terminal_modification().is_none() {
                    diagnostics.push(set_terminal(
                        &mut sequence,
                        modification,
                        true,
                        self.registry,
                    )?);
                }
            } else if modification.is_c_terminal() {
                // Source tests C_TERM || PROTEIN_N_TERM here, a copy/paste slip
                // that leaves a fixed protein C-terminal modification to the
                // residue branch below; this port tests both C-terminal forms.
                if sequence.c_terminal_modification().is_none() {
                    diagnostics.push(set_terminal(
                        &mut sequence,
                        modification,
                        false,
                        self.registry,
                    )?);
                } else {
                    self.warn(format!(
                        "Trying to add a fixed C-term modification from the search_summary to an \
                         already annotated and modified C-terminus of {} ... skipping.",
                        self.current_sequence
                    ));
                }
            } else {
                for index in 0..sequence.len() {
                    if sequence.residue_modification(index)?.is_some() {
                        continue;
                    }
                    let residue = residues
                        .get(index)
                        .copied()
                        .ok_or_else(|| bad("peptide sequence changed while annotating"))?;
                    if declaration.amino_acid.contains(residue) {
                        diagnostics.push(set_residue(
                            &mut sequence,
                            index,
                            modification,
                            self.registry,
                        )?);
                    }
                }
            }
        }
        for diagnostic in diagnostics.into_iter().flatten() {
            self.warn(diagnostic);
        }
        self.current_hit.sequence = sequence;
        let hit = std::mem::take(&mut self.current_hit);
        hit.validate()?;
        self.current_peptide.hits.push(hit);
        Ok(())
    }
}

/// Attach a resolved modification to one residue.
///
/// Returns a diagnostic when the modification could not be attached by name and
/// its mass difference was retained as an anonymous annotation instead. The
/// source stores the resolved `ResidueModification*` directly and so never has
/// to re-resolve it; this port goes through the registry, which restricts a
/// named lookup to specificities valid at the annotated residue.
fn set_residue(
    sequence: &mut AASequence,
    index: usize,
    modification: &ResolvedModification,
    registry: &ModificationsDB,
) -> Result<Option<String>> {
    match modification {
        ResolvedModification::Mass { text, .. } => {
            sequence.set_mass_tag_with_registry(index, text, registry)?;
            Ok(None)
        }
        ResolvedModification::Known(value) => {
            if sequence
                .set_modification_with_registry(index, value.full_id(), registry)
                .is_ok()
            {
                return Ok(None);
            }
            let text = diff_mono_mass_string(value.diff_mono_mass());
            sequence.set_mass_tag_with_registry(index, &text, registry)?;
            Ok(Some(format!(
                "modification '{}' is not valid at the annotated residue; retaining its mass \
                 difference {text} instead",
                value.full_id()
            )))
        }
    }
}

/// Attach a resolved modification to a terminus, with the same fallback.
fn set_terminal(
    sequence: &mut AASequence,
    modification: &ResolvedModification,
    n_terminal: bool,
    registry: &ModificationsDB,
) -> Result<Option<String>> {
    let text = match modification {
        ResolvedModification::Mass { text, .. } => text.clone(),
        ResolvedModification::Known(value) => {
            let named = if n_terminal {
                sequence.set_n_terminal_modification_with_registry(value.full_id(), registry)
            } else {
                sequence.set_c_terminal_modification_with_registry(value.full_id(), registry)
            };
            if named.is_ok() {
                return Ok(None);
            }
            let text = diff_mono_mass_string(value.diff_mono_mass());
            if n_terminal {
                sequence.set_n_terminal_mass_tag_with_registry(&text, registry)?;
            } else {
                sequence.set_c_terminal_mass_tag_with_registry(&text, registry)?;
            }
            return Ok(Some(format!(
                "modification '{}' is not valid at the annotated terminus; retaining its mass \
                 difference {text} instead",
                value.full_id()
            )));
        }
    };
    if n_terminal {
        sequence.set_n_terminal_mass_tag_with_registry(&text, registry)?;
    } else {
        sequence.set_c_terminal_mass_tag_with_registry(&text, registry)?;
    }
    Ok(None)
}

/// Header lookup by absolute mass, as `PepXMLFile::lookupAddFromHeader_`.
///
/// The source allows the full modification tolerance against the header masses
/// as well, and requires the declaration's `aminoacid` string to contain the
/// annotated residue.
fn lookup_from_header(
    declarations: &[HeaderModification],
    mass: f64,
    residue: char,
) -> Option<ResolvedModification> {
    for declaration in declarations {
        if (mass - declaration.mass).abs() < MODIFICATION_TOLERANCE
            && declaration.amino_acid.contains(residue)
        {
            // Only one modification should match, so the search stops here.
            return declaration.resolved.clone();
        }
    }
    None
}

/// Flanking residue of a `peptide_prev_aa`/`peptide_next_aa` attribute.
///
/// pepXML spells a protein terminus `-`, which OpenMS's own marker characters
/// `[` and `]` record losslessly; the writer spells them `-` again. The source
/// stores `-` as the flanking residue character itself, so a stored pepXML
/// round-trips there while an idXML-sourced marker is written out verbatim.
fn flanking(value: &str, before: bool) -> Result<FlankingResidue> {
    // Source takes prev_aa[0] of a possibly empty attribute, which reads the
    // terminating null character; an empty attribute is an error here.
    let code = value
        .chars()
        .next()
        .ok_or_else(|| bad("flanking residue attribute must not be empty"))?;
    if code == '-' {
        return Ok(if before {
            FlankingResidue::NTerminus
        } else {
            FlankingResidue::CTerminus
        });
    }
    FlankingResidue::from_code(code)
}

/// pepXML spelling of a flanking residue, `-` for either protein terminus.
fn flanking_code(residue: FlankingResidue) -> char {
    match residue {
        FlankingResidue::NTerminus | FlankingResidue::CTerminus => '-',
        other => other.code(),
    }
}

fn date_time(date: &DateTime) -> Option<String> {
    if date.is_valid() {
        Some(date.iso_string())
    } else {
        None
    }
}

/// Enzyme lookup with the source's dual-key index.
///
/// `DigestionEnzymeDB` registers every enzyme under its name, its lower-cased
/// name and each synonym, so `name="trypsin"` resolves to `Trypsin`. The Rust
/// [`ProteaseDB`] matches names and synonyms exactly, so the lower-cased key is
/// applied here.
fn find_enzyme(name: &str) -> Option<&'static crate::chemistry::DigestionEnzymeProtein> {
    let database = ProteaseDB::global();
    if let Ok(enzyme) = database.get_enzyme(name) {
        return Some(enzyme);
    }
    database
        .enzymes()
        .iter()
        .find(|enzyme| enzyme.name().to_ascii_lowercase() == name)
}

/// Remove a known file extension, as `FileHandler::stripExtension`.
fn strip_extension(name: &str) -> &str {
    super::file_types::strip_extension(name)
}

// ---------------------------------------------------------------------------
// Public reading entry points
// ---------------------------------------------------------------------------

/// Read a pepXML document with the default options and the global registry.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn read(reader: impl BufRead) -> Result<PepXmlDocument> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read a pepXML document with explicit options and the global registry.
///
/// # Errors
///
/// See [`read_with_registry`].
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<PepXmlDocument> {
    read_with_registry(reader, options, ModificationsDB::global())
}

/// Read a pepXML document, resolving modifications against `registry`.
///
/// The returned sequences own shared handles to the registry's records, so the
/// registry may be dropped immediately after a successful read.
///
/// # Errors
///
/// Returns [`Error::Parse`] for malformed XML, a missing required attribute, a
/// non-finite number, a one-based index of zero, a position outside the peptide
/// or an exceeded resource ceiling; [`Error::Unsupported`] for a byte encoding
/// this port does not decode or for a DTD; [`Error::MissingInformation`] for a
/// modification declaration with neither an amino acid nor a terminus; and
/// [`Error::Io`] when the stream fails.
///
/// A nonempty [`ReadOptions::experiment_name`] that matches no run's `base_name`
/// is [`Error::Parse`], as the source's fatal
/// "Found no experiment with name" parse error.
pub fn read_with_registry(
    reader: impl BufRead,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<PepXmlDocument> {
    if options.max_depth == 0 || options.max_depth > 512 || options.max_elements == 0 {
        return Err(bad("invalid pepXML limits"));
    }
    let mut meter = Meter::new(options);
    let text = document(reader, options, &mut meter)?;
    let mut state = ReaderState::new(options, registry)?;
    meter.spend(text.len().saturating_mul(8), 0)?;
    let mut parser = Reader::from_str(&text);
    parser.config_mut().check_end_names = true;
    parser.config_mut().check_comments = true;
    let mut stack: Vec<String> = Vec::new();
    let mut roots = 0usize;
    let mut elements = 0usize;
    loop {
        let event = parser.read_event().map_err(|e| bad(e.to_string()))?;
        let empty = matches!(&event, Event::Empty(_));
        match event {
            Event::Start(start) | Event::Empty(start) => {
                elements = elements.checked_add(1).ok_or_else(|| meter.limit())?;
                meter.cap(elements, options.max_elements)?;
                meter.cap(stack.len().saturating_add(1), options.max_depth)?;
                if stack.is_empty() {
                    roots += 1;
                    if roots != 1 {
                        return Err(bad("multiple XML roots"));
                    }
                }
                let name = std::str::from_utf8(start.name().0)
                    .map_err(|_| bad("invalid XML element name"))?
                    .to_owned();
                let attrs = attributes(&start, &mut meter)?;
                state.start(&name, &attrs, &mut meter)?;
                if empty {
                    state.end(&name, &mut meter)?;
                } else {
                    meter.spend(name.len().saturating_add(32), name.len().saturating_add(64))?;
                    stack.push(name);
                }
            }
            Event::End(end) => {
                let name = std::str::from_utf8(end.name().0)
                    .map_err(|_| bad("invalid XML element name"))?
                    .to_owned();
                if stack.pop().as_deref() != Some(name.as_str()) {
                    return Err(bad("mismatched XML closing tag"));
                }
                state.end(&name, &mut meter)?;
            }
            Event::Text(text) => {
                if stack.is_empty() && !text.as_ref().iter().copied().all(space) {
                    return Err(bad("text outside XML root"));
                }
            }
            Event::CData(_) | Event::GeneralRef(_) => {
                if stack.is_empty() {
                    return Err(bad("character data outside XML root"));
                }
            }
            Event::DocType(_) => {
                return Err(unsupported("DTD and external entities in pepXML"));
            }
            Event::Decl(_) | Event::Comment(_) | Event::PI(_) => {}
            Event::Eof => {
                if roots != 1 || !stack.is_empty() {
                    return Err(bad("incomplete XML document"));
                }
                break;
            }
        }
    }
    if !state.seen_experiment {
        return Err(bad(format!(
            "Found no experiment with name '{}'",
            options.experiment_name
        )));
    }
    let mut document = std::mem::take(&mut state.document);
    // Clean up duplicate ProteinHits in each ProteinIdentification separately,
    // keeping the first occurrence of every accession.
    for run in &mut document.protein_identifications {
        let mut seen = BTreeSet::new();
        run.hits.retain(|hit| seen.insert(hit.accession.clone()));
        run.validate()?;
    }
    for peptide in &document.peptide_identifications {
        peptide.validate()?;
    }
    Ok(document)
}

/// Read a pepXML file with the default options.
///
/// # Errors
///
/// See [`read_with_registry`]; the file is also opened, so a missing file is
/// [`Error::Io`].
pub fn load(path: impl AsRef<Path>) -> Result<PepXmlDocument> {
    load_with_options(path, &ReadOptions::default())
}

/// Read a pepXML file with explicit options.
///
/// # Errors
///
/// See [`load`].
pub fn load_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> Result<PepXmlDocument> {
    load_with_registry(path, options, ModificationsDB::global())
}

/// Read a pepXML file, resolving modifications against `registry`.
///
/// # Errors
///
/// See [`load`].
pub fn load_with_registry(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    registry: &ModificationsDB,
) -> Result<PepXmlDocument> {
    let reader = super::path_io::open(path.as_ref())?;
    read_with_registry(reader, options, registry)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Full-precision double text, as OpenMS `precisionWrapper`.
fn full(value: f64) -> String {
    crate::param::value::format_float(value, true)
}

/// Default `std::ostream` text at the source's stream precision of 15.
///
/// The source sets `f.precision(writtenDigits<double>(0.0))` once and then
/// streams several values without `precisionWrapper`, which is C's `%.15g`:
/// scientific notation when the decimal exponent is below -4 or at least 15,
/// otherwise fixed, with trailing zeros removed in both cases.
fn general(value: f64) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    if value == 0.0 {
        return "0".into();
    }
    // The rounded scientific form carries the decimal exponent exactly, whereas
    // log10 of a value just below a power of ten rounds the wrong way.
    let scientific = format!("{:.*e}", 14, value);
    let (mantissa, exponent) = match scientific.split_once('e') {
        Some((mantissa, exponent)) => (mantissa, exponent.parse::<i32>().unwrap_or(0)),
        None => (scientific.as_str(), 0),
    };
    if !(-4..15).contains(&exponent) {
        let mantissa = trim_zeros(mantissa);
        let sign = if exponent < 0 { "-" } else { "+" };
        return format!("{mantissa}e{sign}{:02}", exponent.unsigned_abs());
    }
    let decimals = usize::try_from(14 - exponent).unwrap_or(0);
    trim_zeros(&format!("{:.*}", decimals, value))
}

fn trim_zeros(value: &str) -> String {
    if !value.contains('.') {
        return value.to_owned();
    }
    let trimmed = value.trim_end_matches('0');
    trimmed.trim_end_matches('.').to_owned()
}

struct Output {
    bytes: Vec<u8>,
    limit: usize,
}
impl Output {
    fn push(&mut self, text: &str) -> Result<()> {
        if self
            .bytes
            .len()
            .checked_add(text.len())
            .is_none_or(|total| total > self.limit)
        {
            return Err(bad("pepXML output byte limit exceeded"));
        }
        self.bytes.extend_from_slice(text.as_bytes());
        Ok(())
    }
    /// Append an attribute value, escaping what an attribute cannot carry.
    ///
    /// `>` is deliberately left alone: it is legal inside an attribute value and
    /// the source writes modification descriptions such as
    /// `Glu->pyro-Glu (N-term E)` raw. The source escapes nothing at all, which
    /// would produce a malformed document for an accession containing `&` or a
    /// description containing `<`.
    fn escaped(&mut self, value: &str) -> Result<()> {
        if !value.chars().all(xml_char) {
            return Err(bad("invalid XML 1.0 character in output"));
        }
        let mut start = 0;
        for (index, c) in value.char_indices() {
            let replacement = match c {
                '&' => "&amp;",
                '<' => "&lt;",
                '"' => "&quot;",
                '\n' => "&#10;",
                '\r' => "&#13;",
                '\t' => "&#9;",
                _ => continue,
            };
            self.push(&value[start..index])?;
            self.push(replacement)?;
            start = index + c.len_utf8();
        }
        self.push(&value[start..])
    }
    fn attribute(&mut self, name: &str, value: &str) -> Result<()> {
        self.push(" ")?;
        self.push(name)?;
        self.push("=\"")?;
        self.escaped(value)?;
        self.push("\"")
    }
}

/// Base name of a path without directories or a known extension.
fn stem_name(path: &str) -> String {
    let file = path.rsplit(['/', '\\']).next().unwrap_or(path);
    strip_extension(file).to_owned()
}

/// Cut sites and the forbidden following residue of a cleavage expression.
///
/// Source `store` splits `DigestionEnzymeProtein::getRegEx()` on `)`, matches
/// `(.*?)([A-Z]+)(.*?)` against the first part for the cut residues and looks
/// for `!P` in the second. The second part is indexed without a bounds check, so
/// an expression containing no `)` - which a `user-defined` enzyme built from a
/// `<specificity>` element always produces - reads past the end of the vector.
/// This returns an empty forbidden residue for that case.
fn cleavage_sites(regex: &str) -> (String, String) {
    let mut parts = regex.split(')');
    let first = parts.next().unwrap_or("");
    let second = parts.next().unwrap_or("");
    let cut: String = first
        .chars()
        .skip_while(|c| !c.is_ascii_uppercase())
        .take_while(char::is_ascii_uppercase)
        .collect();
    let no_cut = if second.contains("!P") { "P" } else { "" };
    (cut, no_cut.to_owned())
}

/// Source `AASequence::toBracketString()` with its default arguments.
///
/// Integer nominal masses, absolute rather than delta, and no fixed-modification
/// suppression list.
fn bracket_string(sequence: &AASequence) -> Result<String> {
    let mut result = String::new();
    if sequence.is_empty() {
        return Ok(result);
    }
    if let Some(modification) = sequence.n_terminal_modification() {
        let mass = modification_delta(modification)? + terminal_gain(true)?;
        result.push_str(&format!("n[{}]", nominal(mass)?));
    }
    for (index, residue) in sequence.as_str().chars().enumerate() {
        result.push(residue);
        if let Some(modification) = sequence.residue_modification(index)? {
            let mass = internal_modified_residue_mass(residue, modification)?;
            result.push_str(&format!("[{}]", nominal(mass)?));
        }
    }
    if let Some(modification) = sequence.c_terminal_modification() {
        let mass = modification_delta(modification)? + terminal_gain(false)?;
        result.push_str(&format!("c[{}]", nominal(mass)?));
    }
    Ok(result)
}

fn nominal(mass: f64) -> Result<i64> {
    let rounded = mass.round();
    if !rounded.is_finite() || rounded.abs() > 9.0e15 {
        return Err(bad("modification mass is not representable"));
    }
    Ok(rounded as i64)
}

fn modification_delta(modification: &SequenceModification) -> Result<f64> {
    modification.diff_mono_mass()
}

fn internal_modified_residue_mass(
    residue: char,
    modification: &SequenceModification,
) -> Result<f64> {
    if let Some(tag) = modification.mass_tag() {
        if let Some(mass) = tag.residue_mono_mass() {
            return Ok(mass);
        }
    }
    Ok(internal_residue_mass(residue)? + modification.diff_mono_mass()?)
}

/// Write a pepXML document with the default options.
///
/// # Errors
///
/// See [`write_with_options`].
pub fn write(writer: impl Write, document: &PepXmlDocument) -> Result<()> {
    write_with_options(writer, document, &WriteOptions::default())
}

/// Write a pepXML document.
///
/// Reproduces `PepXMLFile::store`: one `<msms_run_summary>` with the first
/// protein identification's search parameters, one `<spectrum_query>` per
/// peptide *hit* and one `<search_hit>` inside it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a record fails validation,
/// [`Error::MissingInformation`] when Percolator results carry no posterior
/// error probability - as the source's `Exception::MissingInformation` - and
/// [`Error::Parse`] when an output ceiling is exceeded or a mass is not
/// representable. Writing fails with [`Error::Io`] when the stream fails.
pub fn write_with_options(
    mut writer: impl Write,
    document: &PepXmlDocument,
    options: &WriteOptions,
) -> Result<()> {
    let bytes = render(document, options)?;
    writer.write_all(&bytes)?;
    writer.flush()?;
    Ok(())
}

/// Store a pepXML document at `path` with the default options.
///
/// # Errors
///
/// See [`store_with_options`].
pub fn store(path: impl AsRef<Path>, document: &PepXmlDocument) -> Result<()> {
    store_with_options(path, document, &WriteOptions::default())
}

/// Store a pepXML document at `path`.
///
/// The output is written to a sibling file and published only after it is
/// complete, so a failure never leaves a truncated pepXML behind. The source
/// writes straight into an `ofstream` and can leave a partial file. When
/// [`WriteOptions::output_name`] is empty the output path's stem is used, as the
/// source does. Bytes are written uncompressed even for a `.gz` or `.bz2`
/// filename, because the source's `ofstream` never compresses either.
///
/// # Errors
///
/// A path that cannot be created or written is [`Error::Io`], where the source
/// throws `Exception::UnableToCreateFile`. For everything else see
/// [`write_with_options`], which runs to completion before the file is touched.
pub fn store_with_options(
    path: impl AsRef<Path>,
    document: &PepXmlDocument,
    options: &WriteOptions,
) -> Result<()> {
    let path = path.as_ref();
    let mut options = options.clone();
    if options.output_name.is_empty() {
        options.output_name = stem_name(&path.to_string_lossy());
    }
    let bytes = render(document, &options)?;
    // Plain bytes even when the filename ends in .gz or .bz2: the source opens a
    // bare ofstream and never compresses, exactly as IdXMLFile::store does.
    super::path_io::write_plain(path, |writer| {
        writer.write_all(&bytes)?;
        Ok(())
    })
}

fn render(document: &PepXmlDocument, options: &WriteOptions) -> Result<Vec<u8>> {
    for run in &document.protein_identifications {
        run.validate()?;
    }
    for peptide in &document.peptide_identifications {
        peptide.validate()?;
    }
    let mut out = Output {
        bytes: Vec::new(),
        limit: options.max_output_bytes,
    };
    let mut search_engine = String::new();
    let mut parameters = SearchParameters::default();
    if let Some(first) = document.protein_identifications.first() {
        parameters = first.search_parameters.clone();
        search_engine = match first.search_engine.as_str() {
            "XTandem" => "X! Tandem".into(),
            "Mascot" => "MASCOT".into(),
            // Comet writes "Comet" in its pep.xml, so passing it through is fine.
            other => other.to_owned(),
        };
    }
    let raw_data = if options.mz_file.is_empty() {
        "mzML".to_owned()
    } else {
        super::file_types::type_by_file_name(&options.mz_file)
            .name()
            .to_owned()
    };
    let mut base_name = if options.mz_file.is_empty() {
        options.output_name.clone()
    } else {
        stem_name(&options.mz_file)
    };
    // mz_name is the 'base_name' attribute input from IDFileConverter and is
    // only needed when it differs from mz_file.
    if !options.mz_name.is_empty() {
        base_name = options.mz_name.clone();
    }
    // A spectrum query name is split on dots, so the base name must not contain
    // one or the charge cannot be read back.
    if base_name.contains('.') {
        base_name = base_name.replace('.', "_");
    }
    out.push("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
    out.push(
        "<msms_pipeline_analysis date=\"2007-12-05T17:49:46\" \
         xmlns=\"http://regis-web.systemsbiology.net/pepXML\" \
         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
         xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/\
pepXML_v117.xsd\" summary_xml=\".xml\">\n",
    )?;
    out.push("<msms_run_summary")?;
    out.attribute("base_name", &base_name)?;
    out.push(" raw_data_type=\"raw\"")?;
    out.attribute("raw_data", &format!(".{raw_data}"))?;
    out.attribute("search_engine", &search_engine)?;
    out.push(">\n")?;
    out.push("\t<sample_enzyme")?;
    out.attribute("name", &parameters.digestion_enzyme.to_ascii_lowercase())?;
    out.push(">\n")?;
    out.push("\t\t<specificity cut=\"")?;
    // The source's digestion_enzyme is one DigestionEnzymeProtein carrying both
    // the name and the cleavage expression. This crate stores the name, and a
    // nonempty digestion_regex only for a custom expression, so a named enzyme's
    // expression is read back from the registry.
    let regex = if parameters.digestion_regex.is_empty() {
        find_enzyme(&parameters.digestion_enzyme)
            .map(|enzyme| enzyme.regex().to_owned())
            .unwrap_or_default()
    } else {
        parameters.digestion_regex.clone()
    };
    if !regex.is_empty() {
        let (cut, no_cut) = cleavage_sites(&regex);
        out.escaped(&cut)?;
        if !no_cut.is_empty() {
            out.push("\" no_cut=\"")?;
            out.escaped(&no_cut)?;
        }
    }
    out.push("\" sense=\"C\"/>\n")?;
    out.push("\t</sample_enzyme>\n")?;
    let mass_type = if parameters.mass_type == PeakMassType::Monoisotopic {
        "monoisotopic"
    } else {
        "average"
    };
    out.push("\t<search_summary")?;
    out.attribute("base_name", &base_name)?;
    out.attribute("search_engine", &search_engine)?;
    out.attribute("precursor_mass_type", mass_type)?;
    out.attribute("fragment_mass_type", mass_type)?;
    out.push(" out_data_type=\"\" out_data=\"\" search_id=\"1\">\n")?;
    out.push("\t\t<search_database")?;
    out.attribute("local_path", &parameters.database)?;
    out.push(" type=\"AA\"/>\n")?;
    write_modification_declarations(&mut out, document)?;
    out.push("\t</search_summary>\n")?;
    if options.peptideprophet_analyzed {
        out.push(
            "\t<analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" \
             id=\"1\"/>\n",
        )?;
    }
    write_queries(
        &mut out,
        document,
        options,
        &base_name,
        &search_engine,
        &parameters,
    )?;
    out.push("</msms_run_summary>\n")?;
    out.push("</msms_pipeline_analysis>\n")?;
    Ok(out.bytes)
}

fn write_modification_declarations(out: &mut Output, document: &PepXmlDocument) -> Result<()> {
    // Only the FIRST hit of every identification is inspected, as in the source,
    // and each set is keyed by the modification's full identifier, so the
    // declarations come out in the source's `std::set<std::string>` order.
    let mut residue_modifications: BTreeMap<String, (char, SequenceModification)> = BTreeMap::new();
    let mut n_terminal: BTreeMap<String, SequenceModification> = BTreeMap::new();
    let mut c_terminal: BTreeMap<String, SequenceModification> = BTreeMap::new();
    for peptide in &document.peptide_identifications {
        let Some(hit) = peptide.hits.first() else {
            continue;
        };
        if !hit.sequence.is_modified() {
            continue;
        }
        if let Some(modification) = hit.sequence.n_terminal_modification() {
            n_terminal.insert(modification.full_id().to_owned(), modification.clone());
        }
        if let Some(modification) = hit.sequence.c_terminal_modification() {
            c_terminal.insert(modification.full_id().to_owned(), modification.clone());
        }
        for (index, residue) in hit.sequence.as_str().chars().enumerate() {
            if let Some(modification) = hit.sequence.residue_modification(index)? {
                residue_modifications.insert(
                    modification.full_id().to_owned(),
                    (residue, modification.clone()),
                );
            }
        }
    }
    for (full_id, (residue, modification)) in &residue_modifications {
        let origin = modification.origin().unwrap_or(*residue);
        let difference = modification.diff_mono_mass()?;
        let mass = internal_modified_residue_mass(*residue, modification)?;
        out.push("\t\t<aminoacid_modification")?;
        out.attribute("aminoacid", &origin.to_string())?;
        out.attribute("massdiff", &full(difference))?;
        out.attribute("mass", &full(mass))?;
        out.push(" variable=\"Y\" binary=\"N\"")?;
        out.attribute("description", full_id)?;
        out.push("/>\n")?;
    }
    for (terminus, modifications) in [("n", &n_terminal), ("c", &c_terminal)] {
        for (full_id, modification) in modifications.iter() {
            let difference = modification.diff_mono_mass()?;
            // Source writes ResidueModification::getMonoMass(), which is zero for
            // every UniMod record because only the mass difference is tabulated.
            let mass = modification
                .known()
                .map_or(0.0, ResidueModification::mono_mass);
            out.push("\t\t<terminal_modification")?;
            out.attribute("terminus", terminus)?;
            out.attribute("massdiff", &full(difference))?;
            out.attribute("mass", &full(mass))?;
            out.push(" variable=\"Y\"")?;
            out.attribute("description", full_id)?;
            out.push(" protein_terminus=\"\"/>\n")?;
        }
    }
    Ok(())
}

fn write_queries(
    out: &mut Output,
    document: &PepXmlDocument,
    options: &WriteOptions,
    base_name: &str,
    search_engine: &str,
    parameters: &SearchParameters,
) -> Result<()> {
    // The scan index is zero-based and the scan number one-based; both are
    // reconstructed from the identification order when no lookup is available.
    let mut count = 0i64;
    let mut written = 0usize;
    for peptide in &document.peptide_identifications {
        if peptide.hits.is_empty() {
            count += 1;
            continue;
        }
        for hit in &peptide.hits {
            written += 1;
            if written > options.max_identifications {
                return Err(bad("pepXML output identification limit exceeded"));
            }
            let precursor_neutral_mass = hit.sequence.mono_mass()?;
            let mut scan_index = count;
            let scan_number;
            if options.lookup.is_empty() {
                // XTandemXMLFile sets "RT_index" for X!Tandem results.
                if let Some(value) = peptide.metadata.get("RT_index") {
                    scan_index = value.as_i64()?;
                }
                scan_number = scan_index
                    .checked_add(1)
                    .ok_or_else(|| bad("scan index overflows"))?;
            } else {
                let reference = peptide.spectrum_reference();
                let index = if peptide.metadata.contains_key("spectrum_reference") {
                    options.lookup.find_by_native_id(&reference)
                } else {
                    peptide.rt.and_then(|rt| options.lookup.find_by_rt(rt))
                };
                let index = index.ok_or_else(|| {
                    Error::InvalidValue(
                        "no spectrum matches this identification's reference or retention time"
                            .into(),
                    )
                })?;
                let meta = options
                    .lookup
                    .get(index)
                    .ok_or_else(|| Error::InvalidValue("spectrum index out of range".into()))?;
                scan_index = i64::try_from(index)
                    .map_err(|_| Error::InvalidValue("spectrum index overflows".into()))?;
                scan_number = meta
                    .scan_number
                    .and_then(|value| i64::try_from(value).ok())
                    .unwrap_or(-1);
            }
            // PeptideProphet requires this exact "spectrum" format or the TPP
            // reports a parsing error. iProphet's InterProphetParser is the
            // reference for the strictly required "spectrum" and "assumed_charge"
            // attributes and the optional "retention_time_sec", "swath_assay" and
            // "experiment_label".
            let mut spectrum_name = format!("{base_name}.{scan_number}.{scan_number}.");
            if options.keep_native_spectrum_name {
                if let Some(value) = peptide.metadata.get("pepxml_spectrum_name") {
                    spectrum_name = value.to_string();
                }
            }
            out.push("\t<spectrum_query")?;
            out.attribute("spectrum", &format!("{spectrum_name}{}", hit.charge))?;
            out.attribute("start_scan", &scan_number.to_string())?;
            out.attribute("end_scan", &scan_number.to_string())?;
            out.attribute("precursor_neutral_mass", &full(precursor_neutral_mass))?;
            out.attribute("assumed_charge", &hit.charge.to_string())?;
            out.attribute("index", &scan_index.to_string())?;
            if let Some(rt) = peptide.rt {
                out.attribute("retention_time_sec", &general(rt))?;
                out.push(" ")?;
            }
            let label = peptide.experiment_label();
            if !label.is_empty() {
                out.attribute("experiment_label", &label)?;
                out.push(" ")?;
            }
            // "swath_assay" is an optional SWATH-MS parameter; the TPP parses it
            // as "xxx:yyy" where yyy indexes the SWATH window.
            if let Some(value) = peptide.metadata.get("swath_assay") {
                out.attribute("swath_assay", &value.to_string())?;
                out.push(" ")?;
            }
            if let Some(value) = peptide.metadata.get("status") {
                out.attribute("status", &value.to_string())?;
                out.push(" ")?;
            }
            out.push(">\n")?;
            out.push("\t<search_result>\n")?;
            write_hit(
                out,
                peptide,
                hit,
                search_engine,
                parameters,
                options,
                precursor_neutral_mass,
            )?;
            out.push("\t</search_result>\n")?;
            out.push("\t</spectrum_query>\n")?;
        }
        count += 1;
    }
    Ok(())
}

fn write_hit(
    out: &mut Output,
    peptide: &PeptideIdentification,
    hit: &PeptideHit,
    search_engine: &str,
    parameters: &SearchParameters,
    options: &WriteOptions,
    precursor_neutral_mass: f64,
) -> Result<()> {
    // The first evidence is written as the leader and the rest as alternatives.
    let leader = hit.evidences.first().cloned().unwrap_or_default();
    let rank = hit
        .rank
        .checked_add(1)
        .ok_or_else(|| bad("hit rank overflows"))?;
    out.push("\t\t<search_hit")?;
    out.attribute("hit_rank", &rank.to_string())?;
    out.attribute("peptide", hit.sequence.as_str())?;
    out.attribute(
        "peptide_prev_aa",
        &flanking_code(leader.aa_before).to_string(),
    )?;
    out.attribute(
        "peptide_next_aa",
        &flanking_code(leader.aa_after).to_string(),
    )?;
    out.attribute("protein", &leader.protein_accession)?;
    out.push(" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\"")?;
    out.attribute("calc_neutral_pep_mass", &full(precursor_neutral_mass))?;
    out.push(" massdiff=\"0.0\"")?;
    let tryptic = matches!(leader.aa_before, FlankingResidue::Residue('R' | 'K'))
        && parameters.digestion_enzyme == "Trypsin";
    let num_tol_term = if tryptic { 2 } else { 1 };
    out.attribute("num_tol_term", &num_tol_term.to_string())?;
    out.push(" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">\n")?;
    for evidence in hit.evidences.iter().skip(1) {
        out.push("\t\t<alternative_protein")?;
        out.attribute("protein", &evidence.protein_accession)?;
        out.attribute("num_tol_term", &num_tol_term.to_string())?;
        if evidence.aa_before != FlankingResidue::Unknown {
            out.attribute(
                "peptide_prev_aa",
                &flanking_code(evidence.aa_before).to_string(),
            )?;
        }
        if evidence.aa_after != FlankingResidue::Unknown {
            out.attribute(
                "peptide_next_aa",
                &flanking_code(evidence.aa_after).to_string(),
            )?;
        }
        out.push("/>\n")?;
    }
    if hit.sequence.is_modified() {
        out.push("\t\t\t<modification_info")?;
        out.attribute("modified_peptide", &bracket_string(&hit.sequence)?)?;
        if let Some(modification) = hit.sequence.n_terminal_modification() {
            let mass = terminal_gain(true)? + modification.diff_mono_mass()?;
            out.attribute("mod_nterm_mass", &full(mass))?;
        }
        if let Some(modification) = hit.sequence.c_terminal_modification() {
            let mass = terminal_gain(false)? + modification.diff_mono_mass()?;
            out.attribute("mod_cterm_mass", &full(mass))?;
        }
        out.push(">\n")?;
        for (index, residue) in hit.sequence.as_str().chars().enumerate() {
            if let Some(modification) = hit.sequence.residue_modification(index)? {
                // Positions are 1-based. The source adds the modification's
                // absolute mono mass - zero for a UniMod record - to the internal
                // residue mass, which already includes the modification.
                let mass = modification
                    .known()
                    .map_or(0.0, ResidueModification::mono_mass)
                    + internal_modified_residue_mass(residue, modification)?;
                out.push("\t\t\t\t<mod_aminoacid_mass")?;
                out.attribute("position", &(index + 1).to_string())?;
                out.attribute("mass", &full(mass))?;
                out.push("/>\n")?;
            }
        }
        out.push("\t\t\t</modification_info>\n")?;
    }
    let mut peptideprophet_written = false;
    for result in &hit.analysis_results {
        out.push("\t\t\t<analysis_result")?;
        out.attribute("analysis", &result.score_type)?;
        out.push(">\n")?;
        let tag = match result.score_type.as_str() {
            "peptideprophet" => {
                peptideprophet_written = true;
                "peptideprophet_result"
            }
            "interprophet" => "interprophet_result",
            _ => {
                peptideprophet_written = true;
                "peptideprophet_result"
            }
        };
        let score = general(result.main_score);
        out.push(&format!("\t\t\t\t<{tag} probability=\"{score}\""))?;
        out.push(&format!(" all_ntt_prob=\"({score},{score},{score})\">\n"))?;
        if !result.sub_scores.is_empty() {
            out.push("\t\t\t\t\t<search_score_summary>\n")?;
            for (name, value) in &result.sub_scores {
                out.push("\t\t\t\t\t\t<parameter")?;
                out.attribute("name", name)?;
                out.attribute("value", &general(*value))?;
                out.push("/>\n")?;
            }
            out.push("\t\t\t\t\t</search_score_summary>\n")?;
        }
        out.push(&format!("\t\t\t\t</{tag}>\n"))?;
        out.push("\t\t\t</analysis_result>\n")?;
    }
    if options.peptideprophet_analyzed && !peptideprophet_written {
        // Deprecated way of writing PeptideProphet results, used only when
        // requested explicitly and when none were written from AnalysisResults.
        let score = general(hit.score);
        out.push("\t\t\t<analysis_result analysis=\"peptideprophet\">\n")?;
        out.push(&format!(
            "\t\t\t<peptideprophet_result probability=\"{score}\" \
             all_ntt_prob=\"({score},{score},{score})\">\n"
        ))?;
        out.push("\t\t\t</peptideprophet_result>\n")?;
        out.push("\t\t\t</analysis_result>\n")?;
    } else {
        write_scores(out, peptide, hit, search_engine)?;
    }
    out.push("\t\t</search_hit>\n")?;
    Ok(())
}

fn meta_text(hit: &PeptideHit, key: &str) -> Result<String> {
    hit.metadata
        .get(key)
        .map(ToString::to_string)
        .ok_or_else(|| Error::MissingInformation(format!("peptide hit has no '{key}' meta value")))
}

fn write_scores(
    out: &mut Output,
    peptide: &PeptideIdentification,
    hit: &PeptideHit,
    search_engine: &str,
) -> Result<()> {
    let has_pep =
        peptide.score_type == "Posterior Error Probability" || peptide.score_type == "pep";
    let mut percolator = false;
    let score = |name: &str, value: &str| {
        format!("\t\t\t<search_score name=\"{name}\" value=\"{value}\"/>\n")
    };
    match search_engine {
        "X! Tandem" => {
            if peptide.score_type == "XTandem" {
                out.push(&score("hyperscore", &general(hit.score)))?;
                let next = match hit.metadata.get("nextscore") {
                    Some(value) => value.to_string(),
                    None => general(hit.score),
                };
                out.push(&score("nextscore", &next))?;
            } else if let Some(value) = hit.metadata.get("XTandem_score") {
                let value = value.to_string();
                out.push(&score("hyperscore", &value))?;
                let next = match hit.metadata.get("nextscore") {
                    Some(next) => next.to_string(),
                    None => value,
                };
                out.push(&score("nextscore", &next))?;
            }
            out.push(&score("expect", &meta_text(hit, "E-Value")?))?;
        }
        "Comet" => {
            for (name, key) in [
                ("xcorr", "MS:1002252"),
                ("deltacn", "MS:1002253"),
                ("deltacnstar", "MS:1002254"),
                ("spscore", "MS:1002255"),
                ("sprank", "MS:1002256"),
                ("expect", "MS:1002257"),
            ] {
                out.push(&score(name, &meta_text(hit, key)?))?;
            }
        }
        "MASCOT" => {
            out.push(&score("expect", &meta_text(hit, "EValue")?))?;
            out.push(&score("ionscore", &general(hit.score)))?;
        }
        "OMSSA" | "MSGFPlus" => {
            out.push(&score("expect", &general(hit.score)))?;
        }
        "Percolator" => {
            for (name, keys) in [
                ("Percolator_score", ["MS:1001492", "Percolator_score"]),
                ("Percolator_qvalue", ["MS:1001491", "Percolator_qvalue"]),
            ] {
                if let Some(value) = keys.iter().find_map(|key| hit.metadata.get(*key)) {
                    out.push(&score(name, &general(value.as_f64()?)))?;
                }
            }
            let mut error_probability = 0.0;
            match ["MS:1001493", "Percolator_PEP"]
                .iter()
                .find_map(|key| hit.metadata.get(*key))
            {
                Some(value) => error_probability = value.as_f64()?,
                None if !has_pep => {
                    return Err(Error::MissingInformation(
                        "Percolator PEP score missing for pepXML export of Percolator results."
                            .into(),
                    ));
                }
                // The posterior error probability is written below instead.
                None => {}
            }
            out.push(&score("Percolator_PEP", &general(error_probability)))?;
            let probability = general(1.0 - error_probability);
            out.push("\t\t\t<analysis_result analysis=\"peptideprophet\">\n")?;
            out.push(&format!(
                "\t\t\t\t<peptideprophet_result probability=\"{probability}\" \
                 all_ntt_prob=\"(0.0000,0.0000,{probability})\"/>\n"
            ))?;
            out.push("\t\t\t</analysis_result>\n")?;
            percolator = true;
        }
        _ => {
            out.push(&score(&peptide.score_type, &general(hit.score)))?;
        }
    }
    // Any search engine with a posterior error probability, including OpenMS's
    // own IDPosteriorErrorProbability, except Percolator, which is done above.
    if has_pep && !percolator {
        out.push(&score(&peptide.score_type, &general(hit.score)))?;
        let probability = general(1.0 - hit.score);
        out.push("\t\t\t<analysis_result analysis=\"peptideprophet\">\n")?;
        out.push(&format!(
            "\t\t\t\t<peptideprophet_result probability=\"{probability}\" \
             all_ntt_prob=\"(0.0000,0.0000,{probability})\"/>\n"
        ))?;
        out.push("\t\t\t</analysis_result>\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general_matches_stream_precision_fifteen() {
        assert_eq!(general(8.0), "8");
        assert_eq!(general(7.7), "7.7");
        assert_eq!(general(0.00091), "0.00091");
        assert_eq!(general(7.9e-05), "7.9e-05");
        assert_eq!(general(1.3653), "1.3653");
        assert_eq!(general(0.0), "0");
        assert_eq!(general(-0.027), "-0.027");
    }

    #[test]
    fn diff_mass_strings_carry_an_explicit_sign() {
        assert_eq!(diff_mono_mass_string(1.0), "+1.0");
        assert_eq!(diff_mono_mass_string(2.5), "+2.5");
        assert_eq!(diff_mono_mass_string(-2.5), "-2.5");
        assert_eq!(diff_mono_mass_string(3.4), "+3.4");
    }

    #[test]
    fn scan_numbers_come_from_a_trailing_assignment() {
        assert_eq!(extract_scan_number("scan=12"), Some(12));
        assert_eq!(
            extract_scan_number("controllerType=0 controllerNumber=1 scan=7"),
            Some(7)
        );
        // The source regexp is anchored at the end of the identifier, so a
        // trailing token wins even when an earlier one is the real scan number.
        assert_eq!(extract_scan_number("scan=12 merged=1"), Some(1));
        assert_eq!(extract_scan_number("frame=5 scan=7 precursor=3"), Some(3));
        assert_eq!(extract_scan_number("spectrum"), None);
        assert_eq!(extract_scan_number("="), None);
        assert_eq!(extract_scan_number("scan12"), None);
        assert_eq!(extract_scan_number(""), None);
        assert_eq!(extract_scan_number("日本語=4"), Some(4));
    }

    #[test]
    fn cleavage_sites_extract_cut_and_forbidden_residues() {
        assert_eq!(
            cleavage_sites("(?<=[KRX])(?!P)"),
            ("KRX".into(), "P".into())
        );
        assert_eq!(cleavage_sites("(?<=[RX])"), ("RX".into(), String::new()));
        // A user-defined enzyme's expression carries no ')' at all; the source
        // reads past its split vector here.
        assert_eq!(cleavage_sites("KR"), ("KR".into(), String::new()));
        assert_eq!(cleavage_sites("()"), (String::new(), String::new()));
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Minimal in-memory peptide-spectrum search engine, ported from
//! `ANALYSIS/ID/SimpleSearchEngineAlgorithm.h`.
//!
//! Searches MS2 spectra against a protein FASTA database and produces protein-
//! and peptide-level identifications. It is a self-contained reference and
//! teaching implementation, not a feature-complete replacement for external
//! engines such as MSGF+, Comet, Sage or MSFragger.
//!
//! The pipeline run by [`SimpleSearchEngineAlgorithm::search`](crate::analysis::id::simple_search_engine::SimpleSearchEngineAlgorithm::search):
//!
//! 1. Load spectra (MS level 2 only) and the protein database.
//! 2. Preprocess MS2 spectra (filtering, deisotoping, decharging) with
//!    [`SimpleSearchEngineAlgorithm::preprocess_spectra`](crate::analysis::id::simple_search_engine::SimpleSearchEngineAlgorithm::preprocess_spectra), using the
//!    configured fragment tolerance.
//! 3. Digest each protein with the configured enzyme, specificity and missed
//!    cleavages, generate every fixed/variable modification variant (capped at
//!    `modifications:variable_max_per_peptide`), and score them against MS2
//!    spectra whose precursor mass matches the candidate within
//!    `precursor:mass_tolerance`, with the isotope corrections of
//!    `precursor:isotopes`.
//! 4. Keep the top-N candidates per spectrum (`report:top_hits`), optionally
//!    with target/decoy results (`decoys`), FDR filtering (`FDR:PSM`, q-value)
//!    and PSM annotations (`annotate:PSM`), in
//!    [`SimpleSearchEngineAlgorithm::post_process_hits`](crate::analysis::id::simple_search_engine::SimpleSearchEngineAlgorithm::post_process_hits).
//!
//! See `docs/SIMPLE_SEARCH_ENGINE_SUPPORT.md` for the API mapping and evidence.

use crate::analysis::false_discovery_rate::FalseDiscoveryRate;
use crate::analysis::id_filter;
use crate::analysis::peptide_indexing::{DecoyRule, MissingDecoyAction, PeptideIndexing};
use crate::analysis::psm_scoring::HyperScore;
use crate::chemistry::theoretical::TheoreticalSpectrumGenerator;
use crate::chemistry::{
    AASequence, C13C12_MASSDIFF_U, DecoyGenerator, DigestionSpecificity, ModificationsDB,
    ModifiedPeptideGenerator, PROTON_MASS_U, Protease, ProteaseDB, ProteaseDigestion,
    ResidueModification,
};
use crate::comparison::{SpectrumAlignment, Tolerance};
use crate::concept::boost_regex::BoostRegex;
use crate::concept::constants::user_param;
use crate::concept::math_functions::RandomShuffler;
use crate::concept::math_functions::ppm;
use crate::concept::parallel::{Threads, map_collect};
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::data_structures::datetime::DateTime;
use crate::format::fasta::{FASTAEntry, FASTAFile};
use crate::format::file_types::type_by_file_name;
use crate::format::{FileHandler, FileType, PeakFileOptions};
use crate::identification::{
    EnzymeTermSpecificity, PeakAnnotation, PeakMassType, PeptideHit, PeptideIdentification,
    ProteinIdentification, SearchParameters,
};
use crate::kernel::MSExperiment;
use crate::math::statistic_functions::median;
use crate::metadata::{ImTypes, IonMobilityFormat, MetaValue};
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::processing::deisotoping::Deisotoper;
use crate::processing::window_mower::{WindowMower, WindowMowerMethod};
use crate::processing::{NLargest, Normalizer, SpectrumFilter, ThresholdMower};
use crate::{Error, MSSpectrum, Result};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::Arc;

/// The source `DefaultParamHandler` name.
pub const SIMPLE_SEARCH_ENGINE_ALGORITHM_NAME: &str = "SimpleSearchEngineAlgorithm";

/// Unmodified peptides scored per parallel batch before their hits are merged.
/// Bounds the hits held in memory at once; the result does not depend on it.
const PEPTIDE_BATCH: usize = 4096;

/// Unit of a mass tolerance, as the source's `"ppm"`/`"Da"` strings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MassToleranceUnit {
    /// Parts per million of the reference mass (source `"ppm"`).
    #[default]
    Ppm,
    /// Dalton (source `"Da"`).
    Da,
}

impl MassToleranceUnit {
    /// The source parameter string.
    pub fn name(self) -> &'static str {
        match self {
            Self::Ppm => "ppm",
            Self::Da => "Da",
        }
    }

    /// Parse the source parameter string.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for anything but `ppm` or `Da`.
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "ppm" => Ok(Self::Ppm),
            "Da" => Ok(Self::Da),
            _ => Err(bad(format!(
                "mass tolerance unit must be ppm or Da, not {name:?}"
            ))),
        }
    }

    fn tolerance(self, value: f64) -> Tolerance {
        match self {
            Self::Ppm => Tolerance::Ppm(value),
            Self::Da => Tolerance::Absolute(value),
        }
    }
}

/// Outcome of [`SimpleSearchEngineAlgorithm::search`], distinguishing
/// recoverable input issues from execution errors (source `ExitCodes`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExitCodes {
    /// Search completed; the outputs contain the result.
    ExecutionOk,
    /// The input contained no usable data; returned for an empty protein
    /// database, which the source's peptide indexer reports as
    /// `DATABASE_EMPTY`.
    InputFileEmpty,
    /// Internal post-condition violated. Not returned by this port: indexing
    /// failures are errors (see [`SimpleSearchEngineAlgorithm::search`]).
    UnexpectedResult,
    /// Caught a generic exception. Not returned by this port, which reports
    /// failures as `Err`.
    UnknownError,
    /// Inconsistent configuration. Not returned by the source either; this
    /// port reports invalid parameters as `Err`.
    IllegalParameters,
}

/// Compact record of one scored peptide candidate against one spectrum
/// (source `AnnotatedHit_`).
///
/// Used instead of [`PeptideHit`] while scoring so that all candidates per
/// spectrum fit in memory; only the top-N entries per spectrum are promoted to
/// [`PeptideHit`] by [`SimpleSearchEngineAlgorithm::post_process_hits`]. The
/// source's field-order padding optimisation has no Rust meaning.
#[derive(Clone, Debug, PartialEq)]
pub struct AnnotatedHit {
    /// Unmodified peptide sequence (the source keeps a view into the digested
    /// protein; this shares one allocation per distinct peptide).
    pub sequence: Arc<str>,
    /// Index of the modification variant, for re-materialisation.
    pub peptide_mod_index: usize,
    /// Main score (higher is better): ln(hyperscore).
    pub score: f64,
    /// Fraction of annotated prefix ions (a/b/c).
    pub prefix_fraction: f32,
    /// Fraction of annotated suffix ions (x/y/z).
    pub suffix_fraction: f32,
    /// Mean absolute fragment mass error (Da, or ppm for a ppm tolerance).
    pub mean_error: f32,
    /// Matched precursor isotope offset (per `precursor:isotopes`).
    pub isotope_error: i32,
    /// Number of matched prefix ions (a/b/c), saturated to `u16` as the source
    /// casts it.
    pub matched_prefix_ions: u16,
    /// Number of matched suffix ions (x/y/z), saturated to `u16`.
    pub matched_suffix_ions: u16,
}

impl AnnotatedHit {
    /// Strict weak ordering of the source's `hasBetterScore`: higher score
    /// first, then lower modification index, then shorter sequence, then the
    /// lexicographically smaller sequence ("shortlex", which the source keeps
    /// deliberately because its idXML output depends on it).
    pub fn has_better_score(a: &Self, b: &Self) -> bool {
        Self::better_first(a, b) == Ordering::Less
    }

    /// [`AnnotatedHit::has_better_score`] as an ordering, best first. Scores
    /// are finite (hyperscores are validated), so `total_cmp` agrees with the
    /// source's `!=`/`>` comparison.
    pub fn better_first(a: &Self, b: &Self) -> Ordering {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.peptide_mod_index.cmp(&b.peptide_mod_index))
            .then_with(|| a.sequence.len().cmp(&b.sequence.len()))
            .then_with(|| a.sequence.cmp(&b.sequence))
    }
}

/// Minimal in-memory peptide-spectrum search engine (source
/// `SimpleSearchEngineAlgorithm`).
///
/// The fields are the source's members, set from parameters by
/// [`SimpleSearchEngineAlgorithm::from_param`] exactly as its
/// `updateMembers_` does; [`SimpleSearchEngineAlgorithm::defaults`] lists
/// them with their descriptions. The source's `ProgressLogger` base is not
/// ported: progress is not reported.
#[derive(Clone, Debug, PartialEq)]
pub struct SimpleSearchEngineAlgorithm {
    /// Precursor mass tolerance value (`precursor:mass_tolerance`, 10).
    pub precursor_mass_tolerance: f64,
    /// Its unit (`precursor:mass_tolerance_unit`, ppm).
    pub precursor_mass_tolerance_unit: MassToleranceUnit,
    /// Minimum precursor charge considered (`precursor:min_charge`, 2).
    pub precursor_min_charge: i32,
    /// Maximum precursor charge considered (`precursor:max_charge`, 5).
    pub precursor_max_charge: i32,
    /// Allowed precursor isotope offsets (`precursor:isotopes`, `[0, 1]`):
    /// 0 is the monoisotopic peak, 1 corrects a misassignment to the first
    /// isotopic peak.
    pub precursor_isotopes: Vec<i32>,
    /// Fragment mass tolerance value (`fragment:mass_tolerance`, 10).
    pub fragment_mass_tolerance: f64,
    /// Its unit (`fragment:mass_tolerance_unit`, ppm).
    pub fragment_mass_tolerance_unit: MassToleranceUnit,
    /// UniMod names of fixed modifications (`modifications:fixed`,
    /// `["Carbamidomethyl (C)"]`).
    pub modifications_fixed: Vec<String>,
    /// UniMod names of variable modifications (`modifications:variable`,
    /// `["Oxidation (M)"]`).
    pub modifications_variable: Vec<String>,
    /// Maximum residues carrying a variable modification per candidate
    /// (`modifications:variable_max_per_peptide`, 2).
    pub modifications_max_variable_mods_per_peptide: usize,
    /// Enzyme name as [`ProteaseDB`] knows it (`enzyme`, Trypsin).
    pub enzyme: String,
    /// Generate target/decoy results (`decoys`, false).
    pub decoys: bool,
    /// q-value threshold for PSM filtering (`FDR:PSM`, 0 = disabled); needs
    /// [`SimpleSearchEngineAlgorithm::decoys`].
    pub fdr_psm: f64,
    /// PSM meta-value annotations to add (`annotate:PSM`, `["ALL"]`).
    pub annotate_psm: Vec<String>,
    /// Minimum peptide length after digestion (`peptide:min_size`, 7).
    pub peptide_min_size: usize,
    /// Maximum peptide length after digestion, 0 = unlimited
    /// (`peptide:max_size`, 40).
    pub peptide_max_size: usize,
    /// Allowed missed cleavages (`peptide:missed_cleavages`, 1).
    pub peptide_missed_cleavages: usize,
    /// Enzyme specificity at both termini (`peptide:enzyme_specificity`,
    /// full).
    pub peptide_enzyme_specificity: DigestionSpecificity,
    /// Optional regular expression; only peptides matching it in full are
    /// searched (`peptide:motif`, empty).
    pub peptide_motif: String,
    /// Top-scoring PSMs reported per spectrum (`report:top_hits`, 1).
    pub report_top_hits: usize,
    /// Worker threads for preprocessing and scoring (native; the source uses
    /// OpenMP's setting). Results are identical for every value.
    pub threads: Threads,
}

impl Default for SimpleSearchEngineAlgorithm {
    fn default() -> Self {
        Self {
            precursor_mass_tolerance: 10.0,
            precursor_mass_tolerance_unit: MassToleranceUnit::Ppm,
            precursor_min_charge: 2,
            precursor_max_charge: 5,
            precursor_isotopes: vec![0, 1],
            fragment_mass_tolerance: 10.0,
            fragment_mass_tolerance_unit: MassToleranceUnit::Ppm,
            modifications_fixed: vec!["Carbamidomethyl (C)".into()],
            modifications_variable: vec!["Oxidation (M)".into()],
            modifications_max_variable_mods_per_peptide: 2,
            enzyme: "Trypsin".into(),
            decoys: false,
            fdr_psm: 0.0,
            annotate_psm: vec!["ALL".into()],
            peptide_min_size: 7,
            peptide_max_size: 40,
            peptide_missed_cleavages: 1,
            peptide_enzyme_specificity: DigestionSpecificity::Full,
            peptide_motif: String::new(),
            report_top_hits: 1,
            threads: Threads::default(),
        }
    }
}

/// The `annotate:PSM` values the source accepts, in its order.
const ANNOTATIONS: [&str; 11] = [
    "ALL",
    user_param::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
    user_param::PRECURSOR_ERROR_PPM_USERPARAM,
    user_param::MATCHED_PREFIX_IONS_FRACTION,
    user_param::MATCHED_SUFFIX_IONS_FRACTION,
    user_param::NUM_MATCHED_PEAKS,
    user_param::MATCHED_PREFIX_IONS,
    user_param::MATCHED_SUFFIX_IONS,
    user_param::LONGEST_PEPTIDE_ION_SEQUENCE,
    user_param::MATCHED_ION_CURRENT,
    user_param::FRAGMENT_ANNOTATION_USERPARAM,
];

const SPECIFICITY_DESCRIPTION: &str = "Enzyme cleavage specificity required for both peptide termini.\n  'full' : both termini must be enzyme-specific (canonical, e.g. tryptic).\n  'semi' : only one terminus needs to be enzyme-specific (semi-tryptic).\n  'none' : no enzyme constraint at either terminus; every substring of length\n           [min_size, max_size] is enumerated. Use this for immunopeptidomics\n           (e.g. HLA peptides 8..12mers). For very large search spaces consider\n           tightening 'peptide:min_size'/'peptide:max_size'.";

impl SimpleSearchEngineAlgorithm {
    /// The source defaults with their descriptions and restrictions, as its
    /// constructor installs them.
    ///
    /// The modification lists are restricted to every UniMod-backed
    /// modification in [`ModificationsDB::global`] and `enzyme` to every
    /// [`ProteaseDB`] name, as the source's `getAllSearchModifications` and
    /// `getAllNames` restrict them. The source's description of
    /// `fragment:mass_tolerance_unit` is the truncated `"Unit of fragment m"`;
    /// it is kept verbatim.
    ///
    /// # Errors
    ///
    /// Propagates [`Param`] construction errors, which the fixed defaults do
    /// not trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this engine in the layout of
    /// [`SimpleSearchEngineAlgorithm::defaults`], carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a count does not fit the source's
    /// 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        let tags: &[String] = &[];
        let units = ["ppm".to_owned(), "Da".to_owned()];
        p.set_value(
            "precursor:mass_tolerance",
            ParamValue::Float(self.precursor_mass_tolerance),
            "+/- tolerance for precursor mass.",
            tags,
        )?;
        p.set_value(
            "precursor:mass_tolerance_unit",
            ParamValue::String(self.precursor_mass_tolerance_unit.name().into()),
            "Unit of precursor mass tolerance +/- to the theoretical one.",
            tags,
        )?;
        p.set_valid_strings("precursor:mass_tolerance_unit", &units)?;
        p.set_value(
            "precursor:min_charge",
            ParamValue::Integer(i64::from(self.precursor_min_charge)),
            "Minimum precursor charge to be considered.",
            tags,
        )?;
        p.set_value(
            "precursor:max_charge",
            ParamValue::Integer(i64::from(self.precursor_max_charge)),
            "Maximum precursor charge to be considered.",
            tags,
        )?;
        p.set_section_description("precursor", "Precursor (Parent Ion) Options")?;
        p.set_value(
            "precursor:isotopes",
            ParamValue::IntegerList(self.precursor_isotopes.clone()),
            "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)",
            tags,
        )?;
        p.set_value(
            "fragment:mass_tolerance",
            ParamValue::Float(self.fragment_mass_tolerance),
            "Fragment mass tolerance +/- to the theoretical one",
            tags,
        )?;
        p.set_value(
            "fragment:mass_tolerance_unit",
            ParamValue::String(self.fragment_mass_tolerance_unit.name().into()),
            "Unit of fragment m",
            tags,
        )?;
        p.set_valid_strings("fragment:mass_tolerance_unit", &units)?;
        p.set_section_description("fragment", "Fragments (Product Ion) Options")?;
        let all_mods = ModificationsDB::global().all_search_modifications();
        p.set_value(
            "modifications:fixed",
            ParamValue::StringList(self.modifications_fixed.clone()),
            "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'",
            tags,
        )?;
        p.set_valid_strings("modifications:fixed", &all_mods)?;
        p.set_value(
            "modifications:variable",
            ParamValue::StringList(self.modifications_variable.clone()),
            "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'",
            tags,
        )?;
        p.set_valid_strings("modifications:variable", &all_mods)?;
        p.set_value(
            "modifications:variable_max_per_peptide",
            int_value(
                self.modifications_max_variable_mods_per_peptide,
                "variable_max_per_peptide",
            )?,
            "Maximum number of residues carrying a variable modification per candidate peptide",
            tags,
        )?;
        p.set_section_description("modifications", "Modifications Options")?;
        let enzymes: Vec<String> = ProteaseDB::global()
            .names()
            .into_iter()
            .map(str::to_owned)
            .collect();
        p.set_value(
            "enzyme",
            ParamValue::String(self.enzyme.clone()),
            "The enzyme used for peptide digestion.",
            tags,
        )?;
        p.set_valid_strings("enzyme", &enzymes)?;
        p.set_value(
            "decoys",
            ParamValue::String(if self.decoys { "true" } else { "false" }.into()),
            "Should decoys be generated?",
            tags,
        )?;
        p.set_valid_strings("decoys", &["true".into(), "false".into()])?;
        p.set_value(
            "annotate:PSM",
            ParamValue::StringList(self.annotate_psm.clone()),
            "Annotations added to each PSM.",
            tags,
        )?;
        let annotations: Vec<String> = ANNOTATIONS.iter().map(|&s| s.to_owned()).collect();
        p.set_valid_strings("annotate:PSM", &annotations)?;
        p.set_section_description("annotate", "Annotation Options")?;
        p.set_value(
            "peptide:min_size",
            int_value(self.peptide_min_size, "min_size")?,
            "Minimum size a peptide must have after digestion to be considered in the search.",
            tags,
        )?;
        p.set_value(
            "peptide:max_size",
            int_value(self.peptide_max_size, "max_size")?,
            "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).",
            tags,
        )?;
        p.set_value(
            "peptide:missed_cleavages",
            int_value(self.peptide_missed_cleavages, "missed_cleavages")?,
            "Number of missed cleavages.",
            tags,
        )?;
        p.set_value(
            "peptide:enzyme_specificity",
            ParamValue::String(self.peptide_enzyme_specificity.name().into()),
            SPECIFICITY_DESCRIPTION,
            tags,
        )?;
        p.set_valid_strings(
            "peptide:enzyme_specificity",
            &["full".into(), "semi".into(), "none".into()],
        )?;
        p.set_value(
            "peptide:motif",
            ParamValue::String(self.peptide_motif.clone()),
            "If set, only peptides that contain this motif (provided as RegEx) will be considered.",
            tags,
        )?;
        p.set_section_description("peptide", "Peptide Options")?;
        p.set_value(
            "report:top_hits",
            int_value(self.report_top_hits, "top_hits")?,
            "Maximum number of top scoring hits per spectrum that are reported.",
            tags,
        )?;
        p.set_section_description("report", "Reporting Options")?;
        p.set_value(
            "FDR:PSM",
            ParamValue::Float(self.fdr_psm),
            "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, set to 0 to disable filtering and report all PSMs with q-values). Requires '-decoys' to be set.",
            tags,
        )?;
        p.set_min_float("FDR:PSM", 0.0)?;
        p.set_max_float("FDR:PSM", 1.0)?;
        p.set_section_description("FDR", "False Discovery Rate control (requires decoys)")?;
        Ok(p)
    }

    /// Build an engine from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults; unknown ones are ignored, as the
    /// source only warns (see
    /// [`SimpleSearchEngineAlgorithm::from_param_with_warnings`]). Duplicate
    /// fixed or variable modifications are made unique, which, as in the
    /// source's `std::set` round trip, also sorts that list; a list without
    /// duplicates keeps its order. The native [`threads`](Self::threads)
    /// keeps its default.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source `Exception::InvalidParameter`), or
    /// when a count parameter is negative. The source stores those in `Size`
    /// members, where a negative value wraps to a huge count; this port refuses
    /// it.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`SimpleSearchEngineAlgorithm::from_param`], also returning the
    /// warnings the source logs: unknown parameters and the duplicate
    /// modification notices of `updateMembers_`.
    ///
    /// # Errors
    ///
    /// As [`SimpleSearchEngineAlgorithm::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(SIMPLE_SEARCH_ENGINE_ALGORITHM_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let ((engine, mut notices), mut warnings) =
            handler.set_parameters_with(parameters, Self::from_complete_param)?;
        warnings.append(&mut notices);
        Ok((engine, warnings))
    }

    /// Source `updateMembers_` over a complete, validated parameter set.
    fn from_complete_param(p: &Param) -> Result<(Self, Vec<String>)> {
        let mut notices = Vec::new();
        let fixed = unique_modifications(
            p.value("modifications:fixed")?.as_string_list()?,
            "Duplicate fixed modification provided. Making them unique.",
            &mut notices,
        );
        let variable = unique_modifications(
            p.value("modifications:variable")?.as_string_list()?,
            "Duplicate variable modification provided. Making them unique.",
            &mut notices,
        );
        let engine = Self {
            precursor_mass_tolerance: p.value("precursor:mass_tolerance")?.to_f64()?,
            precursor_mass_tolerance_unit: MassToleranceUnit::from_name(
                p.value("precursor:mass_tolerance_unit")?.as_str()?,
            )?,
            precursor_min_charge: charge(p, "precursor:min_charge")?,
            precursor_max_charge: charge(p, "precursor:max_charge")?,
            precursor_isotopes: p.value("precursor:isotopes")?.as_integer_list()?.to_vec(),
            fragment_mass_tolerance: p.value("fragment:mass_tolerance")?.to_f64()?,
            fragment_mass_tolerance_unit: MassToleranceUnit::from_name(
                p.value("fragment:mass_tolerance_unit")?.as_str()?,
            )?,
            modifications_fixed: fixed,
            modifications_variable: variable,
            modifications_max_variable_mods_per_peptide: count(
                p,
                "modifications:variable_max_per_peptide",
            )?,
            enzyme: p.value("enzyme")?.as_str()?.to_owned(),
            decoys: p.value("decoys")?.as_str()? == "true",
            fdr_psm: p.value("FDR:PSM")?.to_f64()?,
            annotate_psm: p.value("annotate:PSM")?.as_string_list()?.to_vec(),
            peptide_min_size: count(p, "peptide:min_size")?,
            peptide_max_size: count(p, "peptide:max_size")?,
            peptide_missed_cleavages: count(p, "peptide:missed_cleavages")?,
            peptide_enzyme_specificity: p.value("peptide:enzyme_specificity")?.as_str()?.parse()?,
            peptide_motif: p.value("peptide:motif")?.as_str()?.to_owned(),
            report_top_hits: count(p, "report:top_hits")?,
            threads: Threads::default(),
        };
        Ok((engine, notices))
    }

    /// Search the MS2 spectra in `in_spectra` against the protein database in
    /// `in_db` (source `search`).
    ///
    /// Spectra are loaded with MS level 2 only (mzML, Bruker TDF or Thermo
    /// RAW, as far as [`FileHandler`] reads them) and the FASTA database with
    /// [`FASTAFile`]; then [`SimpleSearchEngineAlgorithm::search_experiment`]
    /// runs. The two outputs are **not** treated symmetrically and the call is
    /// not additive: `protein_ids` is replaced by one fresh run for this
    /// search, whereas new PSMs are **appended** to `peptide_ids`, and the new
    /// run's identifier is stamped onto **every** element of `peptide_ids`,
    /// including pre-existing ones. Pass empty vectors.
    ///
    /// # Errors
    ///
    /// Returns load errors, and every error of
    /// [`SimpleSearchEngineAlgorithm::search_experiment`]. On error both
    /// outputs are unchanged.
    pub fn search(
        &self,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
    ) -> Result<ExitCodes> {
        self.search_reporting(
            in_spectra,
            in_db,
            protein_ids,
            peptide_ids,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`SimpleSearchEngineAlgorithm::search`], reporting progress to
    /// `progress` as the source's `ProgressLogger` base does.
    ///
    /// The sections are the source's: `Filtering spectra...` (0 to 1),
    /// `Generate decoys...` (0 to 1, with `decoys`), `Scoring peptide models
    /// against spectra...` (0 to the number of proteins, target and decoy,
    /// set to 1, 2, ... n once per protein, which is what a single-threaded
    /// source run reports) and `Post-processing PSMs...` (0 to 1). The
    /// per-protein values are reported after the parallel digestion rather
    /// than interleaved with it. Loading, peptide indexing and FDR report
    /// nothing, as in the source, whose helpers there own silent loggers.
    ///
    /// # Errors
    ///
    /// As [`SimpleSearchEngineAlgorithm::search`], and the errors of
    /// `progress`. A section is ended even when its body fails.
    pub fn search_with_progress(
        &self,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
        progress: &mut ProgressLogger,
    ) -> Result<ExitCodes> {
        self.search_reporting(
            in_spectra,
            in_db,
            protein_ids,
            peptide_ids,
            &mut ProgressReporter::new(Some(progress)),
        )
    }

    fn search_reporting(
        &self,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<ExitCodes> {
        let mut options = PeakFileOptions::default();
        options.clear_ms_levels();
        options.add_ms_level(2)?;
        let spectra = FileHandler::load_experiment_with_options(
            in_spectra,
            &[FileType::MzMl, FileType::BrukerTdf, FileType::Raw],
            &options,
        )?;
        let fasta_db = FASTAFile::new().load(in_db)?;
        self.search_experiment_reporting(
            spectra,
            fasta_db,
            in_spectra,
            in_db,
            protein_ids,
            peptide_ids,
            reporter,
        )
    }

    /// The search of [`SimpleSearchEngineAlgorithm::search`] on spectra and a
    /// database already in memory (native; the source only has the
    /// file-based entry point).
    ///
    /// `in_spectra` and `in_db` are only recorded: the former as the primary
    /// MS run path (unless the experiment names a single existing mzML file
    /// itself), the latter as the search database. `spectra` should hold MS2
    /// spectra only, as the loader of [`SimpleSearchEngineAlgorithm::search`]
    /// provides.
    ///
    /// Candidate peptides containing `X`, `B` or `Z` are skipped. A peptide
    /// seen in several proteins is scored once. Hits are collected per spectrum
    /// in the order a single-threaded source run produces them; the top-N
    /// selection uses a stable sort, so candidates that tie on every key of
    /// [`AnnotatedHit::has_better_score`] (possible only for one peptide
    /// matching one spectrum under two isotope offsets, which needs a Dalton
    /// tolerance wider than the isotope spacing) keep that order, where the
    /// source's `partial_sort` leaves it unspecified.
    ///
    /// With `decoys`, each protein gets a `DECOY_`-prefixed reversed copy
    /// (plain reverse for `none` specificity, reversal within enzymatic
    /// peptides otherwise), and targets and decoys are shuffled with the
    /// source's default-seeded `RandomShuffler`.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown enzyme or modification, an invalid
    /// motif, a spectrum whose precursor m/z or charge makes a mass or error
    /// non-finite, and the errors and resource ceilings of the components
    /// (digestion, modification generation, theoretical spectra, HyperScore,
    /// deisotoping, peptide indexing, FDR). The source maps peptide-indexer
    /// failures to [`ExitCodes::UnexpectedResult`] or
    /// [`ExitCodes::UnknownError`]; they are errors here. On error both
    /// outputs are unchanged.
    pub fn search_experiment(
        &self,
        spectra: MSExperiment,
        fasta_db: Vec<FASTAEntry>,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
    ) -> Result<ExitCodes> {
        self.search_experiment_reporting(
            spectra,
            fasta_db,
            in_spectra,
            in_db,
            protein_ids,
            peptide_ids,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`SimpleSearchEngineAlgorithm::search_experiment`] with the progress
    /// sections of [`SimpleSearchEngineAlgorithm::search_with_progress`].
    ///
    /// # Errors
    ///
    /// As [`SimpleSearchEngineAlgorithm::search_experiment`], and the errors of
    /// `progress`.
    #[allow(clippy::too_many_arguments)]
    pub fn search_experiment_with_progress(
        &self,
        spectra: MSExperiment,
        fasta_db: Vec<FASTAEntry>,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
        progress: &mut ProgressLogger,
    ) -> Result<ExitCodes> {
        self.search_experiment_reporting(
            spectra,
            fasta_db,
            in_spectra,
            in_db,
            protein_ids,
            peptide_ids,
            &mut ProgressReporter::new(Some(progress)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn search_experiment_reporting(
        &self,
        mut spectra: MSExperiment,
        fasta_db: Vec<FASTAEntry>,
        in_spectra: &str,
        in_db: &str,
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<ExitCodes> {
        let motif = BoostRegex::new(&self.peptide_motif)?;
        let fixed_modifications = resolve_modifications(&self.modifications_fixed)?;
        let variable_modifications = resolve_modifications(&self.modifications_variable)?;
        let protease = Protease::from_name(&self.enzyme)?;

        spectra.sort_spectra(true)?;
        reporter.section(0, 1, "Filtering spectra...", |_| {
            Self::preprocess_spectra_with_threads(
                &mut spectra,
                self.fragment_mass_tolerance,
                self.fragment_mass_tolerance_unit == MassToleranceUnit::Ppm,
                self.threads,
            )
        })?;
        let precursors = map_precursor_masses_to_scans(
            self.precursor_min_charge,
            self.precursor_max_charge,
            &self.precursor_isotopes,
            self.peptide_min_size,
            &spectra,
        )?;

        let mut fasta_db = fasta_db;
        if self.decoys {
            fasta_db = reporter.section(0, 1, "Generate decoys...", |_| {
                self.append_decoys(fasta_db, protease)
            })?;
        }
        let digestor = ProteaseDigestion {
            enzyme: protease,
            missed_cleavages: self.peptide_missed_cleavages,
            min_length: self.peptide_min_size,
            max_length: (self.peptide_max_size > 0).then_some(self.peptide_max_size),
            specificity: self.peptide_enzyme_specificity,
            ..ProteaseDigestion::default()
        };
        let proteins = progress_value(fasta_db.len())?;
        let (peptide_count, mut annotated_hits) = reporter.section(
            0,
            proteins,
            "Scoring peptide models against spectra...",
            |reporter| {
                self.score_database(
                    &fasta_db,
                    &digestor,
                    &motif,
                    &fixed_modifications,
                    &variable_modifications,
                    &precursors,
                    &spectra,
                    reporter,
                )
            },
        )?;
        let _ = crate::openms_log_info!("Proteins: {}", fasta_db.len());
        let _ = crate::openms_log_info!("Peptides: {}", peptide_count);
        let _ = crate::openms_log_info!("Processed peptides: {}", peptide_count);

        let mut new_proteins = Vec::new();
        let mut new_peptides = peptide_ids.clone();
        reporter.section(0, 1, "Post-processing PSMs...", |_| {
            self.post_process_hits(
                &spectra,
                &mut annotated_hits,
                &mut new_proteins,
                &mut new_peptides,
                in_db,
            )
        })?;
        let (spectra_data, spectra_data_raw) = primary_ms_run_paths(in_spectra, &spectra);
        if let Some(raw) = spectra_data_raw {
            new_proteins[0].raw_ms_run_paths = vec![raw];
        }
        new_proteins[0].primary_ms_run_paths = spectra_data;

        if fasta_db.is_empty() {
            // Source: PeptideIndexing reports DATABASE_EMPTY after the outputs
            // were written.
            *protein_ids = new_proteins;
            *peptide_ids = new_peptides;
            return Ok(ExitCodes::InputFileEmpty);
        }
        let indexer = PeptideIndexing {
            decoy_rule: DecoyRule::Prefix("DECOY_".into()),
            enzyme: Some(protease),
            specificity: Some(self.peptide_enzyme_specificity),
            missing_decoy_action: MissingDecoyAction::Silent,
            ..PeptideIndexing::default()
        };
        indexer.run(&fasta_db, &mut new_proteins, &mut new_peptides)?;

        if self.fdr_psm > 0.0 && self.decoys {
            FalseDiscoveryRate::default().apply_peptides(&mut new_peptides, false)?;
            id_filter::filter_peptides_by_score(&mut new_peptides, self.fdr_psm)?;
            id_filter::remove_decoy_peptide_hits(&mut new_peptides)?;
            id_filter::remove_empty_peptide_identifications(&mut new_peptides)?;
            id_filter::remove_unreferenced_proteins(&mut new_proteins, &new_peptides)?;
        } else if self.fdr_psm > 0.0 {
            let _ = crate::openms_log_warn!(
                "FDR:PSM is set but decoys are disabled. Skipping FDR filtering."
            );
        }
        *protein_ids = new_proteins;
        *peptide_ids = new_peptides;
        Ok(ExitCodes::ExecutionOk)
    }

    /// Digest every protein, claim each distinct peptide in protein order (as
    /// the source's `processed_petides` set does) and score its modification
    /// variants, keeping per-spectrum hit lists pruned to the top N. Reports
    /// the per-protein progress values of the scoring section. Returns the
    /// number of distinct peptides and the hit lists.
    #[allow(clippy::too_many_arguments)]
    fn score_database(
        &self,
        fasta_db: &[FASTAEntry],
        digestor: &ProteaseDigestion,
        motif: &BoostRegex,
        fixed_modifications: &[Arc<ResidueModification>],
        variable_modifications: &[Arc<ResidueModification>],
        precursors: &[PrecursorEntry],
        spectra: &MSExperiment,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<(usize, Vec<Vec<AnnotatedHit>>)> {
        let digests = map_collect(fasta_db, self.threads, |entry| {
            let mut kept = Vec::new();
            for peptide in digestor.digest_unmodified(&entry.sequence)? {
                if peptide.contains(['X', 'B', 'Z']) {
                    continue;
                }
                if !self.peptide_motif.is_empty() && !motif.is_full_match(peptide.as_bytes())? {
                    continue;
                }
                kept.push(peptide.to_owned());
            }
            Ok::<_, Error>(kept)
        });
        let mut seen = BTreeSet::new();
        let mut peptides: Vec<Arc<str>> = Vec::new();
        for (index, digest) in digests.into_iter().enumerate() {
            reporter.set_count(index + 1)?;
            for peptide in digest? {
                if seen.insert(peptide.clone()) {
                    peptides.push(Arc::from(peptide));
                }
            }
        }
        let scorer = Scorer {
            engine: self,
            fixed: fixed_modifications,
            variable: variable_modifications,
            precursors,
            spectra,
            generator: TheoreticalSpectrumGenerator {
                add_first_prefix_ion: true,
                add_metainfo: true,
                ..TheoreticalSpectrumGenerator::default()
            },
            hyperscore: HyperScore {
                tolerance: self
                    .fragment_mass_tolerance_unit
                    .tolerance(self.fragment_mass_tolerance),
            },
        };
        let mut annotated_hits: Vec<Vec<AnnotatedHit>> = vec![Vec::new(); spectra.spectra.len()];
        for batch in peptides.chunks(PEPTIDE_BATCH) {
            for hits in map_collect(batch, self.threads, |peptide| scorer.score(peptide)) {
                for (scan_index, hit) in hits? {
                    let slot = &mut annotated_hits[scan_index];
                    slot.push(hit);
                    // Prevent the vector from growing indefinitely, but do not
                    // shrink it on every push.
                    if slot.len() >= 2 * self.report_top_hits {
                        keep_best(slot, self.report_top_hits);
                    }
                }
            }
        }
        Ok((peptides.len(), annotated_hits))
    }

    /// Append `DECOY_`-prefixed reversed proteins and shuffle targets and
    /// decoys with the source's default-seeded `RandomShuffler`, so equal
    /// target and decoy scores carry no global bias.
    fn append_decoys(
        &self,
        mut fasta_db: Vec<FASTAEntry>,
        protease: Protease,
    ) -> Result<Vec<FASTAEntry>> {
        let generator = DecoyGenerator::with_seed(0);
        let targets = fasta_db.len();
        fasta_db.reserve(targets);
        for index in 0..targets {
            let mut entry = fasta_db[index].clone();
            let protein = AASequence::parse(&entry.sequence)?;
            let decoy = if self.peptide_enzyme_specificity == DigestionSpecificity::None {
                generator.reverse_protein(&protein)?
            } else {
                generator.reverse_peptides(&protein, protease)?
            };
            entry.sequence = decoy.as_str().to_owned();
            entry.identifier = format!("DECOY_{}", entry.identifier);
            fasta_db.push(entry);
        }
        // boost::mt19937_64's default seed.
        let mut shuffler = RandomShuffler::with_seed(5489);
        let mut draws = fasta_db.len().saturating_mul(64).saturating_add(1024);
        shuffler.shuffle_limited(&mut fasta_db, &mut draws)?;
        Ok(fasta_db)
    }

    /// Preprocess MS2 spectra in place: filter, deisotope, decharge (source
    /// `preprocessSpectra_`, protected there).
    ///
    /// Removes peaks below the [`ThresholdMower`] default threshold (0.05; the
    /// source comments "remove 0 intensities" but uses that default), scales
    /// each spectrum's maximum to one, sorts spectra by retention time, then
    /// per spectrum sorts by m/z, deisotopes and converts to charge 1
    /// (charges 1-3, 3-10 isotope peaks, keeping unassigned peaks, with the
    /// source's shared isotope extensions), keeps the 20 most intense peaks per
    /// jumping 100 Th window and the 400 most intense overall, and sorts by
    /// m/z again. Deisotoping is skipped, with a logged warning, when the
    /// fragment tolerance exceeds its 100 ppm / 0.1 Da limit (low-resolution
    /// data, OpenMS#9619).
    ///
    /// # Arguments
    ///
    /// * `fragment_mass_tolerance` - tolerance for deisotoping and decharging.
    /// * `fragment_mass_tolerance_unit_ppm` - `true` for ppm, `false` for Th.
    ///
    /// # Errors
    ///
    /// Propagates the filter errors (invalid spectra, resource ceilings). On
    /// error `exp` is unchanged.
    pub fn preprocess_spectra(
        exp: &mut MSExperiment,
        fragment_mass_tolerance: f64,
        fragment_mass_tolerance_unit_ppm: bool,
    ) -> Result<()> {
        Self::preprocess_spectra_with_threads(
            exp,
            fragment_mass_tolerance,
            fragment_mass_tolerance_unit_ppm,
            Threads::default(),
        )
    }

    fn preprocess_spectra_with_threads(
        exp: &mut MSExperiment,
        fragment_mass_tolerance: f64,
        fragment_mass_tolerance_unit_ppm: bool,
        threads: Threads,
    ) -> Result<()> {
        let mut work = exp.clone();
        ThresholdMower::default().filter_experiment(&mut work)?;
        Normalizer::default().filter_experiment(&mut work)?;
        work.sort_spectra(false)?;
        let tolerance = if fragment_mass_tolerance_unit_ppm {
            Tolerance::Ppm(fragment_mass_tolerance)
        } else {
            Tolerance::Absolute(fragment_mass_tolerance)
        };
        let do_deisotope = Deisotoper::is_tolerance_supported(tolerance);
        if !do_deisotope {
            let _ = crate::openms_log_warn!(
                "[SimpleSearchEngine] Fragment tolerance {} {} exceeds the deisotoping limit (100 ppm / 0.1 Da); skipping MS2 deisotoping (expected for low-resolution data).",
                fragment_mass_tolerance,
                if fragment_mass_tolerance_unit_ppm {
                    "ppm"
                } else {
                    "Da"
                }
            );
        }
        let deisotoper = Deisotoper {
            tolerance,
            min_charge: 1,
            max_charge: 3,
            keep_only_deisotoped: false,
            min_isotope_peaks: 3,
            max_isotope_peaks: 10,
            make_single_charged: true,
            annotate_charge: false,
            annotate_isotope_peak_count: false,
            annotate_features: false,
            use_decreasing_model: true,
            start_intensity_check: 2,
            add_up_intensity: false,
            allow_shared_isotopes: true,
            ..Deisotoper::default()
        };
        let window_mower = WindowMower {
            window_size: 100.0,
            peak_count: 20,
            method: WindowMowerMethod::Jumping,
            ..WindowMower::default()
        };
        let nlargest = NLargest { n: 400 };
        let processed = map_collect(&work.spectra, threads, |spectrum| {
            let mut spectrum = spectrum.clone();
            spectrum.sort_by_position()?;
            if do_deisotope {
                spectrum = deisotoper.deisotope(&spectrum)?.spectrum;
            }
            window_mower.filter_spectrum(&mut spectrum)?;
            nlargest.filter_spectrum(&mut spectrum)?;
            spectrum.sort_by_position()?;
            Ok::<_, Error>(spectrum)
        });
        work.spectra = processed.into_iter().collect::<Result<Vec<_>>>()?;
        *exp = work;
        Ok(())
    }

    /// Materialise the top-N scored candidates per spectrum as
    /// [`PeptideIdentification`]s and one [`ProteinIdentification`] (source
    /// `postProcessHits_`, protected there).
    ///
    /// Keeps the best [`report_top_hits`](Self::report_top_hits) of each list
    /// in `annotated_hits` (one list per spectrum of `exp`), re-applies each
    /// hit's modification variant, adds the `annotate:PSM` meta values, sorts
    /// each identification's hits and appends one identification per spectrum
    /// with hits to `peptide_ids`, in spectrum order. `protein_ids` is replaced
    /// by one run named `SSE_<date time>`, whose identifier is stamped on every
    /// element of `peptide_ids`, with the search settings (tolerances,
    /// modifications, enzyme, `database_name`) and the `extra_features` list.
    /// The source passes its own members as further arguments; this reads
    /// them from `self`.
    ///
    /// Fragment annotations, the fragment error, the longest ion run and the
    /// matched ion current come from a [`SpectrumAlignment`] with its default
    /// 0.3 Th tolerance, whatever the fragment tolerance, against a b/y
    /// spectrum at charges 1 to `min(charge - 1, 2)` (1 for charge 1). The
    /// precursor error corrects the observed m/z by
    /// `isotope_error * PROTON_MASS_U / charge`, as the source does, although
    /// the precursor mass was corrected by the C13 spacing.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `annotated_hits` does not have one
    /// list per spectrum, a spectrum with hits has no precursor, a stored
    /// modification index no longer exists, or a computed value is not finite
    /// (for example a precursor charge of zero), and propagates spectrum
    /// generation and alignment errors. On error the outputs are unchanged.
    pub fn post_process_hits(
        &self,
        exp: &MSExperiment,
        annotated_hits: &mut [Vec<AnnotatedHit>],
        protein_ids: &mut Vec<ProteinIdentification>,
        peptide_ids: &mut Vec<PeptideIdentification>,
        database_name: &str,
    ) -> Result<()> {
        if annotated_hits.len() != exp.spectra.len() {
            return Err(bad("post-processing needs one hit list per spectrum"));
        }
        for hits in annotated_hits.iter_mut() {
            keep_best(hits, self.report_top_hits);
        }
        let flags = AnnotationFlags::new(&self.annotate_psm);
        let fixed = resolve_modifications(&self.modifications_fixed)?;
        let variable = resolve_modifications(&self.modifications_variable)?;
        let indices: Vec<usize> = (0..annotated_hits.len())
            .filter(|&i| !annotated_hits[i].is_empty())
            .collect();
        let hits: &[Vec<AnnotatedHit>] = annotated_hits;
        let identifications = map_collect(&indices, self.threads, |&scan_index| {
            self.identification(
                exp,
                scan_index,
                &hits[scan_index],
                &flags,
                &fixed,
                &variable,
            )
        })
        .into_iter()
        .collect::<Result<Vec<_>>>()?;

        let now = DateTime::now();
        let identifier = format!("SSE_{}", now.get());
        let mut protein = ProteinIdentification {
            date_time: Some(now.iso_string()),
            search_engine: "SimpleSearchEngine".into(),
            search_engine_version: crate::CORE_SDK_VERSION.into(),
            identifier: identifier.clone(),
            ..ProteinIdentification::default()
        };
        let mut output = peptide_ids.clone();
        output.extend(identifications);
        for id in &mut output {
            id.identifier.clone_from(&identifier);
        }
        let mut feature_set = vec!["score".to_owned(), "isotope_error".to_owned()];
        if flags.fragment_error_ppm {
            feature_set.push(user_param::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM.into());
        }
        if flags.prefix_fraction {
            feature_set.push(user_param::MATCHED_PREFIX_IONS_FRACTION.into());
        }
        if flags.suffix_fraction {
            feature_set.push(user_param::MATCHED_SUFFIX_IONS_FRACTION.into());
        }
        let mut parameters = SearchParameters {
            database: database_name.into(),
            charges: format!(
                "{}:{}",
                self.precursor_min_charge, self.precursor_max_charge
            ),
            mass_type: PeakMassType::Monoisotopic,
            fixed_modifications: self.modifications_fixed.clone(),
            variable_modifications: self.modifications_variable.clone(),
            missed_cleavages: u32::try_from(self.peptide_missed_cleavages)
                .map_err(|_| bad("missed cleavages exceed the source's 32-bit field"))?,
            fragment_tolerance: self
                .fragment_mass_tolerance_unit
                .tolerance(self.fragment_mass_tolerance),
            precursor_tolerance: self
                .precursor_mass_tolerance_unit
                .tolerance(self.precursor_mass_tolerance),
            digestion_enzyme: ProteaseDB::global().get_enzyme(&self.enzyme)?.name().into(),
            enzyme_specificity: match self.peptide_enzyme_specificity {
                DigestionSpecificity::Full => EnzymeTermSpecificity::Full,
                DigestionSpecificity::Semi => EnzymeTermSpecificity::Semi,
                DigestionSpecificity::None => EnzymeTermSpecificity::None,
                // `SPEC_UNKNOWN`, `SPEC_NOCTERM` and `SPEC_NONTERM` are values
                // `ProteaseDigestion::digest` refuses; the search never holds one.
                other => {
                    return Err(bad(format!(
                        "unsupported peptide enzyme specificity {}",
                        other.name()
                    )));
                }
            },
            ..SearchParameters::default()
        };
        parameters
            .metadata
            .insert("extra_features".into(), feature_set.join(",").into());
        protein.search_parameters = parameters;
        // The source reads `exp[0]` unguarded; with pre-existing identifications
        // and no spectra that is undefined, and here it is skipped.
        if let Some(first) = exp.spectra.first() {
            if !output.is_empty()
                && output
                    .iter()
                    .all(|id| id.metadata.contains_key(user_param::IM))
            {
                protein.metadata.insert(
                    user_param::IM.into(),
                    first.drift_time_unit_as_string().into(),
                );
            }
        }
        *protein_ids = vec![protein];
        *peptide_ids = output;
        Ok(())
    }

    /// One spectrum's identification from its (already top-N) hits.
    fn identification(
        &self,
        exp: &MSExperiment,
        scan_index: usize,
        hits: &[AnnotatedHit],
        flags: &AnnotationFlags,
        fixed: &[Arc<ResidueModification>],
        variable: &[Arc<ResidueModification>],
    ) -> Result<PeptideIdentification> {
        let spec = &exp.spectra[scan_index];
        let precursor = spec
            .precursors
            .first()
            .ok_or_else(|| bad("a spectrum with hits has no precursor"))?;
        let mz = precursor.mz;
        let charge = precursor.charge;
        let mut pi = PeptideIdentification {
            score_type: "ln(hyperscore)".into(),
            higher_score_better: true,
            rt: Some(spec.rt),
            mz: Some(mz),
            ..PeptideIdentification::default()
        };
        pi.set_spectrum_reference(spec.native_id.clone());
        pi.metadata.insert(
            "scan_index".into(),
            MetaValue::from(
                u32::try_from(scan_index).map_err(|_| bad("scan index exceeds 32 bits"))?,
            ),
        );
        if ImTypes::determine_im_format(spec) == IonMobilityFormat::PerSpectrum {
            pi.metadata
                .insert(user_param::IM.into(), MetaValue::try_from(spec.drift_time)?);
        }
        let generator = ModifiedPeptideGenerator::default();
        for ah in hits {
            let mut ph = PeptideHit {
                charge,
                score: ah.score,
                ..PeptideHit::default()
            };
            ph.metadata
                .insert(user_param::ISOTOPE_ERROR.into(), ah.isotope_error.into());
            let mut aas = AASequence::parse(&ah.sequence)?;
            generator.apply_fixed_modifications(fixed, &mut aas)?;
            let mut variants = Vec::new();
            generator.apply_variable_modifications(
                variable,
                &aas,
                self.modifications_max_variable_mods_per_peptide,
                &mut variants,
                true,
            )?;
            let peptide = variants
                .get(ah.peptide_mod_index)
                .ok_or_else(|| bad("stored modification variant no longer exists"))?
                .clone();

            let mut alignment = Vec::new();
            let mut theoretical = MSSpectrum::default();
            if flags.need_alignment() {
                let generator = TheoreticalSpectrumGenerator {
                    add_metainfo: true,
                    add_first_prefix_ion: true,
                    ..TheoreticalSpectrumGenerator::default()
                };
                let max_frag_z = if charge >= 2 { (charge - 1).min(2) } else { 1 };
                theoretical = generator.generate(
                    &peptide,
                    1,
                    u8::try_from(max_frag_z).map_err(|_| bad("fragment charge out of range"))?,
                    None,
                )?;
                alignment = SpectrumAlignment::default().align(&theoretical, spec)?;
            }
            if flags.fragment_error_ppm {
                let mut errors = alignment
                    .iter()
                    .map(|&(t, e)| Ok(ppm(spec.peaks[e].mz, theoretical.peaks[t].mz)?.abs()))
                    .collect::<Result<Vec<f64>>>()?;
                let median_ppm_error = if errors.is_empty() {
                    0.0
                } else {
                    median(&mut errors)?
                };
                ph.metadata.insert(
                    user_param::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM.into(),
                    MetaValue::try_from(median_ppm_error)?,
                );
            }
            if flags.precursor_error_ppm {
                let theo_mz = peptide.mz(charge)?;
                let observed = mz - f64::from(ah.isotope_error) * PROTON_MASS_U / f64::from(charge);
                ph.metadata.insert(
                    user_param::PRECURSOR_ERROR_PPM_USERPARAM.into(),
                    MetaValue::try_from(ppm(observed, theo_mz)?)?,
                );
            }
            if flags.prefix_fraction {
                ph.metadata.insert(
                    user_param::MATCHED_PREFIX_IONS_FRACTION.into(),
                    MetaValue::try_from(ah.prefix_fraction)?,
                );
            }
            if flags.suffix_fraction {
                ph.metadata.insert(
                    user_param::MATCHED_SUFFIX_IONS_FRACTION.into(),
                    MetaValue::try_from(ah.suffix_fraction)?,
                );
            }
            if flags.num_matched_peaks {
                ph.metadata.insert(
                    user_param::NUM_MATCHED_PEAKS.into(),
                    (i32::from(ah.matched_prefix_ions) + i32::from(ah.matched_suffix_ions)).into(),
                );
            }
            if flags.matched_prefix_ions {
                ph.metadata.insert(
                    user_param::MATCHED_PREFIX_IONS.into(),
                    i32::from(ah.matched_prefix_ions).into(),
                );
            }
            if flags.matched_suffix_ions {
                ph.metadata.insert(
                    user_param::MATCHED_SUFFIX_IONS.into(),
                    i32::from(ah.matched_suffix_ions).into(),
                );
            }
            if flags.fragment_annotations || flags.longest_ion_run || flags.matched_ion_current {
                annotate_alignment(&mut ph, spec, &theoretical, &alignment, flags)?;
            }
            ph.sequence = peptide;
            pi.hits.push(ph);
        }
        pi.sort()?;
        Ok(pi)
    }
}

/// The alignment-based annotations: fragment annotations, matched ion
/// current and the longest consecutive ion ladder.
fn annotate_alignment(
    ph: &mut PeptideHit,
    spec: &MSSpectrum,
    theoretical: &MSSpectrum,
    alignment: &[(usize, usize)],
    flags: &AnnotationFlags,
) -> Result<()> {
    let names = theoretical
        .string_data_arrays
        .first()
        .ok_or_else(|| bad("theoretical spectrum lacks ion names"))?;
    let charges = theoretical
        .integer_data_arrays
        .first()
        .ok_or_else(|| bad("theoretical spectrum lacks ion charges"))?;
    let mut annotations = Vec::with_capacity(alignment.len());
    let (mut prefix, mut suffix) = (Vec::new(), Vec::new());
    let mut matched_ion_current = 0.0_f64;
    for &(t, e) in alignment {
        let name = names
            .data
            .get(t)
            .ok_or_else(|| bad("ion name index out of range"))?;
        if flags.fragment_annotations {
            annotations.push(PeakAnnotation {
                mz: spec.peaks[e].mz,
                intensity: f64::from(spec.peaks[e].intensity),
                annotation: name.clone(),
                charge: *charges
                    .data
                    .get(t)
                    .ok_or_else(|| bad("ion charge index out of range"))?,
            });
        }
        if flags.matched_ion_current {
            matched_ion_current += f64::from(spec.peaks[e].intensity);
        }
        if flags.longest_ion_run && name.len() >= 2 {
            if let Some((is_prefix, ordinal)) = ion_ordinal(name) {
                if is_prefix {
                    prefix.push(ordinal);
                } else {
                    suffix.push(ordinal);
                }
            }
        }
    }
    if flags.fragment_annotations {
        ph.peak_annotations = annotations;
    }
    if flags.matched_ion_current {
        ph.metadata.insert(
            user_param::MATCHED_ION_CURRENT.into(),
            MetaValue::try_from(matched_ion_current)?,
        );
    }
    if flags.longest_ion_run {
        let longest = longest_run(&mut prefix).max(longest_run(&mut suffix));
        ph.metadata.insert(
            user_param::LONGEST_PEPTIDE_ION_SEQUENCE.into(),
            longest.into(),
        );
    }
    Ok(())
}

/// Series and ordinal of an ion name: `"b5"`, `"y3-H2O1+"`, `"c12++"` give
/// (prefix, 5), (suffix, 3), (prefix, 12). Names without a digit after the
/// series letter, or of another series, give `None`. An ordinal beyond `i32`
/// saturates, where the source's `toInt32` throws.
fn ion_ordinal(name: &str) -> Option<(bool, i32)> {
    let bytes = name.as_bytes();
    let is_prefix = match bytes.first()? {
        b'a' | b'b' | b'c' => true,
        b'x' | b'y' | b'z' => false,
        _ => return None,
    };
    let digits = bytes[1..].iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    let ordinal = name[1..=digits].bytes().fold(0_i64, |n, d| {
        (n * 10 + i64::from(d - b'0')).min(i64::from(i32::MAX))
    });
    Some((is_prefix, ordinal as i32))
}

/// Longest run of consecutive distinct ordinals, 0 when empty.
fn longest_run(values: &mut Vec<i32>) -> i32 {
    if values.is_empty() {
        return 0;
    }
    values.sort_unstable();
    values.dedup();
    let (mut best, mut run) = (1, 1);
    for pair in values.windows(2) {
        if i64::from(pair[1]) == i64::from(pair[0]) + 1 {
            run += 1;
            best = best.max(run);
        } else {
            run = 1;
        }
    }
    best
}

/// Which `annotate:PSM` annotations are requested; `ALL` requests every one.
#[derive(Clone, Copy, Debug, Default)]
struct AnnotationFlags {
    precursor_error_ppm: bool,
    fragment_error_ppm: bool,
    prefix_fraction: bool,
    suffix_fraction: bool,
    num_matched_peaks: bool,
    matched_prefix_ions: bool,
    matched_suffix_ions: bool,
    longest_ion_run: bool,
    matched_ion_current: bool,
    fragment_annotations: bool,
}

impl AnnotationFlags {
    fn new(requested: &[String]) -> Self {
        let has = |key: &str| requested.iter().any(|r| r == key);
        if has("ALL") {
            return Self {
                precursor_error_ppm: true,
                fragment_error_ppm: true,
                prefix_fraction: true,
                suffix_fraction: true,
                num_matched_peaks: true,
                matched_prefix_ions: true,
                matched_suffix_ions: true,
                longest_ion_run: true,
                matched_ion_current: true,
                fragment_annotations: true,
            };
        }
        Self {
            precursor_error_ppm: has(user_param::PRECURSOR_ERROR_PPM_USERPARAM),
            fragment_error_ppm: has(user_param::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM),
            prefix_fraction: has(user_param::MATCHED_PREFIX_IONS_FRACTION),
            suffix_fraction: has(user_param::MATCHED_SUFFIX_IONS_FRACTION),
            num_matched_peaks: has(user_param::NUM_MATCHED_PEAKS),
            matched_prefix_ions: has(user_param::MATCHED_PREFIX_IONS),
            matched_suffix_ions: has(user_param::MATCHED_SUFFIX_IONS),
            longest_ion_run: has(user_param::LONGEST_PEPTIDE_ION_SEQUENCE),
            matched_ion_current: has(user_param::MATCHED_ION_CURRENT),
            fragment_annotations: has(user_param::FRAGMENT_ANNOTATION_USERPARAM),
        }
    }

    /// Fragment error, fragment annotations, longest ion run and matched ion
    /// current all need the spectrum alignment.
    fn need_alignment(&self) -> bool {
        self.fragment_error_ppm
            || self.fragment_annotations
            || self.longest_ion_run
            || self.matched_ion_current
    }
}

/// One precursor mass hypothesis: `(mass, scan index, isotope offset)`.
type PrecursorEntry = (f64, usize, i32);

/// Map each MS2 spectrum's (isotope-corrected) precursor mass to its scan
/// index, sorted by mass with insertion order among equal masses, like the
/// source's `std::multimap` (`mapPrecursorMassesToScans`, file-local there).
///
/// Only spectra with exactly one precursor, a charge in
/// `[min_precursor_charge, max_precursor_charge]` and at least
/// `peptide_min_size` peaks are mapped: one entry per isotope offset `i`, with
/// mass `z * mz - z * PROTON_MASS_U - i * C13C12_MASSDIFF_U`.
fn map_precursor_masses_to_scans(
    min_precursor_charge: i32,
    max_precursor_charge: i32,
    precursor_isotopes: &[i32],
    peptide_min_size: usize,
    spectra: &MSExperiment,
) -> Result<Vec<PrecursorEntry>> {
    let mut entries = Vec::new();
    for (scan_index, spectrum) in spectra.spectra.iter().enumerate() {
        if spectrum.precursors.len() != 1 || spectrum.peaks.len() < peptide_min_size {
            continue;
        }
        let precursor = &spectrum.precursors[0];
        let z = precursor.charge;
        if z < min_precursor_charge || z > max_precursor_charge {
            continue;
        }
        for &i in precursor_isotopes {
            let mut mass = f64::from(z) * precursor.mz - f64::from(z) * PROTON_MASS_U;
            if i != 0 {
                mass -= f64::from(i) * C13C12_MASSDIFF_U;
            }
            if !mass.is_finite() {
                return Err(bad("precursor mass is not finite"));
            }
            entries.push((mass, scan_index, i));
        }
    }
    entries.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(entries)
}

/// Everything one peptide's scoring needs, shared across worker threads.
struct Scorer<'a> {
    engine: &'a SimpleSearchEngineAlgorithm,
    fixed: &'a [Arc<ResidueModification>],
    variable: &'a [Arc<ResidueModification>],
    precursors: &'a [PrecursorEntry],
    spectra: &'a MSExperiment,
    generator: TheoreticalSpectrumGenerator,
    hyperscore: HyperScore,
}

impl Scorer<'_> {
    /// Score every modification variant of one unmodified peptide against the
    /// spectra in its precursor window, in the source's loop order.
    fn score(&self, peptide: &Arc<str>) -> Result<Vec<(usize, AnnotatedHit)>> {
        let engine = self.engine;
        let generator = ModifiedPeptideGenerator::default();
        let mut aas = AASequence::parse(peptide)?;
        generator.apply_fixed_modifications(self.fixed, &mut aas)?;
        let mut variants = Vec::new();
        generator.apply_variable_modifications(
            self.variable,
            &aas,
            engine.modifications_max_variable_mods_per_peptide,
            &mut variants,
            true,
        )?;
        let mut hits = Vec::new();
        for (mod_pep_idx, candidate) in variants.iter().enumerate() {
            let mass = candidate.mono_mass()?;
            let (low, high) = match engine.precursor_mass_tolerance_unit {
                MassToleranceUnit::Ppm => {
                    let delta = mass * engine.precursor_mass_tolerance * 1e-6;
                    (mass - delta, mass + delta)
                }
                MassToleranceUnit::Da => (
                    mass - engine.precursor_mass_tolerance,
                    mass + engine.precursor_mass_tolerance,
                ),
            };
            let from = self.precursors.partition_point(|e| e.0 < low);
            let to = self.precursors.partition_point(|e| e.0 <= high).max(from);
            if from == to {
                continue;
            }
            let theo = self.generator.generate(candidate, 1, 1, None)?;
            for &(_, scan_index, isotope) in &self.precursors[from..to] {
                let detail = self
                    .hyperscore
                    .compute_with_detail(&self.spectra.spectra[scan_index], &theo)?;
                if detail.score == 0.0 {
                    continue;
                }
                let length = peptide.len() as f64;
                hits.push((
                    scan_index,
                    AnnotatedHit {
                        sequence: Arc::clone(peptide),
                        peptide_mod_index: mod_pep_idx,
                        score: detail.score,
                        prefix_fraction: (detail.matched_prefix_ions as f64 / length) as f32,
                        suffix_fraction: (detail.matched_suffix_ions as f64 / length) as f32,
                        mean_error: detail.mean_error as f32,
                        isotope_error: isotope,
                        matched_prefix_ions: saturate_u16(detail.matched_prefix_ions),
                        matched_suffix_ions: saturate_u16(detail.matched_suffix_ions),
                    },
                ));
            }
        }
        Ok(hits)
    }
}

/// Sort best first (stable) and keep `n`.
fn keep_best(hits: &mut Vec<AnnotatedHit>, n: usize) {
    hits.sort_by(AnnotatedHit::better_first);
    hits.truncate(n);
}

/// The source casts matched-ion counts to `uint16_t`, which wraps; peptides
/// are far below that, and this saturates instead.
fn saturate_u16(value: usize) -> u16 {
    u16::try_from(value).unwrap_or(u16::MAX)
}

/// The primary MS run annotation of source
/// `ProteinIdentification::setPrimaryMSRunPath({in_spectra}, spectra)`:
/// `(spectra_data, spectra_data_raw)`.
///
/// When the experiment names exactly one source file: an existing mzML file
/// becomes `spectra_data` alone; a RAW file becomes `spectra_data_raw`, and
/// `in_spectra` still becomes `spectra_data`. Otherwise `in_spectra` is
/// `spectra_data`. Paths are assembled as `MSExperiment::getPrimaryMSRunPath`
/// does: `path + separator + name`, with `\` only for a backslash-only path.
fn primary_ms_run_paths(in_spectra: &str, spectra: &MSExperiment) -> (Vec<String>, Option<String>) {
    let mut paths = Vec::new();
    for file in &spectra.settings.source_files {
        if file.path.is_empty() || file.name.is_empty() {
            continue;
        }
        let actual = file.path.strip_prefix("file:///").unwrap_or(&file.path);
        let separator = if actual.contains('\\') && !actual.contains('/') {
            "\\"
        } else {
            "/"
        };
        paths.push(format!("{}{separator}{}", file.path, file.name));
    }
    let mut raw = None;
    if paths.len() == 1 {
        let kind = type_by_file_name(&paths[0]);
        if kind == FileType::MzMl && std::path::Path::new(&paths[0]).exists() {
            return (paths, None);
        }
        if kind == FileType::Raw {
            raw = Some(paths[0].clone());
        }
    }
    (vec![in_spectra.to_owned()], raw)
}

fn resolve_modifications(names: &[String]) -> Result<Vec<Arc<ResidueModification>>> {
    let names: Vec<&str> = names.iter().map(String::as_str).collect();
    ModifiedPeptideGenerator::get_modifications(&names)
}

/// The list itself, or its sorted unique set with a notice when it holds
/// duplicates (source `updateMembers_`).
fn unique_modifications(list: &[String], notice: &str, notices: &mut Vec<String>) -> Vec<String> {
    let unique: BTreeSet<&String> = list.iter().collect();
    if unique.len() == list.len() {
        return list.to_vec();
    }
    let _ = crate::openms_log_warn!("{notice}");
    notices.push(notice.to_owned());
    unique.into_iter().cloned().collect()
}

fn charge(p: &Param, key: &str) -> Result<i32> {
    i32::try_from(p.value(key)?.to_i64()?)
        .map_err(|_| bad(format!("{key} does not fit a 32-bit charge")))
}

fn count(p: &Param, key: &str) -> Result<usize> {
    usize::try_from(p.value(key)?.to_i64()?).map_err(|_| bad(format!("{key} must not be negative")))
}

fn int_value(value: usize, name: &str) -> Result<ParamValue> {
    i32::try_from(value)
        .map(|v| ParamValue::Integer(i64::from(v)))
        .map_err(|_| bad(format!("'{name}' does not fit the source 32-bit parameter")))
}

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

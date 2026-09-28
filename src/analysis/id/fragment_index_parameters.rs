// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` half of `FragmentIndex`: defaults, validation and
//! the `updateMembers_` conversion.

use super::FragmentIndex;
use crate::chemistry::{DigestionSpecificity, ModificationsDB, Protease, ProteaseDB};
use crate::concept::constants::user_param;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

/// The source `DefaultParamHandler` name.
pub const FRAGMENT_INDEX_NAME: &str = "FragmentIndex";

/// Parameters of [`FragmentIndex`], one field per source member set by
/// `updateMembers_`, plus the parameters the source only declares so that a
/// shared `Search:` subtree does not warn.
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentIndexParameters {
    /// Index b ions (`ions:add_b_ions`, true).
    pub add_b_ions: bool,
    /// Index y ions (`ions:add_y_ions`, true).
    pub add_y_ions: bool,
    /// Index a ions (`ions:add_a_ions`, false).
    pub add_a_ions: bool,
    /// Index c ions (`ions:add_c_ions`, false).
    pub add_c_ions: bool,
    /// Index x ions (`ions:add_x_ions`, false).
    pub add_x_ions: bool,
    /// Index z ions (`ions:add_z_ions`, false).
    pub add_z_ions: bool,
    /// Lower-side precursor tolerance, a nonnegative magnitude
    /// (`precursor:mass_tolerance_lower`, 20); the window is
    /// `[-lower, +upper]`.
    pub precursor_mass_tolerance_lower: f64,
    /// Upper-side precursor tolerance magnitude
    /// (`precursor:mass_tolerance_upper`, 20).
    pub precursor_mass_tolerance_upper: f64,
    /// Precursor tolerance in ppm, else Da
    /// (`precursor:mass_tolerance_unit`, ppm).
    pub precursor_mass_tolerance_unit_ppm: bool,
    /// Fragment tolerance (`fragment:mass_tolerance`, 10), kept in `f32` as
    /// the source member is.
    pub fragment_mz_tolerance: f32,
    /// Fragment tolerance in ppm, else Da (`fragment:mass_tolerance_unit`,
    /// ppm).
    pub fragment_mz_tolerance_unit_ppm: bool,
    /// Smallest precursor charge tried when a spectrum's is unknown
    /// (`precursor:min_charge`, 2).
    pub min_precursor_charge: u16,
    /// Largest precursor charge tried (`precursor:max_charge`, 5).
    pub max_precursor_charge: u16,
    /// Smallest indexed fragment m/z (`fragment:min_mz`, 150; an integer
    /// parameter).
    pub fragment_min_mz: i32,
    /// Largest indexed fragment m/z (`fragment:max_mz`, 2000).
    pub fragment_max_mz: i32,
    /// Ions with index at most this are not indexed (`fragment:min_ion_index`,
    /// 2: skips b1/b2/y1/y2; 0 indexes all).
    pub min_ion_index: usize,
    /// UniMod names of fixed modifications (`modifications:fixed`,
    /// `["Carbamidomethyl (C)"]`).
    pub modifications_fixed: Vec<String>,
    /// UniMod names of variable modifications (`modifications:variable`,
    /// `["Oxidation (M)"]`).
    pub modifications_variable: Vec<String>,
    /// Maximum variable modifications per peptide
    /// (`modifications:variable_max_per_peptide`, 2).
    pub modifications_max_variable_mods_per_peptide: usize,
    /// Digestion enzyme (`enzyme`, Trypsin).
    pub enzyme: String,
    /// Missed cleavages (`peptide:missed_cleavages`, 1).
    pub missed_cleavages: usize,
    /// Enzyme specificity (`peptide:enzyme_specificity`, full).
    pub enzyme_specificity: DigestionSpecificity,
    /// Minimal peptide length (`peptide:min_size`, 7).
    pub peptide_min_length: usize,
    /// Maximal peptide length, 0 = unlimited (`peptide:max_size`, 40).
    pub peptide_max_length: usize,
    /// Minimal peptide (M+H)+ (`peptide:min_mass`, 100; an integer parameter).
    pub peptide_min_mass: i32,
    /// Maximal peptide (M+H)+ (`peptide:max_mass`, 9000).
    pub peptide_max_mass: i32,
    /// Candidates with fewer matched peaks are discarded
    /// (`fragment:min_matched_ions`, 5).
    pub min_matched_peaks: u16,
    /// Smallest isotope error of a closed search
    /// (`precursor:isotope_error_min`, 0).
    pub min_isotope_error: i16,
    /// Largest isotope error of a closed search
    /// (`precursor:isotope_error_max`, 2).
    pub max_isotope_error: i16,
    /// Opt into SNES for `peptide:enzyme_specificity` `none` (`snes_enabled`,
    /// false).
    pub snes_enabled: bool,
    /// Largest fragment charge queried (`fragment:max_charge`, 2).
    pub max_fragment_charge: u16,
    /// Candidates kept per spectrum (`scoring:max_candidates_per_spectrum`,
    /// 50).
    pub max_processed_hits: u32,
    /// Declared only (`decoys`, false).
    pub decoys: bool,
    /// Declared only (`annotate:PSM`, `["ALL"]`).
    pub annotate_psm: Vec<String>,
    /// Declared only (`report:top_hits`, 1).
    pub report_top_hits: i64,
    /// Declared only (`peptide:motif`, empty).
    pub peptide_motif: String,
    /// Declared only (`precursor:isotopes`, `[0, 1]`).
    pub precursor_isotopes: Vec<i32>,
}

impl Default for FragmentIndexParameters {
    fn default() -> Self {
        Self {
            add_b_ions: true,
            add_y_ions: true,
            add_a_ions: false,
            add_c_ions: false,
            add_x_ions: false,
            add_z_ions: false,
            precursor_mass_tolerance_lower: 20.0,
            precursor_mass_tolerance_upper: 20.0,
            precursor_mass_tolerance_unit_ppm: true,
            fragment_mz_tolerance: 10.0,
            fragment_mz_tolerance_unit_ppm: true,
            min_precursor_charge: 2,
            max_precursor_charge: 5,
            fragment_min_mz: 150,
            fragment_max_mz: 2000,
            min_ion_index: 2,
            modifications_fixed: vec!["Carbamidomethyl (C)".into()],
            modifications_variable: vec!["Oxidation (M)".into()],
            modifications_max_variable_mods_per_peptide: 2,
            enzyme: "Trypsin".into(),
            missed_cleavages: 1,
            enzyme_specificity: DigestionSpecificity::Full,
            peptide_min_length: 7,
            peptide_max_length: 40,
            peptide_min_mass: 100,
            peptide_max_mass: 9000,
            min_matched_peaks: 5,
            min_isotope_error: 0,
            max_isotope_error: 2,
            snes_enabled: false,
            max_fragment_charge: 2,
            max_processed_hits: 50,
            decoys: false,
            annotate_psm: vec!["ALL".into()],
            report_top_hits: 1,
            peptide_motif: String::new(),
            precursor_isotopes: vec![0, 1],
        }
    }
}

const SPECIFICITY_DESCRIPTION: &str = "Enzyme cleavage specificity required for both peptide termini.\n  'full' : both termini must be enzyme-specific (canonical, e.g. tryptic).\n  'semi' : only one terminus needs to be enzyme-specific (semi-tryptic).\n  'none' : no enzyme constraint at either terminus; every substring of length\n           [min_size, max_size] is enumerated. This is the canonical setting for\n           immunopeptidomics (e.g. HLA peptides 8..12mers). For very large search\n           spaces consider tightening 'peptide:min_size'/'peptide:max_size'.";

const SNES_DESCRIPTION: &str = "[experimental, v1 opt-in] When peptide:enzyme_specificity=none, use mother-peptide indexing (Single-N + Single-C) instead of naïve O(L^2) sub-peptide enumeration. Orders-of-magnitude smaller index and faster search on non-specific workloads (immunopeptidomics). Ignored for specific/semi-specific enzymes. Variable modifications are applied via query-time subset enumeration (v1.1).";

const ANNOTATIONS: [&str; 15] = [
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
    user_param::HYPERSCORE_ZSCORE,
    user_param::LN_NUM_CANDIDATES,
    user_param::MATCHED_ION_CURRENT_FRACTION,
    user_param::COMPLEMENTARY_IONS_FRACTION,
];

impl FragmentIndexParameters {
    /// The source defaults with descriptions and restrictions, as the
    /// `FragmentIndex` constructor installs them (in its order). The
    /// fragment-unit description is the source's truncated
    /// `"Unit of fragment m"`, kept verbatim.
    ///
    /// # Errors
    ///
    /// Propagates [`Param`] construction errors, which the fixed defaults do
    /// not trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// These parameters in the layout of [`FragmentIndexParameters::defaults`]
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a count does not fit the source's
    /// 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut p = Param::new();
        let tags: &[String] = &[];
        let flags = ["true".to_owned(), "false".to_owned()];
        let units = ["ppm".to_owned(), "Da".to_owned()];
        let flag = |v: bool| ParamValue::String(if v { "true" } else { "false" }.into());
        for (key, value, description) in [
            (
                "ions:add_y_ions",
                self.add_y_ions,
                "Add peaks of y-ions to the spectrum",
            ),
            (
                "ions:add_b_ions",
                self.add_b_ions,
                "Add peaks of b-ions to the spectrum",
            ),
            (
                "ions:add_a_ions",
                self.add_a_ions,
                "Add peaks of a-ions to the spectrum",
            ),
            (
                "ions:add_c_ions",
                self.add_c_ions,
                "Add peaks of c-ions to the spectrum",
            ),
            (
                "ions:add_x_ions",
                self.add_x_ions,
                "Add peaks of  x-ions to the spectrum",
            ),
            (
                "ions:add_z_ions",
                self.add_z_ions,
                "Add peaks of z-ions to the spectrum",
            ),
        ] {
            p.set_value(key, flag(value), description, tags)?;
            p.set_valid_strings(key, &flags)?;
        }
        p.set_section_description("ions", "Theoretical ion series toggles")?;
        p.set_value(
            "precursor:mass_tolerance_lower",
            ParamValue::Float(self.precursor_mass_tolerance_lower),
            "Lower-side precursor-mass tolerance (positive magnitude; effective window is [-lower, +upper] around the precursor). When strongly asymmetric, also review precursor:isotope_error_min.",
            tags,
        )?;
        p.set_min_float("precursor:mass_tolerance_lower", 0.0)?;
        p.set_value(
            "precursor:mass_tolerance_upper",
            ParamValue::Float(self.precursor_mass_tolerance_upper),
            "Upper-side precursor-mass tolerance (positive magnitude).",
            tags,
        )?;
        p.set_min_float("precursor:mass_tolerance_upper", 0.0)?;
        p.set_value(
            "precursor:mass_tolerance_unit",
            ParamValue::String(unit(self.precursor_mass_tolerance_unit_ppm).into()),
            "Unit of precursor mass tolerance.",
            tags,
        )?;
        p.set_valid_strings("precursor:mass_tolerance_unit", &units)?;
        p.set_value(
            "fragment:mass_tolerance",
            ParamValue::Float(f64::from(self.fragment_mz_tolerance)),
            "Fragment mass tolerance",
            tags,
        )?;
        p.set_value(
            "fragment:mass_tolerance_unit",
            ParamValue::String(unit(self.fragment_mz_tolerance_unit_ppm).into()),
            "Unit of fragment m",
            tags,
        )?;
        p.set_valid_strings("fragment:mass_tolerance_unit", &units)?;
        p.set_value(
            "precursor:min_charge",
            ParamValue::Integer(i64::from(self.min_precursor_charge)),
            "min precursor charge",
            tags,
        )?;
        p.set_value(
            "precursor:max_charge",
            ParamValue::Integer(i64::from(self.max_precursor_charge)),
            "max precursor charge",
            tags,
        )?;
        p.set_value(
            "fragment:min_mz",
            ParamValue::Integer(i64::from(self.fragment_min_mz)),
            "Minimal fragment mz for database",
            tags,
        )?;
        p.set_value(
            "fragment:max_mz",
            ParamValue::Integer(i64::from(self.fragment_max_mz)),
            "Maximal fragment mz for database",
            tags,
        )?;
        p.set_value(
            "fragment:min_ion_index",
            int(self.min_ion_index, "min_ion_index")?,
            "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.",
            tags,
        )?;
        p.set_min_int("fragment:min_ion_index", 0)?;
        let all_mods = all_search_modifications(ModificationsDB::global());
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
            int(
                self.modifications_max_variable_mods_per_peptide,
                "variable_max_per_peptide",
            )?,
            "Maximum number of residues carrying a variable modification per candidate peptide",
            tags,
        )?;
        let enzymes: Vec<String> = ProteaseDB::global()
            .names()
            .into_iter()
            .map(str::to_owned)
            .collect();
        p.set_value(
            "enzyme",
            ParamValue::String(self.enzyme.clone()),
            "Enzyme for digestion",
            tags,
        )?;
        p.set_valid_strings("enzyme", &enzymes)?;
        p.set_value(
            "peptide:missed_cleavages",
            int(self.missed_cleavages, "missed_cleavages")?,
            "Missed cleavages for digestion",
            tags,
        )?;
        p.set_value(
            "peptide:enzyme_specificity",
            ParamValue::String(self.enzyme_specificity.name().into()),
            SPECIFICITY_DESCRIPTION,
            tags,
        )?;
        p.set_valid_strings(
            "peptide:enzyme_specificity",
            &["full".into(), "semi".into(), "none".into()],
        )?;
        p.set_value(
            "peptide:min_size",
            int(self.peptide_min_length, "min_size")?,
            "Minimal peptide length for database",
            tags,
        )?;
        p.set_value(
            "peptide:max_size",
            int(self.peptide_max_length, "max_size")?,
            "Maximal peptide length for database",
            tags,
        )?;
        p.set_value(
            "peptide:min_mass",
            ParamValue::Integer(i64::from(self.peptide_min_mass)),
            "Minimal peptide mass for database",
            tags,
        )?;
        p.set_value(
            "peptide:max_mass",
            ParamValue::Integer(i64::from(self.peptide_max_mass)),
            "Maximal peptide mass for database",
            tags,
        )?;
        p.set_value(
            "fragment:min_matched_ions",
            ParamValue::Integer(i64::from(self.min_matched_peaks)),
            "Minimal number of matched ions to report a PSM",
            tags,
        )?;
        p.set_value(
            "precursor:isotope_error_min",
            ParamValue::Integer(i64::from(self.min_isotope_error)),
            "Minimum allowed precursor isotope error",
            tags,
        )?;
        p.set_value(
            "precursor:isotope_error_max",
            ParamValue::Integer(i64::from(self.max_isotope_error)),
            "Maximum allowed precursor isotope error",
            tags,
        )?;
        p.set_value(
            "snes_enabled",
            flag(self.snes_enabled),
            SNES_DESCRIPTION,
            tags,
        )?;
        p.set_valid_strings("snes_enabled", &flags)?;
        p.set_value(
            "fragment:max_charge",
            ParamValue::Integer(i64::from(self.max_fragment_charge)),
            "max fragment charge",
            tags,
        )?;
        p.set_value(
            "scoring:max_candidates_per_spectrum",
            ParamValue::Integer(i64::from(self.max_processed_hits)),
            "The number of initial hits for which we calculate a score",
            tags,
        )?;
        p.set_section_description("scoring", "Search/Scoring Limits")?;
        p.set_value(
            "decoys",
            flag(self.decoys),
            "Should decoys be generated?",
            tags,
        )?;
        p.set_valid_strings("decoys", &flags)?;
        p.set_value(
            "annotate:PSM",
            ParamValue::StringList(self.annotate_psm.clone()),
            "Annotations added to each PSM.",
            tags,
        )?;
        let annotations: Vec<String> = ANNOTATIONS.iter().map(|&s| s.to_owned()).collect();
        p.set_valid_strings("annotate:PSM", &annotations)?;
        p.set_value(
            "report:top_hits",
            ParamValue::Integer(self.report_top_hits),
            "Maximum number of top scoring hits per spectrum that are reported.",
            tags,
        )?;
        p.set_section_description("report", "Reporting Options")?;
        p.set_value(
            "peptide:motif",
            ParamValue::String(self.peptide_motif.clone()),
            "If set, only peptides that contain this motif (provided as RegEx) will be considered.",
            tags,
        )?;
        p.set_section_description("peptide", "Peptide Options")?;
        p.set_value(
            "precursor:isotopes",
            ParamValue::IntegerList(self.precursor_isotopes.clone()),
            "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)",
            tags,
        )?;
        Ok(p)
    }

    /// Parameters from a [`Param`], as source `setParameters` followed by
    /// `updateMembers_`: missing entries take their defaults, unknown ones
    /// are ignored (see [`FragmentIndexParameters::from_param_with_warnings`]),
    /// then [`FragmentIndexParameters::validate`] runs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a wrong type, a violated
    /// restriction or a failed validation (the source's
    /// `Exception::InvalidParameter`), and for an integer outside its member
    /// type (a negative count or charge, a charge above 65535, an isotope
    /// error outside `i16`); the source converts those silently and wraps.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`FragmentIndexParameters::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`FragmentIndexParameters::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(FRAGMENT_INDEX_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (params, warnings) =
            handler.set_parameters_with(parameters, Self::from_complete_param)?;
        params.validate()?;
        Ok((params, warnings))
    }

    fn from_complete_param(p: &Param) -> Result<Self> {
        let flag = |key: &str| -> Result<bool> { Ok(p.value(key)?.as_str()? == "true") };
        let float = |key: &str| -> Result<f64> { p.value(key)?.to_f64() };
        let integer = |key: &str| -> Result<i64> { p.value(key)?.to_i64() };
        fn narrow<T: TryFrom<i64>>(value: i64, key: &str) -> Result<T> {
            T::try_from(value).map_err(|_| Error::InvalidValue(format!("{key} is out of range")))
        }
        Ok(Self {
            add_b_ions: flag("ions:add_b_ions")?,
            add_y_ions: flag("ions:add_y_ions")?,
            add_a_ions: flag("ions:add_a_ions")?,
            add_c_ions: flag("ions:add_c_ions")?,
            add_x_ions: flag("ions:add_x_ions")?,
            add_z_ions: flag("ions:add_z_ions")?,
            precursor_mass_tolerance_lower: float("precursor:mass_tolerance_lower")?,
            precursor_mass_tolerance_upper: float("precursor:mass_tolerance_upper")?,
            precursor_mass_tolerance_unit_ppm: p
                .value("precursor:mass_tolerance_unit")?
                .as_str()?
                == "ppm",
            fragment_mz_tolerance: float("fragment:mass_tolerance")? as f32,
            fragment_mz_tolerance_unit_ppm: p.value("fragment:mass_tolerance_unit")?.as_str()?
                == "ppm",
            min_precursor_charge: narrow(integer("precursor:min_charge")?, "precursor:min_charge")?,
            max_precursor_charge: narrow(integer("precursor:max_charge")?, "precursor:max_charge")?,
            fragment_min_mz: narrow(integer("fragment:min_mz")?, "fragment:min_mz")?,
            fragment_max_mz: narrow(integer("fragment:max_mz")?, "fragment:max_mz")?,
            min_ion_index: narrow(integer("fragment:min_ion_index")?, "fragment:min_ion_index")?,
            modifications_fixed: p.value("modifications:fixed")?.as_string_list()?.to_vec(),
            modifications_variable: p
                .value("modifications:variable")?
                .as_string_list()?
                .to_vec(),
            modifications_max_variable_mods_per_peptide: narrow(
                integer("modifications:variable_max_per_peptide")?,
                "modifications:variable_max_per_peptide",
            )?,
            enzyme: p.value("enzyme")?.as_str()?.to_owned(),
            missed_cleavages: narrow(
                integer("peptide:missed_cleavages")?,
                "peptide:missed_cleavages",
            )?,
            enzyme_specificity: p.value("peptide:enzyme_specificity")?.as_str()?.parse()?,
            peptide_min_length: narrow(integer("peptide:min_size")?, "peptide:min_size")?,
            peptide_max_length: narrow(integer("peptide:max_size")?, "peptide:max_size")?,
            peptide_min_mass: narrow(integer("peptide:min_mass")?, "peptide:min_mass")?,
            peptide_max_mass: narrow(integer("peptide:max_mass")?, "peptide:max_mass")?,
            min_matched_peaks: narrow(
                integer("fragment:min_matched_ions")?,
                "fragment:min_matched_ions",
            )?,
            min_isotope_error: narrow(
                integer("precursor:isotope_error_min")?,
                "precursor:isotope_error_min",
            )?,
            max_isotope_error: narrow(
                integer("precursor:isotope_error_max")?,
                "precursor:isotope_error_max",
            )?,
            snes_enabled: flag("snes_enabled")?,
            max_fragment_charge: narrow(integer("fragment:max_charge")?, "fragment:max_charge")?,
            max_processed_hits: narrow(
                integer("scoring:max_candidates_per_spectrum")?,
                "scoring:max_candidates_per_spectrum",
            )?,
            decoys: flag("decoys")?,
            annotate_psm: p.value("annotate:PSM")?.as_string_list()?.to_vec(),
            report_top_hits: integer("report:top_hits")?,
            peptide_motif: p.value("peptide:motif")?.as_str()?.to_owned(),
            precursor_isotopes: p.value("precursor:isotopes")?.as_integer_list()?.to_vec(),
        })
    }

    /// The checks of source `updateMembers_`, plus native ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a precursor tolerance is negative
    /// or not finite, when `lower + upper` is not positive (zero-width
    /// window), when SNES mode is active without both b and y ions (all
    /// source checks), and natively when the fragment tolerance is negative
    /// or not finite or the enzyme is unknown.
    pub fn validate(&self) -> Result<()> {
        let (lower, upper) = (
            self.precursor_mass_tolerance_lower,
            self.precursor_mass_tolerance_upper,
        );
        if !lower.is_finite() || !upper.is_finite() {
            return Err(Error::InvalidValue(
                "precursor:mass_tolerance_lower and mass_tolerance_upper must be finite".into(),
            ));
        }
        if lower < 0.0 || upper < 0.0 {
            return Err(Error::InvalidValue(
                "precursor mass tolerance magnitudes must not be negative".into(),
            ));
        }
        if lower + upper <= 0.0 {
            return Err(Error::InvalidValue(
                "precursor window has zero width (lower + upper must be > 0)".into(),
            ));
        }
        if !self.fragment_mz_tolerance.is_finite() || self.fragment_mz_tolerance < 0.0 {
            return Err(Error::InvalidValue(
                "fragment mass tolerance must be finite and nonnegative".into(),
            ));
        }
        if self.is_snes_mode() && (!self.add_b_ions || !self.add_y_ions) {
            return Err(Error::InvalidValue(
                "SNES mode (snes_enabled=true with enzyme_specificity=none) requires both ions:add_b_ions=true and ions:add_y_ions=true in v1. Additional ion series (a/c/x/z) may be enabled freely for downstream scoring.".into(),
            ));
        }
        Protease::from_name(&self.enzyme)?;
        Ok(())
    }

    /// Whether these tolerances trigger open-search mode (source private
    /// `isOpenSearchMode_`).
    pub fn is_open_search_mode(&self) -> bool {
        FragmentIndex::is_open_search_mode(
            self.precursor_mass_tolerance_lower,
            self.precursor_mass_tolerance_upper,
            self.precursor_mass_tolerance_unit_ppm,
        )
    }

    /// `snes_enabled` together with `none` specificity.
    pub fn is_snes_mode(&self) -> bool {
        self.snes_enabled && self.enzyme_specificity == DigestionSpecificity::None
    }
}

fn unit(ppm: bool) -> &'static str {
    if ppm { "ppm" } else { "Da" }
}

fn int(value: usize, name: &str) -> Result<ParamValue> {
    i32::try_from(value)
        .map(|v| ParamValue::Integer(i64::from(v)))
        .map_err(|_| {
            Error::InvalidValue(format!("'{name}' does not fit the source 32-bit parameter"))
        })
}

/// Source `ModificationsDB::getAllSearchModifications`: full IDs of records
/// with a UniMod record ID, sorted case-insensitively, shorter first.
fn all_search_modifications(database: &ModificationsDB) -> Vec<String> {
    let mut names: Vec<String> = database
        .entries()
        .iter()
        .filter(|record| record.record_id().is_some_and(|id| id > 0))
        .map(|record| record.full_id().to_owned())
        .collect();
    names.sort_by_cached_key(|name| (name.to_lowercase(), name.len()));
    names
}

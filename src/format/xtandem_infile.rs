// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! X! Tandem input (configuration) file writer.
//!
//! Ports `FORMAT/XTandemInfile.h`/`.cpp`. X! Tandem reads its search
//! parameters from `note` elements of a `bioml` document, e.g.
//!
//! ```text
//! <note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
//! <note type="input" label="output, proteins">yes</note>
//! ```
//!
//! The source class derives from `Internal::XMLFile` but only writes: it has
//! no `load`. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::chemistry::{ModificationDefinition, ModificationDefinitionsSet, TermSpecificity};
use crate::format::file_info::text_format::{ostream_g, to_str};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// Unit of a mass tolerance (source `XTandemInfile::ErrorUnit`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ErrorUnit {
    /// Daltons (source `DALTONS = 0`).
    #[default]
    Daltons,
    /// Parts per million (source `PPM`).
    Ppm,
}

/// Mass type of the precursor (source `XTandemInfile::MassType`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum MassType {
    /// Monoisotopic masses (source `MONOISOTOPIC = 0`).
    #[default]
    Monoisotopic,
    /// Average masses (source `AVERAGE`).
    Average,
}

/// The result types X! Tandem accepts for `output, results`.
pub const OUTPUT_RESULTS: [&str; 3] = ["valid", "all", "stochastic"];

/// X! Tandem input file writer.
///
/// Holds the search settings and writes them as an X! Tandem `bioml`
/// configuration with [`XTandemInfile::write`]. The defaults are the
/// source constructor's: fragment tolerance 0.3 Da, precursor tolerance
/// +2.0/-2.0 Da, monoisotopic masses, maximum precursor charge 4, one
/// thread, cleavage site `[KR]|{P}`, one missed cleavage, output results
/// `valid` and a maximum E-value of 0.01.
#[derive(Clone, Debug, PartialEq)]
pub struct XTandemInfile {
    fragment_mass_tolerance: f64,
    precursor_mass_tolerance_plus: f64,
    precursor_mass_tolerance_minus: f64,
    fragment_mass_error_unit: ErrorUnit,
    precursor_mass_error_unit: ErrorUnit,
    fragment_mass_type: MassType,
    precursor_mass_type: MassType,
    max_precursor_charge: i32,
    number_of_threads: u32,
    modifications: ModificationDefinitionsSet,
    input_filename: String,
    output_filename: String,
    taxonomy_file: String,
    taxon: String,
    cleavage_site: String,
    semi_cleavage: bool,
    allow_isotope_error: bool,
    number_of_missed_cleavages: u32,
    default_parameters_file: String,
    output_results: String,
    max_valid_evalue: f64,
}

impl Default for XTandemInfile {
    fn default() -> Self {
        Self {
            fragment_mass_tolerance: 0.3,
            precursor_mass_tolerance_plus: 2.0,
            precursor_mass_tolerance_minus: 2.0,
            fragment_mass_error_unit: ErrorUnit::Daltons,
            precursor_mass_error_unit: ErrorUnit::Daltons,
            fragment_mass_type: MassType::Monoisotopic,
            precursor_mass_type: MassType::Monoisotopic,
            max_precursor_charge: 4,
            number_of_threads: 1,
            modifications: ModificationDefinitionsSet::default(),
            input_filename: String::new(),
            output_filename: String::new(),
            taxonomy_file: String::new(),
            taxon: String::new(),
            cleavage_site: "[KR]|{P}".to_owned(),
            semi_cleavage: false,
            allow_isotope_error: false,
            number_of_missed_cleavages: 1,
            default_parameters_file: String::new(),
            output_results: "valid".to_owned(),
            max_valid_evalue: 0.01,
        }
    }
}

impl XTandemInfile {
    /// A writer with the source constructor's defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the fragment mass tolerance.
    pub fn set_fragment_mass_tolerance(&mut self, tolerance: f64) {
        self.fragment_mass_tolerance = tolerance;
    }

    /// The fragment mass tolerance.
    pub fn fragment_mass_tolerance(&self) -> f64 {
        self.fragment_mass_tolerance
    }

    /// Set the precursor mass tolerance above the precursor mass.
    pub fn set_precursor_mass_tolerance_plus(&mut self, tolerance: f64) {
        self.precursor_mass_tolerance_plus = tolerance;
    }

    /// The precursor mass tolerance above the precursor mass.
    pub fn precursor_mass_tolerance_plus(&self) -> f64 {
        self.precursor_mass_tolerance_plus
    }

    /// Set the precursor mass tolerance below the precursor mass.
    pub fn set_precursor_mass_tolerance_minus(&mut self, tolerance: f64) {
        self.precursor_mass_tolerance_minus = tolerance;
    }

    /// The precursor mass tolerance below the precursor mass.
    pub fn precursor_mass_tolerance_minus(&self) -> f64 {
        self.precursor_mass_tolerance_minus
    }

    /// Set the precursor mass type (source `setPrecursorErrorType`). It is
    /// stored only: the written file has no precursor mass type note.
    pub fn set_precursor_error_type(&mut self, mass_type: MassType) {
        self.precursor_mass_type = mass_type;
    }

    /// The precursor mass type (source `getPrecursorErrorType`).
    pub fn precursor_error_type(&self) -> MassType {
        self.precursor_mass_type
    }

    /// Set the unit of the fragment mass tolerance.
    pub fn set_fragment_mass_error_unit(&mut self, unit: ErrorUnit) {
        self.fragment_mass_error_unit = unit;
    }

    /// The unit of the fragment mass tolerance.
    pub fn fragment_mass_error_unit(&self) -> ErrorUnit {
        self.fragment_mass_error_unit
    }

    /// Set the unit of the precursor mass tolerance.
    pub fn set_precursor_mass_error_unit(&mut self, unit: ErrorUnit) {
        self.precursor_mass_error_unit = unit;
    }

    /// The unit of the precursor mass tolerance.
    pub fn precursor_mass_error_unit(&self) -> ErrorUnit {
        self.precursor_mass_error_unit
    }

    /// Set the number of threads X! Tandem uses.
    pub fn set_number_of_threads(&mut self, threads: u32) {
        self.number_of_threads = threads;
    }

    /// The number of threads X! Tandem uses.
    pub fn number_of_threads(&self) -> u32 {
        self.number_of_threads
    }

    /// Set the fixed and variable modifications.
    pub fn set_modifications(&mut self, mods: ModificationDefinitionsSet) {
        self.modifications = mods;
    }

    /// The fixed and variable modifications.
    pub fn modifications(&self) -> &ModificationDefinitionsSet {
        &self.modifications
    }

    /// Set the path of the X! Tandem result file (`output, path`).
    pub fn set_output_filename(&mut self, output: &str) {
        self.output_filename = output.to_owned();
    }

    /// The path of the X! Tandem result file.
    pub fn output_filename(&self) -> &str {
        &self.output_filename
    }

    /// Set the path of the spectrum file (`spectrum, path`).
    pub fn set_input_filename(&mut self, input_file: &str) {
        self.input_filename = input_file.to_owned();
    }

    /// The path of the spectrum file.
    pub fn input_filename(&self) -> &str {
        &self.input_filename
    }

    /// Set the path of the taxonomy file, which names the sequence database
    /// (`list path, taxonomy information`).
    pub fn set_taxonomy_filename(&mut self, filename: &str) {
        self.taxonomy_file = filename.to_owned();
    }

    /// The path of the taxonomy file.
    pub fn taxonomy_filename(&self) -> &str {
        &self.taxonomy_file
    }

    /// Set the path of the default parameters file
    /// (`list path, default parameters`); an empty path writes no note.
    pub fn set_default_parameters_filename(&mut self, filename: &str) {
        self.default_parameters_file = filename.to_owned();
    }

    /// The path of the default parameters file.
    pub fn default_parameters_filename(&self) -> &str {
        &self.default_parameters_file
    }

    /// Set the taxon looked up in the taxonomy file (`protein, taxon`).
    pub fn set_taxon(&mut self, taxon: &str) {
        self.taxon = taxon.to_owned();
    }

    /// The taxon looked up in the taxonomy file.
    pub fn taxon(&self) -> &str {
        &self.taxon
    }

    /// Set the maximum precursor charge.
    pub fn set_max_precursor_charge(&mut self, max_charge: i32) {
        self.max_precursor_charge = max_charge;
    }

    /// The maximum precursor charge.
    pub fn max_precursor_charge(&self) -> i32 {
        self.max_precursor_charge
    }

    /// Set the number of missed cleavages allowed.
    pub fn set_number_of_missed_cleavages(&mut self, missed_cleavages: u32) {
        self.number_of_missed_cleavages = missed_cleavages;
    }

    /// The number of missed cleavages allowed.
    pub fn number_of_missed_cleavages(&self) -> u32 {
        self.number_of_missed_cleavages
    }

    /// Set the output result type, one of [`OUTPUT_RESULTS`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other value (source
    /// `Exception::FailedAPICall`); the setting is then unchanged.
    pub fn set_output_results(&mut self, result: &str) -> Result<()> {
        if OUTPUT_RESULTS.contains(&result) {
            self.output_results = result.to_owned();
            Ok(())
        } else {
            Err(Error::InvalidValue(format!(
                "Invalid result type provided (must be either all, valid or stochastic).: '{result}'"
            )))
        }
    }

    /// The output result type (`valid`, `all` or `stochastic`).
    pub fn output_results(&self) -> &str {
        &self.output_results
    }

    /// Set the largest E-value a reported hit may have.
    pub fn set_max_valid_evalue(&mut self, value: f64) {
        self.max_valid_evalue = value;
    }

    /// The largest E-value a reported hit may have.
    pub fn max_valid_evalue(&self) -> f64 {
        self.max_valid_evalue
    }

    /// Allow or forbid semi-specific cleavage.
    pub fn set_semi_cleavage(&mut self, semi_cleavage: bool) {
        self.semi_cleavage = semi_cleavage;
    }

    /// Whether semi-specific cleavage is allowed. The source has no getter;
    /// the port adds it for symmetry.
    pub fn semi_cleavage(&self) -> bool {
        self.semi_cleavage
    }

    /// Also consider assigning the precursor to its first or second ¹³C
    /// isotopic peak.
    pub fn set_allow_isotope_error(&mut self, allow_isotope_error: bool) {
        self.allow_isotope_error = allow_isotope_error;
    }

    /// Whether isotope errors are allowed. The source has no getter; the
    /// port adds it for symmetry.
    pub fn allow_isotope_error(&self) -> bool {
        self.allow_isotope_error
    }

    /// Set the cleavage site as an X! Tandem cleavage expression, e.g.
    /// `[KR]|{P}` for trypsin.
    pub fn set_cleavage_site(&mut self, cleavage_site: &str) {
        self.cleavage_site = cleavage_site.to_owned();
    }

    /// The cleavage site expression.
    pub fn cleavage_site(&self) -> &str {
        &self.cleavage_site
    }

    /// Write the configuration to `path`.
    ///
    /// With `ignore_member_parameters`, only the notes fed by the file names
    /// and the fixed output settings are written. With `force_default_mods`,
    /// modifications X! Tandem otherwise covers with its `quick acetyl` and
    /// `quick pyrolidone` options are written explicitly.
    ///
    /// Returns the messages the source logs as warnings or information
    /// (duplicate origins, terminal modifications with a residue, corrected
    /// variable modification masses, implicit N-terminal modifications).
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written (source
    /// `Exception::UnableToCreateFile`), and the modification lookup's error
    /// for a definition without chemistry.
    pub fn write(
        &self,
        path: impl AsRef<Path>,
        ignore_member_parameters: bool,
        force_default_mods: bool,
    ) -> Result<Vec<String>> {
        let (text, messages) = self.render(ignore_member_parameters, force_default_mods)?;
        let mut file = std::fs::File::create(path)?;
        file.write_all(text.as_bytes())?;
        file.flush()?;
        Ok(messages)
    }

    /// Write the configuration to `writer`; see [`XTandemInfile::write`].
    ///
    /// # Errors
    ///
    /// As [`XTandemInfile::write`].
    pub fn write_to(
        &self,
        mut writer: impl Write,
        ignore_member_parameters: bool,
        force_default_mods: bool,
    ) -> Result<Vec<String>> {
        let (text, messages) = self.render(ignore_member_parameters, force_default_mods)?;
        writer.write_all(text.as_bytes())?;
        Ok(messages)
    }

    /// The configuration text and the log messages (source `writeTo_`).
    fn render(
        &self,
        ignore_member_parameters: bool,
        force_default_mods: bool,
    ) -> Result<(String, Vec<String>)> {
        let mut os = String::new();
        let mut messages = Vec::new();
        os.push_str("<?xml version=\"1.0\"?>\n");
        os.push_str("<?xml-stylesheet type=\"text/xsl\" href=\"tandem-input-style.xsl\"?>\n");
        os.push_str("<bioml>\n");

        note(&mut os, "spectrum, path", &self.input_filename);
        note(&mut os, "output, path", &self.output_filename);
        note(
            &mut os,
            "list path, taxonomy information",
            &self.taxonomy_file,
        );
        if !self.default_parameters_file.is_empty() {
            note(
                &mut os,
                "list path, default parameters",
                &self.default_parameters_file,
            );
        }
        // Needed for finding and parsing the results.
        note(&mut os, "output, path hashing", yes_no(false));
        note(&mut os, "output, proteins", yes_no(true));
        note(&mut os, "output, spectra", yes_no(true));
        note(&mut os, "output, sort results by", "spectrum");
        // Required by Percolator to recognise the output file.
        note(&mut os, "output, xsl path", "tandem-style.xsl");
        note(&mut os, "output, parameters", yes_no(true));

        if !ignore_member_parameters {
            let mut force_default_mods = force_default_mods;
            note(
                &mut os,
                "spectrum, fragment monoisotopic mass error",
                &to_str(self.fragment_mass_tolerance),
            );
            note(
                &mut os,
                "spectrum, parent monoisotopic mass error plus",
                &to_str(self.precursor_mass_tolerance_plus),
            );
            note(
                &mut os,
                "spectrum, parent monoisotopic mass error minus",
                &to_str(self.precursor_mass_tolerance_minus),
            );
            note(
                &mut os,
                "spectrum, parent monoisotopic mass isotope error",
                yes_no(self.allow_isotope_error),
            );
            note(
                &mut os,
                "spectrum, fragment monoisotopic mass error units",
                match self.fragment_mass_error_unit {
                    ErrorUnit::Daltons => "Daltons",
                    ErrorUnit::Ppm => "ppm",
                },
            );
            note(
                &mut os,
                "spectrum, parent monoisotopic mass error units",
                match self.precursor_mass_error_unit {
                    ErrorUnit::Ppm => "ppm",
                    ErrorUnit::Daltons => "Daltons",
                },
            );
            note(
                &mut os,
                "spectrum, fragment mass type",
                match self.fragment_mass_type {
                    MassType::Monoisotopic => "monoisotopic",
                    MassType::Average => "average",
                },
            );
            note(
                &mut os,
                "spectrum, maximum parent charge",
                &self.max_precursor_charge.to_string(),
            );
            note(
                &mut os,
                "spectrum, threads",
                &self.number_of_threads.to_string(),
            );
            note(&mut os, "protein, taxon", &self.taxon);
            note(&mut os, "protein, cleavage site", &self.cleavage_site);
            note(
                &mut os,
                "protein, cleavage semi",
                yes_no(self.semi_cleavage),
            );

            // A variable modification of the peptide N-terminus turns the
            // quick acetyl and quick pyrolidone options off (source regex
            // " \(N-term( .)?\)$").
            let var_mods = self.modifications.variable_names();
            if var_mods.iter().any(|name| {
                is_n_term_name(name)
                    && name != ACETYL_N_TERM
                    && name != PYRO_GLU_Q
                    && name != PYRO_GLU_E
            }) {
                force_default_mods = true;
            }
            if !force_default_mods && var_mods.contains(PYRO_GLU_Q) && var_mods.contains(PYRO_GLU_E)
            {
                note(&mut os, "protein, quick pyrolidone", yes_no(true));
                messages.push(format!(
                    "Modifications '{PYRO_GLU_Q}' and '{PYRO_GLU_E}' are handled implicitly by the X! Tandem option 'protein, quick pyrolidone'. Set the 'force' flag in XTandemAdapter to force explicit inclusion of these modifications."
                ));
            }
            if !force_default_mods && var_mods.contains(ACETYL_N_TERM) {
                note(&mut os, "protein, quick acetyl", yes_no(true));
                messages.push(format!(
                    "Modification '{ACETYL_N_TERM}' is handled implicitly by the X! Tandem option 'protein, quick acetyl'. Set the 'force' flag in XTandemAdapter to force explicit inclusion of this modification."
                ));
            }

            let mut affected_origins = BTreeMap::new();
            let fixed = convert_modification_set(
                &self.modifications.fixed_modifications().collect::<Vec<_>>(),
                &mut affected_origins,
                force_default_mods,
                &mut messages,
            )?;
            note(&mut os, "residue, modification mass", &fixed);
            let variable = convert_modification_set(
                &self
                    .modifications
                    .variable_modifications()
                    .collect::<Vec<_>>(),
                &mut affected_origins,
                force_default_mods,
                &mut messages,
            )?;
            note(&mut os, "residue, potential modification mass", &variable);

            note(
                &mut os,
                "scoring, maximum missed cleavage sites",
                &self.number_of_missed_cleavages.to_string(),
            );
            note(&mut os, "output, results", &self.output_results);
            note(
                &mut os,
                "output, maximum valid expectation value",
                &to_str(self.max_valid_evalue),
            );
        }
        os.push_str("</bioml>\n");
        Ok((os, messages))
    }
}

const ACETYL_N_TERM: &str = "Acetyl (N-term)";
const PYRO_GLU_Q: &str = "Gln->pyro-Glu (N-term Q)";
const PYRO_GLU_E: &str = "Glu->pyro-Glu (N-term E)";

/// Source `boost::regex_search(name, " \\(N-term( .)?\\)$")`: the name ends
/// with ` (N-term)` or ` (N-term X)` for any single character `X` other than
/// a line break.
fn is_n_term_name(name: &str) -> bool {
    let Some(body) = name.strip_suffix(')') else {
        return false;
    };
    if body.ends_with(" (N-term") {
        return true;
    }
    let mut chars = body.chars();
    match chars.next_back() {
        Some(c) if c != '\n' && c != '\r' => chars.as_str().ends_with(" (N-term "),
        _ => false,
    }
}

/// Source `writeNote_`: values are written as given, without XML escaping.
fn note(os: &mut String, label: &str, value: &str) {
    os.push_str("\t<note type=\"input\" label=\"");
    os.push_str(label);
    os.push_str("\">");
    os.push_str(value);
    os.push_str("</note>\n");
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no"
    }
}

/// Source `convertModificationSet_`: the modifications as a comma-separated
/// X! Tandem list of `±mass@origin`, with `[`/`]` for a peptide N-/C-terminal
/// modification. A variable modification of an origin that already carries
/// a fixed one has the fixed mass subtracted. `affected_origins` holds the
/// origins of the previous call on entry and this call's on return.
fn convert_modification_set(
    mods: &[&ModificationDefinition],
    affected_origins: &mut BTreeMap<String, f64>,
    force_default_mods: bool,
    messages: &mut Vec<String>,
) -> Result<String> {
    let has_pyroglu_e = mods.iter().any(|m| m.modification_name() == PYRO_GLU_E);
    let has_pyroglu_q = mods.iter().any(|m| m.modification_name() == PYRO_GLU_Q);

    let mut origin_set: BTreeMap<String, f64> = BTreeMap::new();
    let mut xtandem_mods = Vec::new();
    for definition in mods.iter().copied() {
        let name = definition.modification_name();
        if !force_default_mods
            && (name == ACETYL_N_TERM
                || (name == PYRO_GLU_Q && has_pyroglu_e)
                || (name == PYRO_GLU_E && has_pyroglu_q))
        {
            continue;
        }
        let modification = definition.modification()?;
        let mut mod_mass = modification.diff_mono_mass();
        // The source stores 'X' as the origin of a terminal modification
        // of any residue.
        let mut orig = modification.origin().unwrap_or('X').to_string();
        let term = modification.term_specificity();
        if term != TermSpecificity::Anywhere {
            messages.push(format!(
                "Warning: X! Tandem doesn't support modifications with both residue and terminal specificity. Using only terminal specificity for modification '{name}'."
            ));
        }
        match term {
            TermSpecificity::CTerm => orig = "]".to_owned(),
            TermSpecificity::NTerm => orig = "[".to_owned(),
            _ => {}
        }
        if let Some(previous) = origin_set.get(&orig) {
            messages.push(format!(
                "X! Tandem config file: Duplicate modification assignment to origin '{orig}'. X! Tandem will ignore the first modification '{}'!",
                ostream_g(*previous, 6)
            ));
        }
        if let Some(fixed_mass) = affected_origins.get(&orig) {
            messages.push(format!(
                "X! Tandem config file: Fixed modification and variable modification to origin '{orig}' detected. Using corrected mass of {} instead of {}.",
                ostream_g(mod_mass - fixed_mass, 6),
                ostream_g(mod_mass, 6)
            ));
            mod_mass -= fixed_mass;
        }
        // std::map::insert keeps the first value of a duplicate origin.
        origin_set.entry(orig.clone()).or_insert(mod_mass);

        let mut mod_string = if mod_mass >= 0.0 {
            format!("+{}", to_str(mod_mass))
        } else {
            to_str(mod_mass)
        };
        mod_string.push('@');
        mod_string.push_str(&orig);
        xtandem_mods.push(mod_string);
    }
    *affected_origins = origin_set;
    Ok(xtandem_mods.join(","))
}

#[cfg(test)]
mod tests {
    use super::is_n_term_name;

    #[test]
    fn n_term_suffix_matches_the_source_regex() {
        assert!(is_n_term_name("Dimethyl (N-term)"));
        assert!(is_n_term_name("Gln->pyro-Glu (N-term Q)"));
        assert!(!is_n_term_name("Acetyl (Protein N-term)"));
        assert!(!is_n_term_name("Oxidation (M)"));
        assert!(!is_n_term_name("Carbamyl (N-term QQ)"));
        assert!(!is_n_term_name("Ammonium (C-term)"));
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Sequest input file: the `sequest.params` file with its `[SEQUEST]` search
//! settings, the static and differential modifications, and the
//! `[SEQUEST_ENZYME_INFO]` enzyme table.
//!
//! Ports `FORMAT/SequestInfile.h`/`.cpp`. Its modification-line parser
//! `handlePTMs` is the one of `InspectInfile`
//! ([`crate::format::inspect_infile::parse_modification_line`]). See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::data_structures::string_utils::to_float;
use crate::format::file_info::text_format::{ostream_g, to_str_f32};
use crate::format::inspect_infile::parse_modification_line;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// One enzyme of the table: cleavage direction (`0` N- to C-terminal
/// side, `1` C-terminal side), the residues it cuts after, and the residues
/// it does not cut before (`-` for none).
pub type EnzymeInfo = Vec<String>;

/// Sequest input file.
///
/// The defaults are the source constructor's, with the source's seventeen
/// standard enzymes (`AspN` … `Trypsin_Strict`) and enzyme number 0.
/// Equality compares every setting, the enzyme table and the
/// modifications, as the source's `operator==` does.
#[derive(Clone, Debug, PartialEq)]
pub struct SequestInfile {
    enzyme_info: BTreeMap<String, EnzymeInfo>,
    database: String,
    neutral_losses_for_ions: String,
    ion_series_weights: String,
    partial_sequence: String,
    sequence_header_filter: String,
    protein_mass_filter: String,
    precursor_mass_tolerance: f32,
    peak_mass_tolerance: f32,
    match_peak_tolerance: f32,
    ion_cutoff_percentage: f32,
    peptide_mass_unit: usize,
    output_lines: usize,
    enzyme_number: usize,
    max_aa_per_mod_per_peptide: usize,
    max_mods_per_peptide: usize,
    nucleotide_reading_frame: usize,
    max_internal_cleavage_sites: usize,
    match_peak_count: usize,
    match_peak_allowed_error: usize,
    show_fragment_ions: bool,
    print_duplicate_references: bool,
    remove_precursor_near_peaks: bool,
    mass_type_parent: bool,
    mass_type_fragment: bool,
    normalize_xcorr: bool,
    residues_in_upper_case: bool,
    ptm_name_residues_mass_type: BTreeMap<String, Vec<String>>,
}

/// Source `setStandardEnzymeInfo_`: name, direction, cuts after, does not
/// cut before.
const STANDARD_ENZYMES: [(&str, &str, &str, &str); 17] = [
    ("No_Enzyme", "0", "-", "-"),
    ("Trypsin_Strict", "1", "KR", "-"),
    ("Trypsin", "1", "KRLNH", "-"),
    ("Chymotrypsin", "1", "FWYL", "-"),
    ("Chymotrypsin_WYF", "1", "FWY", "-"),
    ("Clostripain", "1", "R", "-"),
    ("Cyanogen_Bromide", "1", "M", "-"),
    ("IodosoBenzoate", "1", "W", "-"),
    ("Proline_Endopept", "1", "P", "-"),
    ("GluC", "1", "E", "-"),
    ("GluC_ED", "1", "ED", "-"),
    ("LysC", "1", "K", "-"),
    ("AspN", "0", "D", "-"),
    ("AspN_DE", "0", "DE", "-"),
    ("Elastase", "1", "ALIV", "P"),
    ("Elastase/Tryp/Chymo", "1", "ALIVKRWFY", "P"),
    ("Trypsin/Chymo", "1", "KRLFWYN", "-"),
];

impl Default for SequestInfile {
    fn default() -> Self {
        Self {
            enzyme_info: STANDARD_ENZYMES
                .iter()
                .map(|(name, direction, cut, no_cut)| {
                    (
                        (*name).to_owned(),
                        vec![
                            (*direction).to_owned(),
                            (*cut).to_owned(),
                            (*no_cut).to_owned(),
                        ],
                    )
                })
                .collect(),
            database: String::new(),
            neutral_losses_for_ions: "0 1 1".to_owned(),
            ion_series_weights: "0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0".to_owned(),
            partial_sequence: String::new(),
            sequence_header_filter: String::new(),
            protein_mass_filter: "0 0".to_owned(),
            precursor_mass_tolerance: 0.0,
            peak_mass_tolerance: 0.0,
            match_peak_tolerance: 0.0,
            ion_cutoff_percentage: 0.0,
            peptide_mass_unit: 0,
            output_lines: 0,
            enzyme_number: 0,
            max_aa_per_mod_per_peptide: 0,
            max_mods_per_peptide: 0,
            nucleotide_reading_frame: 0,
            max_internal_cleavage_sites: 0,
            match_peak_count: 0,
            match_peak_allowed_error: 0,
            show_fragment_ions: true,
            print_duplicate_references: true,
            remove_precursor_near_peaks: false,
            mass_type_parent: false,
            mass_type_fragment: false,
            normalize_xcorr: false,
            residues_in_upper_case: true,
            ptm_name_residues_mass_type: BTreeMap::new(),
        }
    }
}

/// `std::ostream << float` in the default float field, precision 6.
fn stream_float(value: f32) -> String {
    ostream_g(f64::from(value), 6)
}

fn flag(value: bool) -> u8 {
    u8::from(value)
}

/// The (static and differential) mass sums of [`SequestInfile::store`].
#[derive(Default)]
struct ModificationMasses {
    dyn_n_term: f32,
    dyn_c_term: f32,
    stat_n_term: f32,
    stat_c_term: f32,
    stat_n_term_prot: f32,
    stat_c_term_prot: f32,
    stat: BTreeMap<char, f32>,
    dynamic: BTreeMap<char, f32>,
}

fn mass_of(entry: &[String]) -> Result<f32> {
    to_float(&entry[1]).map_err(Error::InvalidValue)
}

impl SequestInfile {
    /// An input file with the source constructor's defaults.
    pub fn new() -> Self {
        Self::default()
    }

    fn modification_masses(&self) -> Result<ModificationMasses> {
        let mut masses = ModificationMasses::default();
        for entry in self.ptm_name_residues_mass_type.values() {
            let (residues, kind) = (entry[0].as_str(), entry[2].as_str());
            match residues {
                "CTERM" => {
                    if kind == "OPT" {
                        masses.dyn_c_term += mass_of(entry)?;
                    }
                    if kind == "FIX" {
                        masses.stat_c_term += mass_of(entry)?;
                    }
                }
                "NTERM" => {
                    if kind == "OPT" {
                        masses.dyn_n_term += mass_of(entry)?;
                    }
                    if kind == "FIX" {
                        masses.stat_n_term += mass_of(entry)?;
                    }
                }
                // Protein-terminal modifications are static whatever their type.
                "CTERM_PROT" => masses.stat_c_term_prot += mass_of(entry)?,
                "NTERM_PROT" => masses.stat_n_term_prot += mass_of(entry)?,
                _ => {
                    let mass = mass_of(entry)?;
                    let target = if kind == "FIX" {
                        &mut masses.stat
                    } else {
                        &mut masses.dynamic
                    };
                    for residue in residues.chars() {
                        *target.entry(residue).or_insert(0.0) += mass;
                    }
                }
            }
        }
        Ok(masses)
    }

    /// The `sequest.params` text (source `store`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a stored modification mass is not a
    /// number (the source's `Exception::ConversionError`).
    pub fn to_params_string(&self) -> Result<String> {
        let masses = self.modification_masses()?;

        // Differential modifications with the same mass change share one
        // entry, in ascending mass order (source std::map<float, String>).
        let mut dyn_mods_masses: Vec<(f32, String)> = Vec::new();
        for (&residue, &mass) in &masses.dynamic {
            match dyn_mods_masses.iter_mut().find(|(m, _)| *m == mass) {
                Some((_, residues)) => residues.push(residue),
                None => {
                    let at = dyn_mods_masses
                        .iter()
                        .position(|(m, _)| mass < *m)
                        .unwrap_or(dyn_mods_masses.len());
                    dyn_mods_masses.insert(at, (mass, residue.to_string()));
                }
            }
        }
        let dyn_mods_string = if dyn_mods_masses.is_empty() {
            "0 X".to_owned()
        } else {
            dyn_mods_masses
                .iter()
                .map(|(mass, residues)| format!("{} {residues}", to_str_f32(*mass)))
                .collect::<Vec<_>>()
                .join(" ")
        };

        let mut s = String::new();
        let mut line = |text: String| {
            s.push_str(&text);
            s.push('\n');
        };
        line("[SEQUEST]".into());
        line(format!("database_name = {}", self.database));
        line(format!(
            "peptide_mass_tolerance = {}",
            stream_float(self.precursor_mass_tolerance)
        ));
        line(format!(
            "peptide_mass_units = {}; 0=amu, 1=mmu, 2=ppm",
            self.peptide_mass_unit
        ));
        line(format!(
            "ion_series = {} {};nABY ABCDVWXYZ",
            self.neutral_losses_for_ions, self.ion_series_weights
        ));
        line(format!(
            "fragment_ion_tolerance = {}",
            stream_float(self.peak_mass_tolerance)
        ));
        line(format!("num_output_lines = {}", self.output_lines));
        line(format!("num_results = {}", self.output_lines));
        line("num_description_lines = 0".into());
        line(format!(
            "show_fragment_ions = {}",
            flag(self.show_fragment_ions)
        ));
        line(format!(
            "print_duplicate_references = {}",
            flag(self.print_duplicate_references)
        ));
        line(format!("enzyme_number = {}", self.enzyme_number));
        line(format!("diff_search_options = {dyn_mods_string}"));
        line(format!(
            "term_diff_search_options = {} {}",
            stream_float(masses.dyn_n_term),
            stream_float(masses.dyn_c_term)
        ));
        line(format!(
            "remove_precursor_peak = {}",
            flag(self.remove_precursor_near_peaks)
        ));
        line(format!(
            "ion_cutoff_percentage = {}",
            stream_float(self.ion_cutoff_percentage)
        ));
        line(format!(
            "protein_mass_filter = {}",
            self.protein_mass_filter
        ));
        line(format!(
            "max_differential_AA_per_mod = {}",
            self.max_aa_per_mod_per_peptide
        ));
        line(format!(
            "max_differential_per_peptide = {}",
            self.max_mods_per_peptide
        ));
        line(format!(
            "nucleotide_reading_frame = {}; 0=protein db, 1-6, 7 = forward three, 8-reverse three, 9=all six",
            self.nucleotide_reading_frame
        ));
        line(format!(
            "mass_type_parent = {}; 0=average masses, 1=monoisotopic masses",
            flag(self.mass_type_parent)
        ));
        line(format!(
            "mass_type_fragment = {}; 0=average masses, 1=monoisotopic masses",
            flag(self.mass_type_fragment)
        ));
        line(format!("normalize_xcorr = {}", flag(self.normalize_xcorr)));
        line(format!(
            "max_internal_cleavage_sites = {}",
            self.max_internal_cleavage_sites
        ));
        line("create_output_files = 1".into());
        line(format!("partial_sequence = {}", self.partial_sequence));
        line(format!(
            "sequence_header_filter = {}",
            self.sequence_header_filter
        ));
        line(format!(
            "match_peak_count = {}; number of auto-detected peaks to try matching (max 5)",
            self.match_peak_count
        ));
        line(format!(
            "match_peak_allowed_error = {}",
            self.match_peak_allowed_error
        ));
        line(format!(
            "match_peak_tolerance = {}",
            stream_float(self.match_peak_tolerance)
        ));
        line(format!(
            "residues_in_upper_case = {}\n\n",
            flag(self.residues_in_upper_case)
        ));
        line(format!(
            "add_Nterm_peptide = {}",
            stream_float(masses.stat_n_term)
        ));
        line(format!(
            "add_Cterm_peptide = {}",
            stream_float(masses.stat_c_term)
        ));
        line(format!(
            "add_Nterm_protein = {}",
            stream_float(masses.stat_n_term_prot)
        ));
        line(format!(
            "add_Cterm_protein = {}\n",
            stream_float(masses.stat_c_term_prot)
        ));
        for (residue, label, comment) in RESIDUE_LINES {
            let mass = masses.stat.get(&residue).copied().unwrap_or(0.0);
            line(format!("{label} = {}; {comment}", stream_float(mass)));
        }
        line(String::new());
        s.push_str(&self.enzyme_info_as_string());
        Ok(s)
    }

    /// Store the `sequest.params` file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written (source
    /// `Exception::UnableToCreateFile`), and the errors of
    /// [`SequestInfile::to_params_string`]; the source creates the file
    /// before it fails there, the port writes nothing.
    pub fn store(&self, path: impl AsRef<Path>) -> Result<()> {
        let text = self.to_params_string()?;
        let mut file = std::fs::File::create(path)?;
        file.write_all(text.as_bytes())?;
        file.flush()?;
        Ok(())
    }

    /// The enzyme table as the `[SEQUEST_ENZYME_INFO]` section: one line
    /// per enzyme in name order, numbered from 0, the columns padded to the
    /// longest name and residue list plus five spaces (source
    /// `getEnzymeInfoAsString`).
    pub fn enzyme_info_as_string(&self) -> String {
        let max_name = self.enzyme_info.keys().map(String::len).max().unwrap_or(0);
        let max_cut = self
            .enzyme_info
            .values()
            .map(|info| info[1].len())
            .max()
            .unwrap_or(0);
        let mut out = String::from("[SEQUEST_ENZYME_INFO]\n");
        for (i, (name, info)) in self.enzyme_info.iter().enumerate() {
            out.push_str(&format!(
                "{i}.  {name}{}{}     {}{}{}\n",
                " ".repeat(max_name + 5 - name.len()),
                info[0],
                info[1],
                " ".repeat(max_cut + 5 - info[1].len()),
                info[2]
            ));
        }
        out
    }

    /// The enzyme table: name to [`EnzymeInfo`] (source `getEnzymeInfo_`).
    pub fn enzyme_info(&self) -> &BTreeMap<String, EnzymeInfo> {
        &self.enzyme_info
    }

    /// Add (or replace) an enzyme and select it: `name`, direction (`0`/`1`),
    /// residues cut after, residues not cut before. A cut-after list with
    /// repeated residues is replaced by its distinct residues in sorted
    /// order; one without repeats is kept as given, as in the source.
    pub fn add_enzyme_info(
        &mut self,
        name: &str,
        direction: &str,
        cuts_after: &str,
        no_cut_before: &str,
    ) {
        let mut distinct: Vec<char> = cuts_after.chars().collect();
        distinct.sort_unstable();
        distinct.dedup();
        let cuts_after = if distinct.len() == cuts_after.chars().count() {
            cuts_after.to_owned()
        } else {
            distinct.into_iter().collect()
        };
        self.enzyme_info.insert(
            name.to_owned(),
            vec![direction.to_owned(), cuts_after, no_cut_before.to_owned()],
        );
        self.enzyme_number = self
            .enzyme_info
            .keys()
            .position(|key| key == name)
            .unwrap_or(self.enzyme_info.len());
    }

    /// Select the enzyme by name. Returns 0 when it is in the table and the
    /// table's size otherwise; the enzyme number is then the table's size,
    /// which names no enzyme (source `setEnzyme`).
    pub fn set_enzyme(&mut self, enzyme_name: &str) -> usize {
        match self.enzyme_info.keys().position(|key| key == enzyme_name) {
            Some(number) => {
                self.enzyme_number = number;
                0
            }
            None => {
                self.enzyme_number = self.enzyme_info.len();
                self.enzyme_info.len()
            }
        }
    }

    /// The number of the selected enzyme in the name-ordered table.
    pub fn enzyme_number(&self) -> usize {
        self.enzyme_number
    }

    /// The name of the selected enzyme.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the enzyme number names no enzyme, after
    /// [`SequestInfile::set_enzyme`] with an unknown name; the source then
    /// dereferences the end of its map.
    pub fn enzyme_name(&self) -> Result<&str> {
        self.enzyme_info
            .keys()
            .nth(self.enzyme_number)
            .map(String::as_str)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "enzyme number {} names no enzyme of the {} in the table",
                    self.enzyme_number,
                    self.enzyme_info.len()
                ))
            })
    }

    /// Parse a modification line into [`SequestInfile::modifications`]; see
    /// [`parse_modification_line`] for the grammar. Residues `CTERM`,
    /// `NTERM`, `CTERM_PROT` and `NTERM_PROT` mark terminal modifications.
    ///
    /// # Errors
    ///
    /// As [`parse_modification_line`]; the modifications are cleared first.
    pub fn handle_ptms(
        &mut self,
        modification_line: &str,
        modifications_file: impl AsRef<Path>,
        monoisotopic: bool,
    ) -> Result<()> {
        self.ptm_name_residues_mass_type.clear();
        parse_modification_line(
            modification_line,
            modifications_file.as_ref(),
            monoisotopic,
            &mut self.ptm_name_residues_mass_type,
        )
    }

    /// The modifications: name to `[residues, mass, type]`.
    pub fn modifications(&self) -> &BTreeMap<String, Vec<String>> {
        &self.ptm_name_residues_mass_type
    }

    /// The sequence database.
    pub fn database(&self) -> &str {
        &self.database
    }
    /// Set the sequence database.
    pub fn set_database(&mut self, database: &str) {
        self.database = database.to_owned();
    }
    /// Whether neutral losses are considered for a-, b- and y-ions, e.g.
    /// `0 1 1`.
    pub fn neutral_losses_for_ions(&self) -> &str {
        &self.neutral_losses_for_ions
    }
    /// Set the neutral-loss switches for a-, b- and y-ions.
    pub fn set_neutral_losses_for_ions(&mut self, value: &str) {
        self.neutral_losses_for_ions = value.to_owned();
    }
    /// The weights of the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series.
    pub fn ion_series_weights(&self) -> &str {
        &self.ion_series_weights
    }
    /// Set the ion series weights (space separated).
    pub fn set_ion_series_weights(&mut self, value: &str) {
        self.ion_series_weights = value.to_owned();
    }
    /// Space-separated sequence parts theoretical spectra must contain.
    pub fn partial_sequence(&self) -> &str {
        &self.partial_sequence
    }
    /// Set the partial sequences.
    pub fn set_partial_sequence(&mut self, value: &str) {
        self.partial_sequence = value.to_owned();
    }
    /// Space-separated texts a protein header must contain (or, after a
    /// tilde, lack).
    pub fn sequence_header_filter(&self) -> &str {
        &self.sequence_header_filter
    }
    /// Set the sequence header filter.
    pub fn set_sequence_header_filter(&mut self, value: &str) {
        self.sequence_header_filter = value.to_owned();
    }
    /// The protein mass filter: minimum and maximum mass, or mass and
    /// tolerance in percent.
    pub fn protein_mass_filter(&self) -> &str {
        &self.protein_mass_filter
    }
    /// Set the protein mass filter.
    pub fn set_protein_mass_filter(&mut self, value: &str) {
        self.protein_mass_filter = value.to_owned();
    }
    /// The fragment peak tolerance.
    pub fn peak_mass_tolerance(&self) -> f32 {
        self.peak_mass_tolerance
    }
    /// Set the fragment peak tolerance.
    pub fn set_peak_mass_tolerance(&mut self, value: f32) {
        self.peak_mass_tolerance = value;
    }
    /// The precursor (peptide) mass tolerance.
    pub fn precursor_mass_tolerance(&self) -> f32 {
        self.precursor_mass_tolerance
    }
    /// Set the precursor mass tolerance.
    pub fn set_precursor_mass_tolerance(&mut self, value: f32) {
        self.precursor_mass_tolerance = value;
    }
    /// The minimum distance between two experimental peaks.
    pub fn match_peak_tolerance(&self) -> f32 {
        self.match_peak_tolerance
    }
    /// Set the match peak tolerance.
    pub fn set_match_peak_tolerance(&mut self, value: f32) {
        self.match_peak_tolerance = value;
    }
    /// The cutoff of the ratio of matched to theoretical peaks.
    pub fn ion_cutoff_percentage(&self) -> f32 {
        self.ion_cutoff_percentage
    }
    /// Set the ion cutoff percentage.
    pub fn set_ion_cutoff_percentage(&mut self, value: f32) {
        self.ion_cutoff_percentage = value;
    }
    /// The peptide mass unit: 0 amu, 1 mmu, 2 ppm.
    pub fn peptide_mass_unit(&self) -> usize {
        self.peptide_mass_unit
    }
    /// Set the peptide mass unit.
    pub fn set_peptide_mass_unit(&mut self, value: usize) {
        self.peptide_mass_unit = value;
    }
    /// The number of peptides reported.
    pub fn output_lines(&self) -> usize {
        self.output_lines
    }
    /// Set the number of peptides reported.
    pub fn set_output_lines(&mut self, value: usize) {
        self.output_lines = value;
    }
    /// The maximum number of residues with the same modification in a
    /// peptide (source `getMaxAAPerModPerPeptide`).
    pub fn max_aa_per_mod_per_peptide(&self) -> usize {
        self.max_aa_per_mod_per_peptide
    }
    /// Set the maximum number of residues with the same modification.
    pub fn set_max_aa_per_mod_per_peptide(&mut self, value: usize) {
        self.max_aa_per_mod_per_peptide = value;
    }
    /// The maximum number of modifications in a peptide.
    pub fn max_mods_per_peptide(&self) -> usize {
        self.max_mods_per_peptide
    }
    /// Set the maximum number of modifications in a peptide.
    pub fn set_max_mods_per_peptide(&mut self, value: usize) {
        self.max_mods_per_peptide = value;
    }
    /// The nucleotide reading frame: 0 a protein database, 1–3 forward and
    /// 4–6 reverse translations starting at the first, second or third
    /// base, 7 all forward, 8 all reverse, 9 all six.
    pub fn nucleotide_reading_frame(&self) -> usize {
        self.nucleotide_reading_frame
    }
    /// Set the nucleotide reading frame.
    pub fn set_nucleotide_reading_frame(&mut self, value: usize) {
        self.nucleotide_reading_frame = value;
    }
    /// The maximum number of internal cleavage sites.
    pub fn max_internal_cleavage_sites(&self) -> usize {
        self.max_internal_cleavage_sites
    }
    /// Set the maximum number of internal cleavage sites.
    pub fn set_max_internal_cleavage_sites(&mut self, value: usize) {
        self.max_internal_cleavage_sites = value;
    }
    /// The number of most abundant peaks matched against theoretical ones.
    pub fn match_peak_count(&self) -> usize {
        self.match_peak_count
    }
    /// Set the match peak count.
    pub fn set_match_peak_count(&mut self, value: usize) {
        self.match_peak_count = value;
    }
    /// How many of those peaks may fail to match.
    pub fn match_peak_allowed_error(&self) -> usize {
        self.match_peak_allowed_error
    }
    /// Set the match peak allowed error.
    pub fn set_match_peak_allowed_error(&mut self, value: usize) {
        self.match_peak_allowed_error = value;
    }
    /// Whether fragment ions are shown.
    pub fn show_fragment_ions(&self) -> bool {
        self.show_fragment_ions
    }
    /// Set whether fragment ions are shown.
    pub fn set_show_fragment_ions(&mut self, value: bool) {
        self.show_fragment_ions = value;
    }
    /// Whether every protein containing a found peptide is reported.
    pub fn print_duplicate_references(&self) -> bool {
        self.print_duplicate_references
    }
    /// Set whether every protein containing a found peptide is reported.
    pub fn set_print_duplicate_references(&mut self, value: bool) {
        self.print_duplicate_references = value;
    }
    /// Whether peaks near (15 amu) the precursor are removed.
    pub fn remove_precursor_near_peaks(&self) -> bool {
        self.remove_precursor_near_peaks
    }
    /// Set whether peaks near the precursor are removed.
    pub fn set_remove_precursor_near_peaks(&mut self, value: bool) {
        self.remove_precursor_near_peaks = value;
    }
    /// The parent mass type. The header documents 0 as monoisotopic, while
    /// the written comment reads `0=average masses, 1=monoisotopic masses`.
    pub fn mass_type_parent(&self) -> bool {
        self.mass_type_parent
    }
    /// Set the parent mass type.
    pub fn set_mass_type_parent(&mut self, value: bool) {
        self.mass_type_parent = value;
    }
    /// The fragment mass type (see [`SequestInfile::mass_type_parent`]).
    pub fn mass_type_fragment(&self) -> bool {
        self.mass_type_fragment
    }
    /// Set the fragment mass type.
    pub fn set_mass_type_fragment(&mut self, value: bool) {
        self.mass_type_fragment = value;
    }
    /// Whether normalised XCorr values are reported.
    pub fn normalize_xcorr(&self) -> bool {
        self.normalize_xcorr
    }
    /// Set whether normalised XCorr values are reported.
    pub fn set_normalize_xcorr(&mut self, value: bool) {
        self.normalize_xcorr = value;
    }
    /// Whether residues are in upper case.
    pub fn residues_in_upper_case(&self) -> bool {
        self.residues_in_upper_case
    }
    /// Set whether residues are in upper case.
    pub fn set_residues_in_upper_case(&mut self, value: bool) {
        self.residues_in_upper_case = value;
    }
}

/// The static residue modification lines of the source's `store`, in its
/// order: residue, parameter name, comment.
const RESIDUE_LINES: [(char, &str, &str); 24] = [
    (
        'G',
        "add_G_Glycine",
        "added to G - avg.  57.0519, mono.  57.02146",
    ),
    (
        'A',
        "add_A_Alanine",
        "added to A - avg.  71.0788, mono.  71.03711",
    ),
    (
        'S',
        "add_S_Serine",
        "added to S - avg.  87.0782, mono.  87.03203",
    ),
    (
        'P',
        "add_P_Proline",
        "added to P - avg.  97.1167, mono.  97.05276",
    ),
    (
        'V',
        "add_V_Valine",
        "added to V - avg.  99.1326, mono.  99.06841",
    ),
    (
        'T',
        "add_T_Threonine",
        "added to T - avg. 101.1051, mono. 101.04768",
    ),
    (
        'C',
        "add_C_Cysteine",
        "added to C - avg. 103.1388, mono. 103.00919",
    ),
    (
        'L',
        "add_L_Leucine",
        "added to L - avg. 113.1594, mono. 113.08406",
    ),
    (
        'I',
        "add_I_Isoleucine",
        "added to I - avg. 113.1594, mono. 113.08406",
    ),
    (
        'X',
        "add_X_LorI",
        "added to X - avg. 113.1594, mono. 113.08406",
    ),
    (
        'N',
        "add_N_Asparagine",
        "added to N - avg. 114.1038, mono. 114.04293",
    ),
    (
        'O',
        "add_O_Ornithine",
        "added to O - avg. 114.1472, mono  114.07931",
    ),
    (
        'B',
        "add_B_avg_NandD",
        "added to B - avg. 114.5962, mono. 114.53494",
    ),
    (
        'D',
        "add_D_Aspartic_Acid",
        "added to D - avg. 115.0886, mono. 115.02694",
    ),
    (
        'Q',
        "add_Q_Glutamine",
        "added to Q - avg. 128.1307, mono. 128.05858",
    ),
    (
        'K',
        "add_K_Lysine",
        "added to K - avg. 128.1741, mono. 128.09496",
    ),
    (
        'Z',
        "add_Z_avg_QandE",
        "added to Z - avg. 128.6231, mono. 128.55059",
    ),
    (
        'E',
        "add_E_Glutamic_Acid",
        "added to E - avg. 129.1155, mono. 129.04259",
    ),
    (
        'M',
        "add_M_Methionine",
        "added to M - avg. 131.1926, mono. 131.04049",
    ),
    (
        'H',
        "add_H_Histidine",
        "added to H - avg. 137.1411, mono. 137.05891",
    ),
    (
        'F',
        "add_F_Phenylalanine",
        "added to F - avg. 147.1766, mono. 147.06841",
    ),
    (
        'R',
        "add_R_Arginine",
        "added to R - avg. 156.1875, mono. 156.10111",
    ),
    (
        'Y',
        "add_Y_Tyrosine",
        "added to Y - avg. 163.1760, mono. 163.06333",
    ),
    (
        'W',
        "add_W_Tryptophan",
        "added to W - avg. 186.2132, mono. 186.07931",
    ),
];

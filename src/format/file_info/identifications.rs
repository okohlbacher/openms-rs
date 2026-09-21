// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo summary of idXML and mzIdentML identifications
//! (`FORMAT/FileInfo.cpp:1312-1470`, `:1990-1996`, `:2105-2107`, `:2373-2376`).
//!
//! The identification branch of the report, shared by idXML and mzIdentML. The
//! runs are loaded through [`FileHandler::load_identifications`](crate::format::FileHandler::load_identifications) for idXML and
//! through [`crate::format::mzidentml::load`] for mzIdentML, and the branch
//! then writes:
//!
//! 1. three TSV lines — database, database version and taxonomy — taken from
//!    the *first* protein identification run. The source writes these before
//!    anything else, and before it has established that a first run exists;
//!    see *What this port refuses*;
//! 2. the search engines, one line each, deduplicated and ordered as the
//!    source's `set<pair<string, string>>` orders them;
//! 3. the run, protein-hit and non-redundant protein-hit counts;
//! 4. the matched-spectrum, peptide-sequence, PSM-per-spectrum, peptide-hit,
//!    modified-top-hit and non-redundant peptide-hit counts;
//! 5. one line of modification counts, when there are any, which does *not*
//!    end in a newline.
//!
//! Two numbers on those lines are deliberately coarse in the source and are
//! reproduced as they are:
//!
//! - `PSMs / spectrum` is an integer division of a `Size` by an `int`, so it
//!   truncates. The structured [`IdentInfo::psms_per_spectrum`](crate::format::file_info::model::IdentInfo::psms_per_spectrum) carries the
//!   real ratio, which is what the source's own `Result` records;
//! - the average peptide length is `Math::round` of the mean of the hit
//!   lengths, streamed at the report's precision rather than through
//!   `StringUtils::toStr`, while the modified-top-hit percentage goes the
//!   other way and is rendered by `StringUtils::toStr`.
//!
//! Modification counting distinguishes the two identities the source uses: a
//! terminal modification is counted under `getId()`, which is empty for a
//! user-defined (mass-only) modification, and a residue modification under
//! `getFullId()`.
//!
//! # The mzIdentML fall-through
//!
//! `-m`, `-p` and `-s` each have an `IDXML` arm but none for `MZIDENTML`, so an
//! mzIdentML input falls into their trailing `else //peaks` arm
//! (`FORMAT/FileInfo.cpp:2005`, `:2115`, `:2384`) and is reported off the
//! `MSExperiment` that this branch never filled. The result is the peak-file
//! metadata layout with every field empty, the peak-file data-processing arm
//! finding an empty experiment, and a peak-file `Intensities:` statistics block
//! over no values. That is what the reference build prints, so it is what this
//! port prints.
//!
//! # What this port refuses
//!
//! Two out-of-bounds `std::vector` accesses, both of which the reference build
//! answers with a segmentation fault. Lead decision D1 refuses exactly there,
//! before any of the report is written:
//!
//! - `FORMAT/FileInfo.cpp:1336-1341` reads `id_data.proteins[0]` unconditionally,
//!   while the structured block at `:1451` guards the very same access with
//!   `if (!id_data.proteins.empty())`. A file with no identification run at all
//!   reaches the unguarded read. This crate refuses such a file, but one layer
//!   earlier: both readers reject a run-less document themselves — idXML with
//!   `idXML needs at least one IdentificationRun` (`idxml.rs:310`) and
//!   mzIdentML with `mzIdentML has no SpectrumIdentification element`
//!   (`mzidentml.rs:1596-1598`) — so the guard below is a defence with no
//!   reachable input rather than the refusal that is measured;
//! - `FORMAT/FileInfo.cpp:1354` reads `temp_hits[0]`, the first element of the
//!   `getHits()` reference taken at `:1352`, behind the `:1347` guard
//!   `if (!id_data.peptides[i].empty())`. But `PeptideIdentification::empty()`
//!   (`PeptideIdentification.cpp:210-217`) tests for a default-constructed
//!   object rather than for an empty hit list: a score type, an identifier, a
//!   non-zero significance threshold or `higher_score_better == false` each
//!   make it false on their own. A hit-less identification carrying any of
//!   those reaches the unguarded read.
//!
//! See `docs/FILE_INFO_A7_SUPPORT.md` for the evidence and the native
//! differences.

#![cfg(feature = "idxml")]

use super::model::{FileInfoResult, IdentInfo, Options};
use super::report::{
    ReportStream, statistics_buffer, summarize, write_meta_title, write_processing,
    write_processing_title, write_statistics_title, write_summary_text,
};
use super::text_format::{WRITTEN_DIGITS_F32, to_str};
use crate::chemistry::{AASequence, SequenceModification};
use crate::concept::math_functions::round;
use crate::format::FileHandler;
use crate::format::file_types::FileType;
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::MSExperiment;
use crate::math::statistic_functions::mean;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// The three fields of an identification document the branch consumes.
struct IdData {
    proteins: Vec<ProteinIdentification>,
    peptides: Vec<PeptideIdentification>,
    /// `IdXMLFile::load`'s third output. `FileHandler::loadIdentifications`
    /// does not fill it, so an mzIdentML input leaves it empty, which is what
    /// the `-m` arm would print if mzIdentML reached it.
    identifier: String,
}

/// Load the identifications and write the branch, `-m`, `-p` and `-s`.
pub(crate) fn report(
    path: &Path,
    in_type: FileType,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let data = load(path, in_type)?;
    // FORMAT/FileInfo.cpp:1336-1341 reads proteins[0] before any emptiness test. No
    // input reaches this arm of the guard: `load` above refuses a run-less
    // idXML and a run-less mzIdentML in the readers themselves, so the tool
    // exits 3 on such a file, never 6. It stays as the branch's own defence.
    let first_run = data.proteins.first().ok_or_else(|| {
        Error::InvalidValue(
            "FileInfo identification branch: the file holds no protein identification run, and \
             the source reads id_data.proteins[0] unguarded (FORMAT/FileInfo.cpp:1336-1341) while its \
             own structured block guards the same access (FORMAT/FileInfo.cpp:1451)"
                .into(),
        )
    })?;
    let parameters = &first_run.search_parameters;
    os_tsv
        .text("general: database\t")
        .text(&parameters.database)
        .text("\ngeneral: database version\t")
        .text(&parameters.database_version)
        .text("\ngeneral: taxonomy\t")
        .text(&parameters.taxonomy)
        .text("\n");

    let mut spectrum_count = 0_u64;
    let mut peptide_hit_count = 0_u64;
    let mut average_peptide_hits = 0_u64;
    let mut modified_peptide_count = 0_u64;
    let mut peptides: BTreeSet<String> = BTreeSet::new();
    let mut peptides_ignore_mods: BTreeSet<String> = BTreeSet::new();
    let mut modification_counts: BTreeMap<String, u64> = BTreeMap::new();
    let mut peptide_length: Vec<f64> = Vec::new();

    for (index, identification) in data.peptides.iter().enumerate() {
        if source_is_empty(identification) {
            continue;
        }
        spectrum_count += 1;
        let hits = count(identification.hits.len())?;
        average_peptide_hits = average_peptide_hits
            .checked_add(hits)
            .ok_or_else(|| overflow("FileInfo peptide hit count overflows 64 bits"))?;
        peptide_hit_count = peptide_hit_count
            .checked_add(hits)
            .ok_or_else(|| overflow("FileInfo peptide hit count overflows 64 bits"))?;
        // FORMAT/FileInfo.cpp:1354 reads temp_hits[0] behind a guard that does not
        // test the hit list.
        let top = identification.hits.first().ok_or_else(|| {
            Error::InvalidValue(format!(
                "FileInfo identification branch: peptide identification #{index} carries no hit \
                 while PeptideIdentification::empty() is false, and the source reads \
                 getHits()[0] unguarded (FORMAT/FileInfo.cpp:1347-1354)"
            ))
        })?;
        if top.sequence.is_modified() {
            modified_peptide_count += 1;
            count_modifications(&top.sequence, &mut modification_counts)?;
        }
        for hit in &identification.hits {
            peptides.insert(hit.sequence.to_string());
            peptides_ignore_mods.insert(hit.sequence.as_str().to_owned());
            peptide_length.push(residue_count_as_u16(hit.sequence.len()));
        }
    }

    let mut runs_count = 0_u64;
    let mut protein_hit_count = 0_u64;
    let mut proteins_seen: BTreeSet<&str> = BTreeSet::new();
    let mut search_engines: BTreeSet<(&str, &str)> = BTreeSet::new();
    for run in &data.proteins {
        runs_count += 1;
        protein_hit_count = protein_hit_count
            .checked_add(count(run.hits.len())?)
            .ok_or_else(|| overflow("FileInfo protein hit count overflows 64 bits"))?;
        for hit in &run.hits {
            proteins_seen.insert(hit.accession.as_str());
        }
        search_engines.insert((
            run.search_engine.as_str(),
            run.search_engine_version.as_str(),
        ));
    }

    // FORMAT/FileInfo.cpp:1396-1399: a single zero keeps Math::mean off an empty range.
    if peptide_length.is_empty() {
        peptide_length.push(0.0);
    }
    let average_length = round(mean(&peptide_length)?);

    os.text("Search Engine(s):\n");
    for (engine, version) in &search_engines {
        os.text("  ")
            .text(engine)
            .text(" (version: ")
            .text(version)
            .text(")\n");
    }
    os.text("Number of:\n");
    os.text("  runs:                       ")
        .value(runs_count)
        .text("\n");
    os.text("  protein hits:               ")
        .value(protein_hit_count)
        .text("\n");
    os.text("  non-redundant protein hits: ")
        .value(proteins_seen.len())
        .text("\n");
    os.text("  (only hits that differ in the accession)\n");
    os.text("\n");
    os.text("  matched spectra:    ")
        .value(spectrum_count)
        .text("\n");
    os.text("  peptide sequences:  ")
        .value(peptides_ignore_mods.len())
        .text("\n");
    os.text("  PSMs / spectrum (ignoring unidentified spectra):    ")
        .value(average_peptide_hits / divisor(spectrum_count))
        .text("\n");
    os.text("  peptide hits:               ")
        .value(peptide_hit_count)
        .text(" (avg. length: ")
        .double(average_length)
        .text(")\n");
    os.text("  modified top-hits:          ")
        .value(modified_peptide_count)
        .text("/")
        .value(spectrum_count)
        .text(&modified_percentage(modified_peptide_count, spectrum_count))
        .text("\n");
    os.text("  non-redundant peptide hits: ")
        .value(peptides.len())
        .text("\n");
    os.text("  (only hits that differ in sequence and/or modifications)\n");
    for (index, (name, count)) in modification_counts.iter().enumerate() {
        if index == 0 {
            os.text("  Modification count (top-hits only): ");
        } else {
            os.text(", ");
        }
        os.text(name).text(" ").value(count);
    }

    for (engine, version) in &search_engines {
        os_tsv
            .text("general: search engine\t")
            .text(engine)
            .text("\t(version: ")
            .text(version)
            .text(")\n");
    }
    os_tsv
        .text("general: num. of runs\t")
        .value(runs_count)
        .text("\n");
    os_tsv
        .text("general: num. of protein hits\t")
        .value(protein_hit_count)
        .text("\n");
    os_tsv
        .text("general: num. of non-redundant protein hits (only hits that differ in the accession)\t")
        .value(proteins_seen.len())
        .text("\n");
    os_tsv
        .text("general: num. of matched spectra\t")
        .value(spectrum_count)
        .text("\n");
    os_tsv
        .text("general: num. of peptide hits\t")
        .value(peptide_hit_count)
        .text("\n");
    os_tsv
        .text("general: num. of modified top-hits\t")
        .value(modified_peptide_count)
        .text("\n");
    // The source's trailing space before the tab is part of the label.
    os_tsv
        .text(
            "general: num. of non-redundant peptide hits (only hits that differ in sequence \
               and/or modifications): \t",
        )
        .value(peptides.len())
        .text("\n");

    result.ident = Some(IdentInfo {
        db_name: parameters.database.clone(),
        db_version: parameters.database_version.clone(),
        taxonomy: parameters.taxonomy.clone(),
        search_engines: search_engines
            .iter()
            .map(|(engine, version)| format!("{engine} (version: {version})"))
            .collect(),
        num_runs: runs_count,
        protein_hits: protein_hit_count,
        non_redundant_protein_hits: count(proteins_seen.len())?,
        matched_spectra: spectrum_count,
        peptide_hits: peptide_hit_count,
        // The real ratio, which the CLI text deliberately does not print.
        psms_per_spectrum: as_double(average_peptide_hits) / as_double(divisor(spectrum_count)),
        avg_peptide_length: mean(&peptide_length)?,
        non_redundant_peptides: count(peptides.len())?,
        modified_tophits: modified_peptide_count,
        modification_counts,
    });

    write_trailing_sections(in_type, &data, options, os, os_tsv, result)
}

/// The `-m`, `-p` and `-s` sections. Only idXML has an arm of its own in each;
/// mzIdentML falls through to the peak-file arm, which runs off an experiment
/// this branch never loaded and therefore reports an empty one.
fn write_trailing_sections(
    in_type: FileType,
    data: &IdData,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let is_idxml = in_type == FileType::IdXml;
    if options.meta {
        if is_idxml {
            write_meta_title(os);
            os.text("Document ID: ").text(&data.identifier).text("\n\n");
            // FORMAT/FileInfo.cpp:1994-1995 streams the value with no trailing newline.
            os_tsv.text("meta: document ID\t").text(&data.identifier);
        } else {
            // FORMAT/FileInfo.cpp:2005-2081: the peak-file arm, over the MSExperiment
            // this branch never loaded.
            super::peaks::write_meta(&MSExperiment::default(), os, os_tsv);
        }
    }
    if options.processing {
        write_processing_title(os);
        // Both arms end with no data processing: the idXML arm is empty and
        // the peak-file arm finds `exp.empty()`.
        write_processing(os, os_tsv, &[], result);
    }
    if options.statistics {
        write_statistics_title(os);
        if !is_idxml {
            // FORMAT/FileInfo.cpp:2384-2434 over an empty experiment: no MS-level-1
            // peak contributes an intensity and `meta_names` is empty, so the
            // arm writes one all-zero block at writtenDigits<float>() and
            // nothing to the TSV report.
            let mut values = statistics_buffer(0)?;
            let stats = summarize(&mut values)?;
            os.set_precision(WRITTEN_DIGITS_F32);
            os.text("Intensities:\n");
            write_summary_text(os, &stats);
            os.text("\n");
        }
    }
    Ok(())
}

/// `IdXMLFile::load` for idXML, `FileHandler::loadIdentifications` restricted
/// to `{MZIDENTML}` for mzIdentML.
fn load(path: &Path, in_type: FileType) -> Result<IdData> {
    if in_type == FileType::MzIdentMl {
        let detected = FileHandler::get_type(path)?;
        if detected != FileType::MzIdentMl {
            return Err(Error::InvalidValue(format!(
                "{} is not an allowed input format",
                detected.name()
            )));
        }
        let document = crate::format::mzidentml::load(path)?;
        return Ok(IdData {
            proteins: document.protein_identifications,
            peptides: document.peptide_identifications,
            identifier: String::new(),
        });
    }
    let document = FileHandler::load_identifications(path, &[FileType::IdXml])?;
    Ok(IdData {
        proteins: document.protein_identifications,
        peptides: document.peptide_identifications,
        identifier: document.document_id,
    })
}

/// `PeptideIdentification::empty()` (`PeptideIdentification.cpp:210-217`): a
/// default-constructed object, not an empty hit list.
fn source_is_empty(identification: &PeptideIdentification) -> bool {
    identification.identifier.is_empty()
        && identification.hits.is_empty()
        && identification.significance_threshold == 0.0
        && identification.score_type.is_empty()
        && identification.higher_score_better
}

/// `FORMAT/FileInfo.cpp:1353-1372`: the C-terminal modification, then the N-terminal
/// one, then every modified residue in order.
///
/// A terminal modification is counted under `getId()`, which
/// `AASequence::getNTerminalModificationName` and its C-terminal twin return;
/// that is the empty string for a user-defined mass-only modification, which
/// carries a full identifier but no identifier
/// (`ResidueModification.cpp:593`, `:631`, `:671`). A residue modification is
/// counted under `getFullId()` instead.
fn count_modifications(
    sequence: &AASequence,
    modification_counts: &mut BTreeMap<String, u64>,
) -> Result<()> {
    if let Some(modification) = sequence.c_terminal_modification() {
        *modification_counts
            .entry(terminal_name(modification).to_owned())
            .or_insert(0) += 1;
    }
    if let Some(modification) = sequence.n_terminal_modification() {
        *modification_counts
            .entry(terminal_name(modification).to_owned())
            .or_insert(0) += 1;
    }
    for index in 0..sequence.len() {
        if let Some(modification) = sequence.residue_modification(index)? {
            *modification_counts
                .entry(modification.full_id().to_owned())
                .or_insert(0) += 1;
        }
    }
    Ok(())
}

/// `ResidueModification::getId()`: the registry identifier, and the empty
/// string for a user-defined mass-only modification.
fn terminal_name(modification: &SequenceModification) -> &str {
    match modification.known() {
        Some(known) => known.name(),
        None => "",
    }
}

/// `std::max(1, (Int)spectrum_count)` as the divisor of the printed ratio.
fn divisor(spectrum_count: u64) -> u64 {
    spectrum_count.max(1)
}

/// `(spectrum_count > 0 ? std::string(" (") + Math::round(m * 1000.0 / s) / 10 + "%)" : "")`.
///
/// The percentage goes through `StringUtils::toStr(double)`, not through the
/// report's stream precision, because the source builds a `std::string` here.
fn modified_percentage(modified: u64, spectrum_count: u64) -> String {
    if spectrum_count == 0 {
        return String::new();
    }
    let percent = round(as_double(modified) * 1000.0 / as_double(spectrum_count)) / 10.0;
    format!(" ({}%)", to_str(percent))
}

/// `(uint16_t)temp_hits[j].getSequence().size()`, then promoted to `double` by
/// `Math::mean`. The cast keeps the source's 16-bit truncation for a peptide
/// of more than 65535 residues.
#[expect(
    clippy::cast_possible_truncation,
    reason = "the source casts the residue count to uint16_t in exactly this place"
)]
fn residue_count_as_u16(residues: usize) -> f64 {
    f64::from(residues as u16)
}

#[expect(
    clippy::cast_precision_loss,
    reason = "the source converts its Size counts to double in exactly this place"
)]
fn as_double(value: u64) -> f64 {
    value as f64
}

fn count(value: usize) -> Result<u64> {
    u64::try_from(value).map_err(|_| overflow("FileInfo identification count overflows 64 bits"))
}

fn overflow(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

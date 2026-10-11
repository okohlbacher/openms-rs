// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Export of isobaric consensus maps into the IBSpectra format used by the R
//! package isobar to load quantification results, ported from
//! `FORMAT/IBSpectraFile.h`. See `docs/IBSPECTRA_FILE_SUPPORT.md`.

use crate::chemistry::{AASequence, PeptideFragmentType};
use crate::format::file_info::text_format::to_str;
use crate::kernel::ConsensusMap;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// Intensities at or below this value are written as 0.0 "to avoid numerical
/// problems when quantifying" (the source's local `intensity_threshold`).
const INTENSITY_THRESHOLD: f64 = 0.00001;

/// Reporter-ion channel centres of the isobaric methods the writer can guess,
/// transcribed from the source's `ItraqFourPlexQuantitationMethod.cpp`,
/// `TMTSixPlexQuantitationMethod.cpp` and `ItraqEightPlexQuantitationMethod.cpp`
/// constructors (core `bc9cc12`, lines 23-26, 23-28 and 25-32). Only the
/// centres are needed here; the `ANALYSIS/QUANTITATION` method classes are not
/// part of this port.
const ITRAQ_4PLEX: [f64; 4] = [114.1112, 115.1082, 116.1116, 117.1149];
const TMT_6PLEX: [f64; 6] = [
    126.127725, 127.124760, 128.134433, 129.131468, 130.141141, 131.138176,
];
const ITRAQ_8PLEX: [f64; 8] = [
    113.1078, 114.1112, 115.1082, 116.1116, 117.1149, 118.1120, 119.1153, 121.1220,
];

/// Implements the export of consensus maps into the IBSpectra format used by
/// isobar to load quantification results (source `IBSpectraFile`).
///
/// The source class has no members; its copy constructor and assignment
/// operator are [`Clone`] and [`Copy`] here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IBSpectraFile;

/// Holds all identification information of one output line (the source's
/// file-local `IdCSV`), with the source's defaults for an unidentified feature.
#[derive(Clone, Debug)]
struct IdCsv {
    accession: String,
    peptide: String,
    modif: String,
    charge: i32,
    theo_mass: f64,
    exp_mass: f64,
    parent_intens: f64,
    retention_time: f64,
    spectrum: String,
    search_engine: &'static str,
}

impl Default for IdCsv {
    fn default() -> Self {
        Self {
            accession: "UNIDENTIFIED_PROTEIN".into(),
            peptide: "UNIDENTIFIED_PEPTIDE".into(),
            modif: String::new(),
            charge: 0,
            theo_mass: -1.0,
            exp_mass: -1.0,
            parent_intens: -1.0,
            retention_time: -1.0,
            spectrum: String::new(),
            search_engine: "open-ms-generic",
        }
    }
}

impl IdCsv {
    fn push_fields(&self, line: &mut Vec<String>) {
        line.push(self.accession.clone());
        line.push(self.peptide.clone());
        line.push(self.modif.clone());
        line.push(self.charge.to_string());
        line.push(to_str(self.theo_mass));
        line.push(to_str(self.exp_mass));
        line.push(to_str(self.parent_intens));
        line.push(to_str(self.retention_time));
        line.push(self.spectrum.clone());
        line.push(self.search_engine.to_owned());
    }
}

impl IBSpectraFile {
    /// A writer (source default constructor).
    pub fn new() -> Self {
        Self
    }

    /// Writes the contents of `cm` into the file named by `filename`.
    ///
    /// The whole table is built first and the file is written only when that
    /// succeeded, as the source builds a `TextFile` and stores it at the end.
    /// The layout is that of [`Self::write`].
    ///
    /// # Errors
    ///
    /// As [`Self::write`], and [`Error::Io`] when the file cannot be created
    /// or written (the source's `Exception::UnableToCreateFile`).
    pub fn store(&self, filename: impl AsRef<Path>, cm: &ConsensusMap) -> Result<()> {
        let text = self.table(cm)?;
        std::fs::write(filename.as_ref(), text)?;
        Ok(())
    }

    /// Writes the IBSpectra table of `cm` to `writer`.
    ///
    /// The isobaric method is guessed from the map: the experiment type must be
    /// `labeled_MS2` or `itraq`, and the number of column headers chooses
    /// iTRAQ 4-plex (4), TMT 6-plex (6) or iTRAQ 8-plex (8).
    ///
    /// The tab-separated header is `accession`, `peptide`, `modif`, `charge`,
    /// `theo.mass`, `exp.mass`, `parent.intens`, `retention.time`, `spectrum`,
    /// `search.engine`, then `X<c>_mass` and `X<c>_ions` for every channel
    /// centre `c` truncated to an integer.
    ///
    /// Each consensus feature yields lines as follows. Without a peptide
    /// identification, or when the map has no protein identification, one
    /// line for `UNIDENTIFIED_PROTEIN` / `UNIDENTIFIED_PEPTIDE` is written,
    /// because the quantification still helps normalisation. Otherwise the
    /// first hit of the first peptide identification gives one line per protein
    /// accession of its evidences (in sorted order) that the first protein
    /// identification of the map contains; an accession it lacks is skipped
    /// (the source reports it on standard error; this port does not), and a hit
    /// without accessions writes nothing. Non-unique peptides are written, as
    /// the source's fixed `allow_non_unique = true` asks. A feature of
    /// intensity 0 writes nothing.
    ///
    /// A line holds the accession, the unmodified peptide, the modification
    /// string (the N-terminal modification name, then `:` and the modification
    /// name of each residue, then `:` and the C-terminal name if there is one),
    /// the hit charge, the theoretical mass (full monoisotopic mass plus
    /// `charge` protons), the consensus m/z, intensity and RT, the consensus
    /// unique ID and `open-ms-generic`, then the channel centres and the
    /// channel intensities. A channel's intensity is that of the handle whose
    /// m/z truncates to the same integer as the centre (the last such handle in
    /// handle order), 0.0 when it is at most 1e-5 or no handle matches. Numbers
    /// use the source's `StringUtils::toStr`.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] with the source's message when the map does
    ///   not hold the result of an isobaric quantification experiment ("Given
    ///   ConsensusMap does not hold any isobaric quantification data.") or its
    ///   column count matches no method ("Could not guess isobaric
    ///   quantification data from ConsensusMap due to non-matching number of
    ///   input maps."); the source throws `Exception::InvalidParameter`.
    /// - [`Error::InvalidValue`] when the first peptide identification of an
    ///   identified feature has no hit, which the source indexes without a
    ///   check; and the errors of the theoretical mass for a sequence with a
    ///   residue of unknown mass.
    /// - [`Error::Io`] when writing fails.
    pub fn write(&self, mut writer: impl Write, cm: &ConsensusMap) -> Result<()> {
        let text = self.table(cm)?;
        writer.write_all(text.as_bytes())?;
        writer.flush()?;
        Ok(())
    }

    /// Guesses the isobaric method from the map (source private
    /// `guessExperimentType_`), returning its channel centres.
    fn guess_experiment_type(cm: &ConsensusMap) -> Result<&'static [f64]> {
        if cm.experiment_type != "labeled_MS2" && cm.experiment_type != "itraq" {
            return Err(Error::InvalidValue(
                "Given ConsensusMap does not hold any isobaric quantification data.".into(),
            ));
        }
        // The map count is taken as the approximation.
        match cm.column_headers.len() {
            4 => Ok(&ITRAQ_4PLEX),
            6 => Ok(&TMT_6PLEX),
            8 => Ok(&ITRAQ_8PLEX),
            _ => Err(Error::InvalidValue(
                "Could not guess isobaric quantification data from ConsensusMap due to non-matching number of input maps."
                    .into(),
            )),
        }
    }

    /// The header of the file for the given channels (source private
    /// `constructHeader_`).
    fn construct_header(channels: &[f64]) -> Vec<String> {
        let mut header: Vec<String> = [
            "accession",
            "peptide",
            "modif",
            "charge",
            "theo.mass",
            "exp.mass",
            "parent.intens",
            "retention.time",
            "spectrum",
            "search.engine",
        ]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
        for center in channels {
            header.push(format!("X{}_mass", *center as i32));
        }
        for center in channels {
            header.push(format!("X{}_ions", *center as i32));
        }
        header
    }

    /// The modification string of `sequence` (source private
    /// `getModifString_`).
    fn modif_string(sequence: &AASequence) -> Result<String> {
        let mut modif = sequence
            .n_terminal_modification()
            .map_or("", |m| m.name())
            .to_owned();
        for index in 0..sequence.len() {
            modif.push(':');
            if let Some(modification) = sequence.residue_modification(index)? {
                modif.push_str(modification.name());
            }
        }
        if let Some(c_term) = sequence.c_terminal_modification() {
            if !c_term.name().is_empty() {
                modif.push(':');
                modif.push_str(c_term.name());
            }
        }
        Ok(modif)
    }

    fn table(&self, cm: &ConsensusMap) -> Result<String> {
        let channels = Self::guess_experiment_type(cm)?;
        let proteins = cm.protein_identifications.first();

        let mut text = Self::construct_header(channels).join("\t");
        text.push('\n');

        for feature in &cm.features {
            let base = &feature.base;
            let mut entries = Vec::new();
            match (base.peptide_identifications.first(), proteins) {
                (Some(identification), Some(proteins)) => {
                    let hit = identification.hits.first().ok_or_else(|| {
                        Error::InvalidValue(
                            "IBSpectra export: a consensus peptide identification has no hit"
                                .into(),
                        )
                    })?;
                    for accession in hit.protein_accessions() {
                        let Some(protein) = proteins.find_hit(accession) else {
                            // "Protein referenced in peptide not found..."
                            continue;
                        };
                        entries.push(IdCsv {
                            accession: protein.accession.clone(),
                            peptide: hit.sequence.as_str().to_owned(),
                            modif: Self::modif_string(&hit.sequence)?,
                            charge: hit.charge,
                            theo_mass: hit
                                .sequence
                                .mono_mass_for(PeptideFragmentType::Full, hit.charge)?,
                            ..IdCsv::default()
                        });
                    }
                }
                // Unidentified hits are stored anyway, because the quantification
                // is still helpful for normalization.
                _ => entries.push(IdCsv::default()),
            }

            // Skip features with 0 intensity.
            if base.intensity == 0.0 {
                continue;
            }

            let mut intensities: BTreeMap<i32, f64> = BTreeMap::new();
            for handle in feature.handles() {
                let intensity = f64::from(handle.intensity);
                intensities.insert(
                    handle.mz as i32,
                    if intensity > INTENSITY_THRESHOLD {
                        intensity
                    } else {
                        0.0
                    },
                );
            }
            for mut entry in entries {
                entry.parent_intens = f64::from(base.intensity);
                entry.retention_time = base.rt;
                entry.spectrum = base.unique_id.to_string();
                entry.exp_mass = base.mz;
                let mut line = Vec::with_capacity(10 + 2 * channels.len());
                entry.push_fields(&mut line);
                for center in channels {
                    line.push(to_str(*center));
                }
                for center in channels {
                    let value = intensities.get(&(*center as i32)).copied().unwrap_or(0.0);
                    line.push(to_str(value));
                }
                text.push_str(&line.join("\t"));
                text.push('\n');
            }
        }
        Ok(text)
    }
}

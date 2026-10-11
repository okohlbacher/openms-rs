// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Identification data as a Parquet directory bundle (`.idparquet`).
//!
//! Ports `FORMAT/PSMArrowIO.h` and `PSMArrowIO.cpp`. An idparquet bundle is a
//! directory holding four Parquet files:
//!
//! * `psms.parquet`: the peptide identifications, lossless, in the internal
//!   PSM layout of [`QpxFile::export_to_arrow`](crate::format::qpx_file::QpxFile::export_to_arrow);
//! * `proteins.parquet`: the protein hits;
//! * `protein_groups.parquet`: the protein and indistinguishable groups;
//! * `search_params.parquet`: the run-level parameters, with the
//!   modification definitions the identifications use.
//!
//! The three protein-side files are those of
//! [`ProteinIdentificationArrowIO`](crate::format::protein_identification_arrow_io::ProteinIdentificationArrowIO).
//! All four are required on read. See `docs/ARROW_PARQUET_SUPPORT.md` for the
//! member mapping and the evidence.

use crate::chemistry::ModificationsDB;
use crate::format::arrow_io_helpers::{ParquetWriteConfig, write_table_to_parquet};
use crate::format::modification_definitions;
use crate::format::parquet_file::ParquetFile;
use crate::format::protein_identification_arrow_io::ProteinIdentificationArrowIO;
use crate::format::qpx_file::{QpxFile, QpxOutput};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::{Error, Result};
use std::path::Path;

/// The PSM table's file name in a bundle (source `kPSMs`).
pub const PSMS: &str = "psms.parquet";
/// The protein table's file name in a bundle (source `kProteins`).
pub const PROTEINS: &str = "proteins.parquet";
/// The protein-group table's file name in a bundle (source `kProteinGroups`).
pub const PROTEIN_GROUPS: &str = "protein_groups.parquet";
/// The search-parameter table's file name in a bundle (source `kSearchParams`).
pub const SEARCH_PARAMS: &str = "search_params.parquet";

/// Read and write OpenMS identification data as a Parquet bundle
/// (`.idparquet`); source `PSMArrowIO`, a class of static members.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PSMArrowIO;

impl PSMArrowIO {
    /// Export protein and peptide identifications to an idparquet directory
    /// bundle, resolving the search parameters' modification names in the
    /// global [`ModificationsDB`].
    ///
    /// [`export_to_parquet_with_registry`](Self::export_to_parquet_with_registry)
    /// with [`ModificationsDB::global`].
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn export_to_parquet(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        dir: impl AsRef<Path>,
        export_all_psms: bool,
        config: &ParquetWriteConfig,
    ) -> Result<Vec<String>> {
        Self::export_to_parquet_with_registry(
            proteins,
            peptides,
            dir,
            export_all_psms,
            config,
            ModificationsDB::global(),
        )
    }

    /// Export protein and peptide identifications to an idparquet directory
    /// bundle; source `exportToParquet(protein_identifications,
    /// peptide_identifications, dir, export_all_psms = true, config = {})`.
    ///
    /// Writes (and overwrites) the four canonical files inside `dir`; other
    /// files in `dir` are left untouched. If `dir` does not exist it is
    /// created (only `dir` itself; the parent must already exist). In the
    /// source's order:
    ///
    /// 1. the run identifiers must be unique, checked before anything is
    ///    created, so no partial bundle is left behind;
    /// 2. `psms.parquet`: [`QpxFile::export_to_arrow`] (all hits per
    ///    identification with `export_all_psms`, which the source defaults
    ///    to `true` so that the bundle round-trips losslessly; else the best
    ///    hit per spectrum), written by [`write_table_to_parquet`];
    /// 3. `proteins.parquet`, `protein_groups.parquet` and
    ///    `search_params.parquet` by
    ///    [`ProteinIdentificationArrowIO`]'s exports, the last with the
    ///    modification definitions [`modification_definitions::collect`]
    ///    gathers from the runs and identifications (search-parameter names
    ///    resolved in `registry`), encoded per run by
    ///    [`modification_definitions::encode_by_run`].
    ///
    /// Returns the PSM export's warnings.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's message when two runs share
    /// an identifier (the source's `Exception::InvalidValue`, thrown before
    /// any file is opened) and, with the source's logged message
    /// "PSMArrowIO: target path exists and is not a directory: `dir`", when
    /// `dir` is a regular file; [`Error::Io`] when `dir` cannot be created
    /// ("PSMArrowIO: failed to create directory: `dir`"); and the errors of
    /// the exports and writes, which the source logs before returning
    /// `false`. A failure after the first file leaves the files already
    /// written, as in the source.
    pub fn export_to_parquet_with_registry(
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
        dir: impl AsRef<Path>,
        export_all_psms: bool,
        config: &ParquetWriteConfig,
        registry: &ModificationsDB,
    ) -> Result<Vec<String>> {
        let dir = dir.as_ref();
        // Mirror XMLHandler::checkUniqueIdentifiers_: fail before any file is
        // opened, so a partial bundle is never left behind.
        ProteinIdentificationArrowIO::check_unique_identifiers(proteins)?;
        ensure_directory(dir)?;

        let QpxOutput {
            value: psm_table,
            warnings,
        } = QpxFile::export_to_arrow(proteins, peptides, export_all_psms)?;
        write_table_to_parquet(&psm_table, dir.join(PSMS), config)?;

        ProteinIdentificationArrowIO::export_proteins_to_parquet(
            proteins,
            dir.join(PROTEINS),
            config,
        )?;
        ProteinIdentificationArrowIO::export_protein_groups_to_parquet(
            proteins,
            dir.join(PROTEIN_GROUPS),
            config,
        )?;
        let definitions = modification_definitions::collect(proteins, peptides, registry)?;
        let encoded = modification_definitions::encode_by_run(proteins, &definitions)?;
        ProteinIdentificationArrowIO::export_search_params_to_parquet(
            proteins,
            dir.join(SEARCH_PARAMS),
            config,
            &encoded,
        )?;
        Ok(warnings)
    }

    /// Import protein and peptide identifications from an idparquet
    /// directory bundle, resolving modifications against a copy of the
    /// global [`ModificationsDB`] that the bundle's definitions extend and
    /// that is then dropped.
    ///
    /// [`import_from_parquet_with_registry`](Self::import_from_parquet_with_registry)
    /// with that copy.
    ///
    /// # Errors
    ///
    /// As the registry form.
    pub fn import_from_parquet(
        dir: impl AsRef<Path>,
        proteins: &mut Vec<ProteinIdentification>,
        peptides: &mut Vec<PeptideIdentification>,
    ) -> Result<Vec<String>> {
        let mut registry = ModificationsDB::global().clone();
        Self::import_from_parquet_with_registry(dir, proteins, peptides, &mut registry)
    }

    /// Import protein and peptide identifications from an idparquet
    /// directory bundle; source `importFromParquet(dir,
    /// protein_identifications, peptide_identifications)`.
    ///
    /// All four canonical files must be present in `dir`. The runs come from
    /// the three protein-side files
    /// ([`ProteinIdentificationArrowIO::import_from_parquet_with_registry`],
    /// which registers the stored modification definitions in `registry`,
    /// as the source registers them process-wide), then the identifications
    /// from `psms.parquet` ([`QpxFile::import_from_arrow_with_registry`],
    /// with the stored run identifiers as join keys). Only then, as the idXML
    /// reader does on every load (`IdXMLFile.cpp:530`), every run gets a
    /// fresh identifier
    /// ([`ProteinIdentificationArrowIO::synthesize_run_identifiers`]) and the
    /// identifications are re-stamped in lock-step
    /// ([`ProteinIdentificationArrowIO::apply_run_identifier_rename`]).
    ///
    /// On success `proteins` and `peptides` hold exactly the bundle's data
    /// (the source swaps its results in); returns the imports' warnings.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] with the source's logged messages when `dir`
    /// is not a directory ("PSMArrowIO: not a directory: `dir`") or one of
    /// the four files is missing ("PSMArrowIO: missing required file:
    /// `dir`/`file`"), and the errors of the reads and imports. `proteins`,
    /// `peptides` and `registry` are unchanged on error (the source leaves
    /// its outputs unchanged too, but may already have registered
    /// definitions).
    pub fn import_from_parquet_with_registry(
        dir: impl AsRef<Path>,
        proteins: &mut Vec<ProteinIdentification>,
        peptides: &mut Vec<PeptideIdentification>,
        registry: &mut ModificationsDB,
    ) -> Result<Vec<String>> {
        let dir = dir.as_ref();
        if !dir.is_dir() {
            return Err(Error::InvalidValue(format!(
                "PSMArrowIO: not a directory: {}",
                dir.display()
            )));
        }
        for sub in [PSMS, PROTEINS, PROTEIN_GROUPS, SEARCH_PARAMS] {
            if !dir.join(sub).exists() {
                return Err(Error::InvalidValue(format!(
                    "PSMArrowIO: missing required file: {}",
                    dir.join(sub).display()
                )));
            }
        }

        let mut staged = registry.clone();
        let mut tmp_proteins = Vec::new();
        let mut tmp_peptides = Vec::new();
        let mut warnings = ProteinIdentificationArrowIO::import_from_parquet_with_registry(
            dir.join(PROTEINS),
            dir.join(PROTEIN_GROUPS),
            dir.join(SEARCH_PARAMS),
            &mut tmp_proteins,
            &mut staged,
        )?;
        let psm_table = ParquetFile::read_table(dir.join(PSMS))?;
        warnings.extend(QpxFile::import_from_arrow_with_registry(
            &psm_table,
            &mut tmp_proteins,
            &mut tmp_peptides,
            &mut staged,
        )?);

        // IdXMLFile.cpp:530 parity: fresh run identifiers on every load,
        // synthesized only after all four tables were joined on the stored
        // ones; the identifications follow in lock-step.
        let (rename, rename_warnings) =
            ProteinIdentificationArrowIO::synthesize_run_identifiers(&mut tmp_proteins);
        warnings.extend(rename_warnings);
        ProteinIdentificationArrowIO::apply_run_identifier_rename(&rename, &mut tmp_peptides);

        *proteins = tmp_proteins;
        *peptides = tmp_peptides;
        *registry = staged;
        Ok(warnings)
    }
}

/// The source's `ensureDirectory_`: an existing directory is used, a missing
/// one is created (not its parents), anything else is refused.
fn ensure_directory(dir: &Path) -> Result<()> {
    if dir.exists() {
        if !dir.is_dir() {
            return Err(Error::InvalidValue(format!(
                "PSMArrowIO: target path exists and is not a directory: {}",
                dir.display()
            )));
        }
        return Ok(());
    }
    std::fs::create_dir(dir).map_err(|error| {
        Error::Io(std::io::Error::new(
            error.kind(),
            format!(
                "PSMArrowIO: failed to create directory: {}: {error}",
                dir.display()
            ),
        ))
    })
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The conventions every search-engine adapter shares: the native form of the
//! `OpenMS4-cli` package's `APPLICATIONS/SearchEngineBase.h` and
//! `SearchEngineBase.cpp` (cli `c19e494`).
//!
//! The source class sits between `TOPPExternalToolBase` and the adapters
//! (`MSFraggerAdapter`, `MSGFPlusAdapter`, `SageAdapter` at topp `174b576`)
//! and bundles: the `-in` spectra check ([`SearchEngineBase::get_rawfile_name`]),
//! the `-database` resolution ([`SearchEngineBase::get_db_filename`]), and the
//! optional post-search reindexing with `PeptideIndexing`, whose parameters
//! [`SearchEngineBase::register_peptide_indexing_parameter`] registers. Its
//! external-process runner is the one every tool has,
//! [`ToolContext::run_external_process`].
//!
//! Not ported: `reindex_`. It sets up `PeptideIndexing` from the tool's
//! `PeptideIndexing:` parameters (`setParameters` on the indexer's
//! `getParameters()` defaults) and maps the indexer's `ExitCodes`; the native
//! [`PeptideIndexing`](crate::analysis::peptide_indexing::PeptideIndexing)
//! has neither a parameter surface nor exit codes yet, so there is nothing
//! for `reindex_` to call. See `docs/TOPP_CLI_SUPPORT.md` (*Tool base
//! classes*).
//!
//! The source's constructor forwards the tool name, description, `official`,
//! citations and `toolhandler_test` to `TOPPBase`; here those are the tool's
//! [`Tool`](super::Tool) implementation, and the deleted default and copy
//! constructors and the defaulted destructor have no counterpart.

use super::console::log_warning;
use super::context::ToolContext;
use super::outcome::ToolError;
use super::parameter::ExitCode;
use super::spec::ToolSpec;
use crate::format::file_handler::FileHandler;
use crate::format::file_types::FileType;
use crate::param::{Param, ParamValue};
use crate::system::file::{self, FileContext};
use crate::{Error, Result};
use std::io::Write;

/// The `PeptideIndexing` parameters the adapters register as advanced
/// (`SearchEngineBase.cpp:147-149`), in the source's order.
pub const PEPTIDE_INDEXING_ADVANCED: [&str; 12] = [
    "decoy_string",
    "decoy_string_position",
    "missing_decoy_action",
    "enzyme:name",
    "enzyme:specificity",
    "write_protein_sequence",
    "write_protein_description",
    "keep_unreferenced_proteins",
    "unmatched_action",
    "aaa_max",
    "mismatches_max",
    "IL_equivalent",
];

/// How many recognised spectra `getRawfileName` inspects: the default of
/// `MzMLFile::getCentroidInfo`'s `first_n_spectra_only`.
pub const RAWFILE_INSPECTED_SPECTRA: usize = 10;

/// The shared part of the search-engine adapters (source class
/// `SearchEngineBase`).
///
/// The source's members are protected methods of the adapter object; here
/// they are associated functions taking the run context, so an adapter's
/// [`Tool`](super::Tool) implementation calls them from its registration and
/// its body.
#[derive(Clone, Copy, Debug, Default)]
pub struct SearchEngineBase;

impl SearchEngineBase {
    /// Resolve `-in` and check that its spectra suit a search engine (source
    /// `getRawfileName`, `SearchEngineBase.cpp:31-89`), returning the `-in`
    /// path as given.
    ///
    /// The file's type decides ([`FileHandler::get_type`]):
    ///
    /// * **mzML**: the first [`RAWFILE_INSPECTED_SPECTRA`] recognised spectra
    ///   of MS level `ms_level` are classified
    ///   ([`centroid_info_with_options`](crate::format::mzml::centroid_info_with_options)
    ///   restricted to that level). With no spectrum of the level, the run ends
    ///   as the source's `FileEmpty` does. With profile spectra, or with no
    ///   centroided one, the run ends as the source's `IllegalArgument` does,
    ///   unless `-force` is set, in which case a warning is written instead
    ///   and the check goes on. Spectra of unknown type are not judged: that
    ///   type only arises for (almost) empty spectra without metadata.
    /// * **MGF** and **Bruker TDF**: accepted without inspection — MGF is
    ///   centroided by definition and timsTOF data inherently centroided.
    /// * anything else: the warning `Warning: make sure that MS<level> spectra
    ///   in '<in>' are centroided. Otherwise the results may be undefined!`.
    ///
    /// The warnings go to `err` as the source's `OPENMS_LOG_WARN` writes them
    /// (yellow on a terminal, not in the `-log` file). The source writes the
    /// last one without a line break, so its text stays in the log stream's
    /// buffer and appears in front of the next warning line, or at the end of
    /// the process; this writes it as a line of its own.
    ///
    /// # Errors
    ///
    /// * [`ExitCode::InputFileEmpty`], `Error: File empty (the file 'Error: No
    ///   MS<level> spectra in input file.' is empty)`: no spectrum of the
    ///   level (source `FileEmpty`, whose message frames the text it is given
    ///   as a file name);
    /// * [`ExitCode::UnknownError`], `Error: Unexpected internal error (Error:
    ///   Profile data provided but centroided MS<level> spectra required. …)`
    ///   or `(Error: No centroided MS<level> spectra were found, but are
    ///   required. …)` without `-force` (source `IllegalArgument`, which
    ///   reaches the `BaseException` arm);
    /// * the accessor error for an unregistered `in` or `force`, the type
    ///   detection's and the mzML reader's errors, which the lifecycle maps;
    ///   an mzML input without the `mzml` feature is [`Error::Unsupported`];
    /// * [`Error::Io`] when `err` fails.
    pub fn get_rawfile_name(
        ctx: &ToolContext,
        err: &mut dyn Write,
        ms_level: i32,
    ) -> std::result::Result<String, ToolError> {
        let inputfile_name = ctx.string("in")?.to_owned();
        match FileHandler::get_type(&inputfile_name)? {
            FileType::MzMl => {
                Self::check_mzml_centroided(ctx, err, &inputfile_name, ms_level)?;
            }
            FileType::Mgf | FileType::BrukerTdf => {}
            _ => {
                log_warning(
                    err,
                    &format!(
                        "Warning: make sure that MS{ms_level} spectra in '{inputfile_name}' are centroided. Otherwise the results may be undefined!"
                    ),
                )?;
            }
        }
        Ok(inputfile_name)
    }

    /// The mzML branch of [`get_rawfile_name`](Self::get_rawfile_name).
    #[cfg(feature = "mzml")]
    fn check_mzml_centroided(
        ctx: &ToolContext,
        err: &mut dyn Write,
        path: &str,
        ms_level: i32,
    ) -> std::result::Result<(), ToolError> {
        use crate::format::mzml::{
            CentroidInfoLimits, LoadOptions, ReadOptions, centroid_info_with_options,
        };
        let mut options = LoadOptions::default();
        options.scientific.set_ms_levels(&[ms_level])?;
        let centroid_info = centroid_info_with_options(
            path,
            RAWFILE_INSPECTED_SPECTRA,
            &options,
            &ReadOptions::default(),
            CentroidInfoLimits::default(),
        )?;
        // The source looks the level up as an unsigned key, so a negative
        // level finds nothing, as it does here.
        let level_info = u32::try_from(ms_level)
            .ok()
            .and_then(|level| centroid_info.get(&level));
        let Some(level_info) = level_info else {
            return Err(ToolError::caught(
                ExitCode::InputFileEmpty,
                format!(
                    "Error: File empty (the file 'Error: No MS{ms_level} spectra in input file.' is empty)"
                ),
            ));
        };
        let force = ctx.flag("force")?;
        if level_info.count_profile > 0 {
            if force {
                log_warning(
                    err,
                    "Warning: Profile data found, but centroid MS spectra required. Since '-force' flag is in effect, we will continue, but results are likely bogus.",
                )?;
            } else {
                return Err(ToolError::unexpected(format!(
                    "Error: Profile data provided but centroided MS{ms_level} spectra required. To enforce processing (unwise!) of the data enable the -force flag (results will be bogus!)."
                )));
            }
        }
        if level_info.count_centroided == 0 {
            if force {
                log_warning(
                    err,
                    &format!(
                        "Warning: No centroided MS{ms_level} were found, but are required. Since '-force' flag is in effect, we will continue, but results might be bogus."
                    ),
                )?;
            } else {
                return Err(ToolError::unexpected(format!(
                    "Error: No centroided MS{ms_level} spectra were found, but are required. To enforce processing of the data enable the -force flag (results will likely be bogus!)."
                )));
            }
        }
        Ok(())
    }

    /// The mzML branch without the `mzml` feature.
    #[cfg(not(feature = "mzml"))]
    fn check_mzml_centroided(
        _ctx: &ToolContext,
        _err: &mut dyn Write,
        _path: &str,
        _ms_level: i32,
    ) -> std::result::Result<(), ToolError> {
        Err(Error::Unsupported(
            "checking the spectra of an mzML input needs the mzml feature".into(),
        )
        .into())
    }

    /// Resolve the search database (source `getDBFilename`,
    /// `SearchEngineBase.cpp:91-99`): `db` when it is not empty, otherwise
    /// the `-database` value; a name that is not a readable file is looked up
    /// in the `id_db_dir` directories of `OpenMS.ini` and then in the shared
    /// data directory ([`FileContext::find_database`]).
    ///
    /// `db` stands in for `-database` for special database formats (the
    /// source names OMSSA's). `files` is the context the source's `File`
    /// functions read from the process: `OpenMS.ini` and the data path.
    ///
    /// When the lookup finds the name, the source's `File::findDatabase`
    /// writes `Augmenting database name '<db>' with path given in
    /// 'OpenMS.ini:id_db_dir'. Full name is now: '<full>'` to its info log,
    /// which is `out` here. The name returned is the file as given when it
    /// exists (readable or not), and otherwise the lexically normalised
    /// directory entry that holds it.
    ///
    /// # Errors
    ///
    /// When the lookup fails, the source's `findDatabase` writes `Input
    /// database '<db>' not found (<what>). Make sure it exists (and check
    /// 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!` to its
    /// error log (`err`, not the `-log` file; the parenthesis is the
    /// source's) and rethrows `FileNotFound`, whose catch block ends the run
    /// with [`ExitCode::InputFileNotFound`] and `Error: File not found
    /// (<what>)`. `<what>` is `the file '<db>' could not be found`, followed
    /// for a non-blank name by the source's hint naming the shared-data
    /// directory and where it was resolved from. Also the accessor error for
    /// an unregistered `database`, the errors of reading `OpenMS.ini` and of
    /// resolving the data path, [`Error::InvalidValue`] for a resolved path
    /// that is not UTF-8, and [`Error::Io`] when `out` or `err` fails.
    pub fn get_db_filename(
        ctx: &ToolContext,
        out: &mut dyn Write,
        err: &mut dyn Write,
        files: &FileContext,
        db: &str,
    ) -> std::result::Result<String, ToolError> {
        let db_name = if db.is_empty() {
            ctx.string("database")?.to_owned()
        } else {
            db.to_owned()
        };
        if file::readable(&db_name) {
            return Ok(db_name);
        }
        match files.find_database(&db_name) {
            Ok(full) => {
                let full = full
                    .to_str()
                    .ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "the database path found for '{db_name}' is not UTF-8"
                        ))
                    })?
                    .to_owned();
                writeln!(
                    out,
                    "Augmenting database name '{db_name}' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '{full}'"
                )
                .map_err(Error::Io)?;
                Ok(full)
            }
            Err(Error::Io(error)) if error.kind() == std::io::ErrorKind::NotFound => {
                let what = Self::database_not_found_text(files, &db_name)?;
                super::console::log_error(
                    err,
                    &format!(
                        "Input database '{db_name}' not found ({what}). Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!"
                    ),
                )
                .map_err(Error::Io)?;
                Err(ToolError::file_not_found(what))
            }
            Err(error) => Err(error.into()),
        }
    }

    /// The `what()` of the `FileNotFound` that `File::find` throws for
    /// `filename` (`File.cpp:584-587` and `:616-624`, message
    /// `Exception.cpp:166-167`).
    fn database_not_found_text(files: &FileContext, filename: &str) -> Result<String> {
        if filename.bytes().all(|byte| byte.is_ascii_whitespace()) {
            return Ok(format!("the file '{filename}' could not be found"));
        }
        let data = files.get_openms_data_path()?;
        let source = files.get_openms_data_path_source()?;
        Ok(format!(
            "the file '{filename}' could not be found. OpenMS searched its shared-data directory '{}' (via {source}). If this is a wrong or outdated OpenMS installation, reinstall the Core SDK with its shared data. OPENMS_DATA_PATH is an explicit override that takes precedence over installation discovery, so unset or correct it if it points to a stale location",
            data.display()
        ))
    }

    /// Register `-reindex` and the `PeptideIndexing` parameters (source
    /// `registerPeptideIndexingParameter_`, `SearchEngineBase.cpp:139-157`).
    ///
    /// `-reindex` is a string choice of `true` and `false`, default `true`:
    /// "Recalculate peptide to protein association using OpenMS. Annotates
    /// target-decoy information." Then `peptide_indexing_parameter` — the
    /// indexer's defaults, which an adapter may have adjusted (for example for
    /// non-tryptic digestion) — gets `missing_decoy_action` set to `warn`,
    /// replacing that entry with a bare one as the source's
    /// `Param::setValue` does (its description, restrictions and tags are
    /// dropped), and the twelve entries of [`PEPTIDE_INDEXING_ADVANCED`]
    /// tagged `advanced`. The whole tree is registered under the section
    /// `PeptideIndexing:`, so that its names cannot collide with the tool's
    /// own ([`ToolSpec::register_full_param`]).
    ///
    /// # Errors
    ///
    /// The registration errors of [`ToolSpec`], and the parameter tree's
    /// error when one of the twelve entries is missing from
    /// `peptide_indexing_parameter` (the source's `Param::addTag` throws
    /// `ElementNotFound` then).
    pub fn register_peptide_indexing_parameter(
        spec: &mut ToolSpec,
        mut peptide_indexing_parameter: Param,
    ) -> Result<()> {
        spec.register_string_option(
            "reindex",
            "<choice>",
            "true",
            "Recalculate peptide to protein association using OpenMS. Annotates target-decoy information.",
            false,
            false,
        )?;
        spec.set_valid_strings("reindex", &["true", "false"])?;

        peptide_indexing_parameter.set_value(
            "missing_decoy_action",
            ParamValue::String("warn".into()),
            "",
            &[],
        )?;
        for name in PEPTIDE_INDEXING_ADVANCED {
            peptide_indexing_parameter.add_tag(name, "advanced")?;
        }
        let mut combined = Param::new();
        combined.insert("PeptideIndexing:", &peptide_indexing_parameter)?;
        spec.register_full_param(&combined)
    }
}

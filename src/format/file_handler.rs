// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared experiment dispatch for TOPP-style pipelines, the native side of
//! `FORMAT/FileHandler.h`.
//!
//! Format recognition is separate from the available adapters, whose
//! documented representation limits apply. [`FileHandler::get_type`] ports the
//! source's name-then-content type detection, and the `*_with_options` loaders
//! hand the source's `PeakFileOptions` and `FeatureFileOptions` to the native
//! readers. `docs/MZML_MOBILITY_SUPPORT.md` records the API mapping and the
//! executed C++ comparison.

use super::file_types::{consistent_output_type, type_by_file_name};
use super::{FileType, PeakFileOptions, dta, mgf};
use crate::concept::progress_logger::ProgressLogger;
use crate::kernel::{
    ChromatogramConversionLimits, ChromatogramTools, MSChromatogram, MSSpectrum,
};
use crate::metadata::{ChecksumType, ScanMode, SourceFile};
use crate::{Error, MSExperiment, Result};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

/// Bytes read for bounded content recognition; [`type_by_content`] examines at
/// most this many.
const SNIFF_BYTES: u64 = 65_536;

/// The arguments of source `FileHandler::loadExperiment` beyond the file name,
/// the allowed types and the log type, together with the handler's
/// `PeakFileOptions` (`FileHandler.h:164-179`).
///
/// [`FileHandler::load_experiment_with_progress`] takes these; the default is
/// the source's default call: default options, no rewriting, no hash.
#[derive(Clone, Debug, Default)]
pub struct ExperimentLoadOptions {
    /// The handler's `getOptions()`, handed to the readers the source hands
    /// them to: mzML, mzXML and mzData take all they consume, DTA2D its three
    /// ranges, and DTA, MGF and MS2 none.
    pub peaks: PeakFileOptions,
    /// The native mzML reader's limits and source-compatibility switches, with
    /// no source counterpart; the strict default unless a caller opts in.
    #[cfg(feature = "mzml")]
    pub mzml: super::mzml::ReadOptions,
    /// Source `rewrite_source_file`: replace the experiment's source files by
    /// one record naming the loaded file, as
    /// [`FileHandler::load_experiment_with_progress`] describes.
    pub rewrite_source_file: bool,
    /// Source `compute_hash`: with `rewrite_source_file`, record the loaded
    /// file's SHA-1 ([`FileHandler::compute_file_hash`]) in that record.
    pub compute_hash: bool,
}

/// Native reader/writer dispatch; unsupported formats return explicit errors.
pub struct FileHandler;

impl FileHandler {
    /// Whether an identification adapter is compiled into this crate.
    pub fn can_read_identifications(kind: FileType) -> bool {
        kind == FileType::IdXml && cfg!(feature = "idxml")
    }
    /// Whether an identification writer is compiled into this crate; the same
    /// adapter set as [`FileHandler::can_read_identifications`].
    pub fn can_write_identifications(kind: FileType) -> bool {
        Self::can_read_identifications(kind)
    }

    /// Identification stream dispatch; remaining identification formats are errors.
    #[cfg(feature = "idxml")]
    pub fn read_identifications(
        reader: impl BufRead,
        kind: FileType,
    ) -> Result<super::idxml::IdXmlDocument> {
        if kind != FileType::IdXml {
            return Err(Error::Unsupported(format!(
                "native {} identification input",
                kind.name()
            )));
        }
        super::idxml::read(reader)
    }

    /// Write an identification document to an explicitly typed stream.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] for every kind but idXML, the only native
    /// identification writer, and the writer's error otherwise.
    #[cfg(feature = "idxml")]
    pub fn write_identifications(
        writer: impl Write,
        document: &super::idxml::IdXmlDocument,
        kind: FileType,
    ) -> Result<()> {
        if kind != FileType::IdXml {
            return Err(Error::Unsupported(format!(
                "native {} identification output",
                kind.name()
            )));
        }
        super::idxml::write(writer, document)
    }

    /// Known extensions take precedence; unknown extensions use bounded content
    /// recognition. An empty allowed list accepts every compiled adapter.
    #[cfg(feature = "idxml")]
    pub fn load_identifications(
        path: impl AsRef<Path>,
        allowed: &[FileType],
    ) -> Result<super::idxml::IdXmlDocument> {
        let reader = typed_input(path.as_ref(), allowed, FileType::IdXml)?;
        Self::read_identifications(reader, FileType::IdXml)
    }

    /// Source identification output selection: a single allowed format supplies
    /// an unknown suffix; known suffixes must belong to a nonempty allowed list.
    /// IdXMLFile stores plain bytes even when a compression suffix is present.
    #[cfg(feature = "idxml")]
    pub fn store_identifications(
        path: impl AsRef<Path>,
        document: &super::idxml::IdXmlDocument,
        allowed: &[FileType],
    ) -> Result<()> {
        let path = path.as_ref();
        let mut kind = type_by_file_name(filename(path)?);
        if kind == FileType::Unknown && allowed.len() == 1 {
            kind = allowed[0];
        }
        if !allowed.is_empty() && !allowed.contains(&kind) {
            return Err(Error::InvalidValue(format!(
                "{} is not an allowed output format",
                kind.name()
            )));
        }
        if !Self::can_write_identifications(kind) {
            return Err(Error::Unsupported(format!(
                "native {} identification output",
                kind.name()
            )));
        }
        super::idxml::store(path, document)
    }

    /// Whether a native experiment reader exists for `kind`: DTA, MGF, DTA2D
    /// and MS2 always, mzML, mzXML and mzData with the `mzml` feature.
    ///
    /// This is native availability in this build, not the source capabilities
    /// that [`FileType::source_properties`] lists.
    pub fn can_read_experiment(kind: FileType) -> bool {
        matches!(
            kind,
            FileType::Dta | FileType::Mgf | FileType::Dta2d | FileType::Ms2
        ) || (matches!(kind, FileType::MzMl | FileType::MzXml | FileType::MzData)
            && cfg!(feature = "mzml"))
    }

    /// Whether a native experiment writer exists for `kind`: DTA, MGF, DTA2D
    /// and MS2 always, mzML with the `mzml` feature. mzXML and mzData are read
    /// but not yet written through this dispatch.
    pub fn can_write_experiment(kind: FileType) -> bool {
        matches!(
            kind,
            FileType::Dta | FileType::Mgf | FileType::Dta2d | FileType::Ms2
        ) || (kind == FileType::MzMl && cfg!(feature = "mzml"))
    }

    /// Read an explicitly typed, uncompressed stream with the native adapter.
    pub fn read_experiment(reader: impl BufRead, kind: FileType) -> Result<MSExperiment> {
        match kind {
            FileType::Dta => Ok(MSExperiment {
                spectra: vec![dta::read(reader)?],
                ..Default::default()
            }),
            FileType::Mgf => mgf::read(reader),
            FileType::Dta2d => super::dta2d::read(reader),
            FileType::Ms2 => super::ms2::read(reader),
            #[cfg(feature = "mzml")]
            FileType::MzMl => super::mzml::read(reader),
            #[cfg(feature = "mzml")]
            FileType::MzXml => super::mzxml::read(reader),
            #[cfg(feature = "mzml")]
            FileType::MzData => super::mzdata::read(reader),
            _ => Err(unsupported(kind, "reading")),
        }
    }

    /// Write a typed stream. The selected adapter validates before emitting data;
    /// underlying I/O failures may leave partial output on a caller-owned stream.
    pub fn write_experiment(
        writer: impl Write,
        experiment: &MSExperiment,
        kind: FileType,
    ) -> Result<()> {
        match kind {
            FileType::Dta => {
                if experiment.settings.has_transport_metadata() {
                    return Err(Error::Unsupported(
                        "DTA cannot store experiment settings".into(),
                    ));
                }
                experiment.validate()?;
                if experiment.spectra.len() != 1 || !experiment.chromatograms.is_empty() {
                    return Err(Error::Unsupported("DTA experiment output requires exactly one spectrum and no experiment metadata or chromatograms".into()));
                }
                dta::write(writer, &experiment.spectra[0])
            }
            FileType::Mgf => mgf::write(writer, experiment),
            FileType::Dta2d => super::dta2d::write(writer, experiment),
            FileType::Ms2 => super::ms2::write(writer, experiment),
            #[cfg(feature = "mzml")]
            FileType::MzMl => super::mzml::write(writer, experiment),
            _ => Err(unsupported(kind, "writing")),
        }
    }

    /// Tries to determine the file type, by name first and then by content.
    ///
    /// Source `FileHandler::getType` (`FileHandler.cpp:204-251`). Trailing `/`
    /// and `\` separators are stripped first, because shell completion appends
    /// one to directory formats such as Bruker `.d`. A recognised extension
    /// wins without touching the file system, so a missing file with a known
    /// extension still has that type.
    ///
    /// A Bruker TDF name (`.d`, including `.d.zip`) yields
    /// [`FileType::Unknown`]. The source returns that unless it was built
    /// `WITH_OPENTIMS`, and then only when the marker files exist; this crate
    /// has no TDF reader, and the product-sdk oracle is built without OpenTIMS
    /// too, so the port follows that branch.
    ///
    /// Every other unknown name is resolved by [`type_by_content`] over at most
    /// 64 KiB of the file, opened through the loaders' gzip/bzip2 detection.
    /// The source reads five lines (8 KiB when compressed) plus 512 lines for
    /// imzML markers; the byte bound is native and makes the work independent
    /// of line length.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a non-UTF-8 path, and the I/O error
    /// when content recognition cannot open or read the file, where the source
    /// throws `Exception::FileNotFound`. A gzip or bzip2 input without the
    /// `file-compression` feature, or a ZIP container, is
    /// [`Error::Unsupported`]; the source looks inside all three.
    pub fn get_type(path: impl AsRef<Path>) -> Result<FileType> {
        let normalized = filename(path.as_ref())?.trim_end_matches(['/', '\\']);
        let kind = type_by_file_name(normalized);
        if kind == FileType::BrukerTdf {
            return Ok(FileType::Unknown);
        }
        if kind != FileType::Unknown {
            return Ok(kind);
        }
        let mut preview = Vec::new();
        super::path_io::open(Path::new(normalized))?
            .take(SNIFF_BYTES)
            .read_to_end(&mut preview)?;
        Ok(type_by_content(&preview))
    }

    /// Load by extension, falling back to bounded content recognition only for
    /// unknown extensions. An empty allowed list accepts every available adapter.
    /// Gzip and bzip2 containers use the optional `file-compression` feature.
    /// mzML keeps the stream reader's input order; use
    /// [`FileHandler::load_experiment_with_options`] for source loading options.
    ///
    /// An mzML input then has its SRM spectra moved into chromatograms, as
    /// every source `loadExperiment` does (`FileHandler.cpp:910`); see
    /// [`FileHandler::load_experiment_with_progress`].
    pub fn load_experiment(path: impl AsRef<Path>, allowed: &[FileType]) -> Result<MSExperiment> {
        let path = path.as_ref();
        let kind = type_by_file_name(path.to_str().unwrap_or(""));
        let mut document = crate::metadata::DocumentIdentifier::new();
        document.set_loaded_file_path(filename(path)?)?;
        document.set_loaded_file_type(path)?;
        let mut result = load_stream(super::path_io::open(path)?, kind, allowed)?;
        result.settings.document.loaded_file_path = document.loaded_file_path;
        result.settings.document.loaded_file_type = document.loaded_file_type;
        Ok(result)
    }

    /// Load an experiment with explicit source `PeakFileOptions`.
    ///
    /// Source `FileHandler::loadExperiment` after `getOptions() = options`
    /// (`FileHandler.cpp:849-1052`), as FeatureFinderCentroided loads MS1
    /// spectra with positive intensities. The type comes from the file name,
    /// else from bounded content recognition, and must be in a nonempty
    /// `allowed` list; an empty list accepts every native adapter. The options
    /// reach the readers the source hands them to:
    ///
    /// - mzML: every option, through `mzml::read_with_load_options` with the
    ///   default, strict `mzml::ReadOptions`: default XML and binary limits and
    ///   no source-compatibility switch. `load_experiment_with_read_options`
    ///   takes explicit ones. Unlike [`FileHandler::load_experiment`] this
    ///   applies the source's scientific defaults, including sorting peaks by
    ///   m/z.
    /// - mzXML and mzData: every option their readers consume
    ///   (`docs/MZXML_SUPPORT.md`, `docs/MZDATA_SUPPORT.md`), with default
    ///   resource limits.
    /// - DTA2D: the retention-time, m/z and intensity ranges, the only options
    ///   `DTA2DFile::load` consumes.
    /// - DTA, MGF and MS2: none, exactly as the source ignores them there.
    ///
    /// Ranges are half-open like `DRange::encloses`: a value equal to the
    /// minimum is kept and one equal to the maximum is dropped, so the
    /// FeatureFinderCentroided range `[f64::MIN_POSITIVE, f64::MAX)` removes
    /// zero and negative intensities. The source keeps a NaN intensity, because
    /// both of its comparisons fail; the native mzML and DTA2D readers reject
    /// non-finite values before any filter runs.
    ///
    /// The source's mzML branch finally moves SRM spectra into chromatograms
    /// (`ChromatogramTools::convertSpectraToChromatograms`, `:910`), and so
    /// does this loader, as [`FileHandler::load_experiment_with_progress`]
    /// describes. No progress is reported and the source files are kept as
    /// read; that entry point takes the source's log type,
    /// `rewrite_source_file` and `compute_hash`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the detected type is not in a
    /// nonempty `allowed` list, where the source throws
    /// `Exception::ParseError`; [`Error::Unsupported`] for a type without a
    /// native adapter; and the reader's error for malformed input, options the
    /// reader cannot execute, DTA2D ranges without finite ordered endpoints, or
    /// exceeded limits.
    pub fn load_experiment_with_options(
        path: impl AsRef<Path>,
        allowed: &[FileType],
        options: &PeakFileOptions,
    ) -> Result<MSExperiment> {
        let options = ExperimentLoadOptions {
            peaks: options.clone(),
            ..Default::default()
        };
        load_experiment_from(path.as_ref(), allowed, &options, None)
    }

    /// Load an experiment as source `FileHandler::loadExperiment` does with
    /// all of its arguments (`FileHandler.cpp:849-1052`).
    ///
    /// The type is detected, checked against a nonempty `allowed` list and
    /// read with `options.peaks` exactly as for
    /// [`FileHandler::load_experiment_with_options`]. Then:
    ///
    /// - **Log type.** The source hands each file object only
    ///   `setLogType(log)`. Here the reader reports to a copy of `logger` made
    ///   as [`ProgressLogger::clone`] makes one: a fresh backend of its type
    ///   that shares its clock, nesting and GUI factory. mzML, mzXML, mzData
    ///   and DTA2D make their readers' calls (`mzml::load_with_progress`,
    ///   `mzxml::load_with_progress`, `mzdata::load_with_progress`,
    ///   `dta2d::read_with_progress`). DTA and MS2 make none, as `DTAFile` and
    ///   `MS2File::load` make none in the source. MGF makes none either, where
    ///   the source's `MascotGenericFile::load` reports its section: the
    ///   native MGF dispatch is a different reader, and the support document
    ///   records the gap. A `.dta2d` file that cannot be opened fails before
    ///   its section starts, where the source starts it first.
    /// - **SRM conversion.** An mzML input then has its SRM spectra moved into
    ///   chromatograms: `ChromatogramTools().convertSpectraToChromatograms(exp,
    ///   true)` (`:910`). Spectra whose scan mode is
    ///   [`ScanMode::SelectedReactionMonitoring`] are removed, and those with
    ///   one precursor and at least one peak become one chromatogram per
    ///   precursor and product m/z. Each removed SRM spectrum that cannot be
    ///   converted writes the source's `ChromatogramTools: need exactly one
    ///   precursor ...` line to the warning log. An experiment without SRM
    ///   spectra is untouched. The conversion's own output is bounded by
    ///   [`ChromatogramConversionLimits`]' defaults, raised only by what moving
    ///   the already loaded records costs.
    /// - **`rewrite_source_file`** (`:1014-1051`). The first source file is
    ///   kept, or a default one when there is none, with the warning `No source
    ///   file annotated.`; more than one warns `Expecting a single source file
    ///   in mzML. Found N will take only first one for rewriting.` Its name
    ///   becomes the basename of `path` after trailing `/` and `\` are
    ///   stripped, its path `file://` followed by the directory of
    ///   [`crate::system::file::absolute_path`] of that name (`file:///` when
    ///   the directory does not start with `/`), and its file type the mzML
    ///   term of the detected type ([`FileType::mzml_name`]). Everything else
    ///   the record carried, its checksum included, stays. It then replaces all
    ///   source files.
    /// - **`compute_hash`**: that record's checksum becomes
    ///   [`FileHandler::compute_file_hash`] of `path` with type
    ///   [`ChecksumType::Sha1`].
    ///
    /// Warnings go to the global warning log
    /// ([`crate::concept::log_stream::LogLevel::Warn`]), as the source's go to
    /// `OPENMS_LOG_WARN`.
    ///
    /// # Errors
    ///
    /// As [`FileHandler::load_experiment_with_options`], plus the progress
    /// calls' errors. [`Error::InvalidValue`] when `compute_hash` is set
    /// without `rewrite_source_file`, the source's debug-build precondition
    /// (a release build ignores the flag then), or when the absolute path is
    /// not UTF-8. [`Error::Unsupported`] for `compute_hash` without the `mzml`
    /// feature, which carries the SHA-1 implementation. The SRM conversion's
    /// errors for non-finite values or exceeded ceilings; the source has
    /// neither check.
    pub fn load_experiment_with_progress(
        path: impl AsRef<Path>,
        allowed: &[FileType],
        options: &ExperimentLoadOptions,
        logger: &ProgressLogger,
    ) -> Result<MSExperiment> {
        load_experiment_from(path.as_ref(), allowed, options, Some(logger))
    }

    /// The SHA-1 of the file's bytes as 40 lowercase hexadecimal digits.
    ///
    /// Source `FileHandler::computeFileHash` (`FileHandler.cpp:768-778`),
    /// which hashes the file as stored, so a compressed file's digest is that
    /// of its compressed bytes. The source returns the empty string when the
    /// file cannot be opened, and so does this; it also returns the empty
    /// string when a read fails part way, where the source returns the digest
    /// of the bytes read so far.
    #[cfg(feature = "mzml")]
    pub fn compute_file_hash(path: impl AsRef<Path>) -> String {
        use sha1::{Digest, Sha1};
        let Ok(mut file) = std::fs::File::open(path) else {
            return String::new();
        };
        let mut hasher = Sha1::new();
        let mut buffer = [0u8; 8192];
        loop {
            match file.read(&mut buffer) {
                Ok(0) => break,
                Ok(read) => hasher.update(&buffer[..read]),
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(_) => return String::new(),
            }
        }
        hasher
            .finalize()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    /// Load an experiment with explicit source `PeakFileOptions` and explicit
    /// mzML reader options.
    ///
    /// As [`FileHandler::load_experiment_with_options`], except that an mzML
    /// input is read with `read` instead of the strict default
    /// [`crate::format::mzml::ReadOptions`]; the other formats ignore it. A
    /// TOPP tool that reproduces the source's `FileHandler::loadExperiment`
    /// passes its source-compatibility switches here, such as
    /// [`crate::format::mzml::ReadOptions::source_dangling_references`], while
    /// every library default stays strict (decision D10 of the early TOPP
    /// bundle).
    ///
    /// # Errors
    ///
    /// As [`FileHandler::load_experiment_with_options`].
    #[cfg(feature = "mzml")]
    pub fn load_experiment_with_read_options(
        path: impl AsRef<Path>,
        allowed: &[FileType],
        options: &PeakFileOptions,
        read: &super::mzml::ReadOptions,
    ) -> Result<MSExperiment> {
        let options = ExperimentLoadOptions {
            peaks: options.clone(),
            mzml: read.clone(),
            ..Default::default()
        };
        load_experiment_from(path.as_ref(), allowed, &options, None)
    }

    /// Store to a sibling temporary file, replacing the destination only after
    /// successful serialization and flush. Invalid data never truncates an
    /// existing output. `requested` can select a format for an unknown suffix.
    pub fn store_experiment(
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
        requested: Option<FileType>,
    ) -> Result<()> {
        let path = path.as_ref();
        let path_text = path
            .to_str()
            .ok_or_else(|| Error::InvalidValue("output filename must be UTF-8".into()))?;
        let kind = consistent_output_type(path_text, requested.map(FileType::name).unwrap_or(""));
        if !Self::can_write_experiment(kind) {
            return Err(unsupported(kind, "writing"));
        }
        super::path_io::write(path, |writer| {
            Self::write_experiment(writer, experiment, kind)
        })
    }

    /// Read a native feature map. Other feature-map formats remain explicit errors.
    pub fn read_feature_map(
        reader: impl BufRead,
        kind: FileType,
    ) -> Result<crate::kernel::FeatureMap> {
        #[cfg(feature = "featurexml")]
        if kind == FileType::FeatureXml {
            return super::featurexml::read(reader);
        }
        let _ = reader;
        Err(Error::Unsupported(format!(
            "native {} feature-map input",
            kind.name()
        )))
    }

    /// Write a feature map to an explicitly typed stream.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] unless `kind` is featureXML and the
    /// `featurexml` feature is enabled, and the writer's error otherwise.
    pub fn write_feature_map(
        writer: impl Write,
        map: &crate::kernel::FeatureMap,
        kind: FileType,
    ) -> Result<()> {
        #[cfg(feature = "featurexml")]
        if kind == FileType::FeatureXml {
            return super::featurexml::write(writer, map);
        }
        let _ = (writer, map);
        Err(Error::Unsupported(format!(
            "native {} feature-map output",
            kind.name()
        )))
    }

    /// Load a feature map by extension, else by bounded content recognition,
    /// with default reader options, and record the loaded path and type.
    ///
    /// An empty `allowed` list accepts any detected type, but only featureXML
    /// has a native reader.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the type is not in a nonempty
    /// `allowed` list, [`Error::Unsupported`] for any type but featureXML, and
    /// the reader's error for malformed input.
    pub fn load_feature_map(
        path: impl AsRef<Path>,
        allowed: &[FileType],
    ) -> Result<crate::kernel::FeatureMap> {
        let path = path.as_ref();
        let reader = typed_input(path, allowed, FileType::FeatureXml)?;
        let mut map = Self::read_feature_map(reader, FileType::FeatureXml)?;
        map.loaded_file_path = filename(path)?.into();
        map.loaded_file_type = FileType::FeatureXml;
        Ok(map)
    }

    /// Load a feature map with explicit source `FeatureFileOptions`.
    ///
    /// Source `FileHandler::loadFeatures` after `getFeatOptions()` is set
    /// (`FileHandler.cpp:1253-1280`); FileInfo turns convex hulls and
    /// subordinates off this way (`FORMAT/FileInfo.cpp:1080-1084`). The type is
    /// detected as for [`FileHandler::load_feature_map`] and must be
    /// featureXML, the only native feature-map reader; the options then reach
    /// that reader unchanged, with its default resource limits. The source call
    /// also takes an allowed-type list, which cannot select anything else here.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] for any detected type but featureXML,
    /// where the source would also dispatch TSV, SpecArray, Kroenik and OMS
    /// maps, and the reader's error for malformed input or exceeded limits.
    #[cfg(feature = "featurexml")]
    pub fn load_feature_map_with_options(
        path: impl AsRef<Path>,
        options: &super::featurexml::FeatureFileOptions,
    ) -> Result<crate::kernel::FeatureMap> {
        let path = path.as_ref();
        let reader = typed_input(path, &[], FileType::FeatureXml)?;
        let mut map = super::featurexml::read_with_options(
            reader,
            &super::featurexml::ReadOptions {
                feature_options: options.clone(),
                ..Default::default()
            },
        )?;
        map.loaded_file_path = filename(path)?.into();
        map.loaded_file_type = FileType::FeatureXml;
        Ok(map)
    }

    /// Atomically store a feature map. `requested` selects the type for an
    /// unknown suffix and must agree with a known one, as source
    /// `getConsistentOutputfileType`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] unless the resolved type is featureXML
    /// with the `featurexml` feature, and the writer's or I/O error otherwise;
    /// an existing destination is left unchanged on failure.
    pub fn store_feature_map(
        path: impl AsRef<Path>,
        map: &crate::kernel::FeatureMap,
        requested: Option<FileType>,
    ) -> Result<()> {
        let path = path.as_ref();
        let kind =
            consistent_output_type(filename(path)?, requested.map(FileType::name).unwrap_or(""));
        super::path_io::write(path, |writer| Self::write_feature_map(writer, map, kind))
    }

    /// Read a consensus map from an explicitly typed stream.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] unless `kind` is consensusXML and the
    /// `consensusxml` feature is enabled, and the reader's error otherwise.
    pub fn read_consensus_map(
        reader: impl BufRead,
        kind: FileType,
    ) -> Result<crate::kernel::ConsensusMap> {
        #[cfg(feature = "consensusxml")]
        if kind == FileType::ConsensusXml {
            return super::consensusxml::read(reader);
        }
        let _ = reader;
        Err(Error::Unsupported(format!(
            "native {} consensus-map input",
            kind.name()
        )))
    }

    /// Write a consensus map to an explicitly typed stream.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] unless `kind` is consensusXML and the
    /// `consensusxml` feature is enabled, and the writer's error otherwise.
    pub fn write_consensus_map(
        writer: impl Write,
        map: &crate::kernel::ConsensusMap,
        kind: FileType,
    ) -> Result<()> {
        #[cfg(feature = "consensusxml")]
        if kind == FileType::ConsensusXml {
            return super::consensusxml::write(writer, map);
        }
        let _ = (writer, map);
        Err(Error::Unsupported(format!(
            "native {} consensus-map output",
            kind.name()
        )))
    }

    /// Load a consensus map by extension, else by bounded content recognition,
    /// and record the loaded path and type.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the type is not in a nonempty
    /// `allowed` list, [`Error::Unsupported`] for any type but consensusXML,
    /// and the reader's error for malformed input.
    pub fn load_consensus_map(
        path: impl AsRef<Path>,
        allowed: &[FileType],
    ) -> Result<crate::kernel::ConsensusMap> {
        let path = path.as_ref();
        let reader = typed_input(path, allowed, FileType::ConsensusXml)?;
        let mut map = Self::read_consensus_map(reader, FileType::ConsensusXml)?;
        map.loaded_file_path = filename(path)?.into();
        map.loaded_file_type = FileType::ConsensusXml;
        Ok(map)
    }

    /// Atomically store a consensus map, resolving the type as
    /// [`FileHandler::store_feature_map`] does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] unless the resolved type is consensusXML
    /// with the `consensusxml` feature, and the writer's or I/O error
    /// otherwise; an existing destination is left unchanged on failure.
    pub fn store_consensus_map(
        path: impl AsRef<Path>,
        map: &crate::kernel::ConsensusMap,
        requested: Option<FileType>,
    ) -> Result<()> {
        let path = path.as_ref();
        let kind =
            consistent_output_type(filename(path)?, requested.map(FileType::name).unwrap_or(""));
        super::path_io::write(path, |writer| Self::write_consensus_map(writer, map, kind))
    }
}

fn unsupported(kind: FileType, operation: &str) -> Error {
    Error::Unsupported(format!("native {} experiment {operation}", kind.name()))
}

/// The three `PeakFileOptions` ranges `DTA2DFile::load` consumes, as the
/// DTA2D reader's half-open ranges.
fn dta2d_options(options: &PeakFileOptions) -> super::dta2d::ReadOptions {
    let range =
        |present: bool, range: crate::kernel::NumericRange| present.then_some(range.min..range.max);
    super::dta2d::ReadOptions {
        rt_range: range(options.has_rt_range(), options.rt_range()),
        mz_range: range(options.has_mz_range(), options.mz_range()),
        intensity_range: range(options.has_intensity_range(), options.intensity_range()),
        ..Default::default()
    }
}

fn load_stream(
    mut reader: impl BufRead,
    kind: FileType,
    allowed: &[FileType],
) -> Result<MSExperiment> {
    let (kind, mut experiment) = if kind != FileType::Unknown {
        check_allowed(kind, allowed)?;
        (kind, FileHandler::read_experiment(reader, kind)?)
    } else {
        let mut preview = Vec::new();
        // The same bytes are replayed to the parser. This also works with nonseekable
        // decompression streams and caps both sniffing memory and decompression work.
        reader
            .by_ref()
            .take(SNIFF_BYTES)
            .read_to_end(&mut preview)?;
        let kind = type_by_content(&preview);
        check_allowed(kind, allowed)?;
        let experiment =
            FileHandler::read_experiment(BufReader::new(preview.as_slice().chain(reader)), kind)?;
        (kind, experiment)
    };
    if kind == FileType::MzMl {
        convert_srm_spectra(&mut experiment)?;
    }
    Ok(experiment)
}

fn check_allowed(kind: FileType, allowed: &[FileType]) -> Result<()> {
    if !allowed.is_empty() && !allowed.contains(&kind) {
        return Err(Error::InvalidValue(format!(
            "{} is not an allowed input format",
            kind.name()
        )));
    }
    Ok(())
}

/// Bounded content recognition for common SDK interchange types. This is a
/// heuristic, not validation; authoritative syntax checks belong to the reader.
/// It examines at most 64 KiB (five lines normally, 512 for IMS markers).
pub fn type_by_content(bytes: &[u8]) -> FileType {
    let bytes = &bytes[..bytes.len().min(65_536)];
    if bytes.starts_with(b"\x89PNG\r\n\x1a\n") {
        return FileType::Png;
    }
    let text = String::from_utf8_lossy(bytes);
    let lines: Vec<_> = text.lines().take(512).map(str::trim).collect();
    let first = lines.iter().take(5).copied().collect::<Vec<_>>().join(" ");
    if first.contains("<mzML") {
        let imaging = lines.iter().any(|s| {
            [
                "Imaging MS Ontology",
                "IMS:1000050",
                "IMS:1000030",
                "IMS:1000080",
            ]
            .iter()
            .any(|m| s.contains(m))
        });
        return if imaging {
            FileType::ImzMl
        } else {
            FileType::MzMl
        };
    }
    for (marker, kind) in [
        ("<mzXML", FileType::MzXml),
        ("<mzData", FileType::MzData),
        ("<MzIdentML", FileType::MzIdentMl),
        ("<MzQualityMLType", FileType::QcMl),
        (
            "http://regis-web.systemsbiology.net/pepXML",
            FileType::PepXml,
        ),
        (
            "http://regis-web.systemsbiology.net/protXML",
            FileType::ProtXml,
        ),
        ("<featureMap", FileType::FeatureXml),
        ("<IdXML", FileType::IdXml),
        ("<consensusXML", FileType::ConsensusXml),
        ("<TrafoXML", FileType::TransformationXml),
        ("<GelML", FileType::GelMl),
        ("<TraML", FileType::TraMl),
        ("<MSResponse", FileType::OmssaXml),
        ("<mascot_search_results", FileType::MascotXml),
    ] {
        if first.contains(marker) {
            return kind;
        }
    }
    if first.contains("<PARAMETERS") {
        return if first.contains("<NODE name=\"info\"")
            && first.contains("<ITEM name=\"num_vertices\"")
        {
            FileType::Toppas
        } else {
            FileType::Ini
        };
    }
    if first.starts_with('{') {
        return FileType::Json;
    }
    for line in lines.iter().take(5) {
        if line.starts_with('>') {
            return FileType::Fasta;
        }
        if !line.starts_with('#') {
            break;
        }
    }
    if lines
        .iter()
        .take(5)
        .any(|s| s.starts_with("Num peaks: ") || (s.starts_with("Name: ") && s.contains('/')))
    {
        return FileType::Msp;
    }
    let numbers: Vec<_> = lines
        .iter()
        .take(5)
        .skip(1)
        .flat_map(|s| s.split_whitespace())
        .collect();
    if numbers
        .iter()
        .all(|s| s.parse::<f64>().is_ok_and(f64::is_finite))
    {
        if numbers.len() == 8 {
            return FileType::Dta;
        }
        if numbers.len() == 12 {
            return FileType::Dta2d;
        }
    }
    if lines
        .iter()
        .take(5)
        .any(|s| *s == "BEGIN IONS" || *s == "FORMAT=Mascot generic")
    {
        return FileType::Mgf;
    }
    if first.starts_with('H') && first.contains("CreationDate") {
        return FileType::Ms2;
    }
    for (marker, kind) in [
        ("MTD\tmzTab-version", FileType::MzTab),
        (
            "scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription",
            FileType::Tsv,
        ),
    ] {
        if lines.iter().take(5).any(|line| line.contains(marker)) {
            return kind;
        }
    }
    if let Some(line) = lines.first() {
        // Leading indentation is already removed by the bounded line reader.
        if line.contains("m/z\t     rt(min)\t       snr\t      charge\t   intensity") {
            return FileType::Peplist;
        }
        if line.contains("File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\tBest Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications") {
            return FileType::Kroenik;
        }
        if line.starts_with("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds") {
            return FileType::Psms;
        }
    }
    FileType::Unknown
}

fn filename(path: &Path) -> Result<&str> {
    path.to_str()
        .ok_or_else(|| Error::InvalidValue("filename must be UTF-8".into()))
}

/// The body of [`FileHandler::load_experiment_with_options`],
/// [`FileHandler::load_experiment_with_read_options`] and, with a logger,
/// [`FileHandler::load_experiment_with_progress`].
fn load_experiment_from(
    path: &Path,
    allowed: &[FileType],
    options: &ExperimentLoadOptions,
    logger: Option<&ProgressLogger>,
) -> Result<MSExperiment> {
    let text = filename(path)?;
    if options.compute_hash && !options.rewrite_source_file {
        return Err(Error::InvalidValue(
            "compute_hash requires rewrite_source_file: the checksum is stored in the rewritten source file".into(),
        ));
    }
    if options.compute_hash && cfg!(not(feature = "mzml")) {
        return Err(Error::Unsupported(
            "source-file SHA-1 checksums need the `mzml` feature".into(),
        ));
    }
    let mut document = crate::metadata::DocumentIdentifier::new();
    document.set_loaded_file_path(text)?;
    document.set_loaded_file_type(path)?;
    let (kind, reader) = detect(path)?;
    check_allowed(kind, allowed)?;
    // Source `f.setLogType(log)` on a fresh file object: the reader reports to
    // a copy of the caller's logger with a fresh backend of the same type.
    let mut file_logger = logger.cloned();
    let mut result = match kind {
        #[cfg(feature = "mzml")]
        FileType::MzMl => {
            let load = super::mzml::LoadOptions {
                scientific: options.peaks.clone(),
                ..Default::default()
            };
            match file_logger.as_mut() {
                Some(logger) => {
                    drop(reader);
                    super::mzml::load_with_progress(path, &load, &options.mzml, logger)?
                }
                None => super::mzml::read_with_load_options(reader, &load, &options.mzml)?,
            }
        }
        #[cfg(feature = "mzml")]
        FileType::MzXml => {
            let read = super::mzxml::ReadOptions {
                peaks: options.peaks.clone(),
                ..Default::default()
            };
            match file_logger.as_mut() {
                Some(logger) => {
                    drop(reader);
                    super::mzxml::load_with_progress(path, &read, logger)?
                }
                None => super::mzxml::read_with_options(reader, &read)?,
            }
        }
        #[cfg(feature = "mzml")]
        FileType::MzData => {
            let limits = super::mzdata::ReadLimits::default();
            match file_logger.as_mut() {
                Some(logger) => {
                    drop(reader);
                    super::mzdata::load_with_progress(path, &options.peaks, &limits, logger)?
                }
                None => super::mzdata::read_with_options(reader, &options.peaks, &limits)?,
            }
            .experiment
        }
        FileType::Dta2d => {
            let read = dta2d_options(&options.peaks);
            match file_logger.as_mut() {
                Some(logger) => super::dta2d::read_with_progress(reader, &read, logger)?,
                None => super::dta2d::read_with_options(reader, &read)?,
            }
        }
        FileType::Dta | FileType::Mgf | FileType::Ms2 => {
            FileHandler::read_experiment(reader, kind)?
        }
        _ => return Err(unsupported(kind, "reading")),
    };
    if kind == FileType::MzMl {
        convert_srm_spectra(&mut result)?;
    }
    result.settings.document.loaded_file_path = document.loaded_file_path;
    result.settings.document.loaded_file_type = document.loaded_file_type;
    if options.rewrite_source_file {
        rewrite_source_file(&mut result, text, kind, options.compute_hash)?;
    }
    Ok(result)
}

/// Source `ChromatogramTools().convertSpectraToChromatograms<PeakMap>(exp,
/// true)`, the last step of every mzML load (`FileHandler.cpp:910`,
/// `KERNEL/ChromatogramTools.h:110-215`).
///
/// An experiment without SRM spectra is left alone, which is what the source
/// conversion does with it too, and costs one scan. Otherwise the conversion
/// runs with the default [`ChromatogramConversionLimits`] for its own output,
/// raised by the cost of scanning and moving the records already loaded, so
/// that a large experiment with a few SRM spectra is not refused for its size.
fn convert_srm_spectra(experiment: &mut MSExperiment) -> Result<()> {
    let is_srm = |spectrum: &&MSSpectrum| {
        spectrum.instrument_settings.scan_mode == ScanMode::SelectedReactionMonitoring
    };
    if !experiment.spectra.iter().any(|spectrum| is_srm(&spectrum)) {
        return Ok(());
    }
    // ChromatogramTools.h:162-166: every selected SRM spectrum without exactly
    // one precursor and at least one peak is skipped with this warning, in
    // experiment order, before any chromatogram is added.
    for spectrum in experiment.spectra.iter().filter(is_srm) {
        if spectrum.precursors.len() != 1 || spectrum.peaks.is_empty() {
            let _ = crate::concept::log_stream::log_message(
                crate::concept::log_stream::LogLevel::Warn,
                format_args!(
                    "ChromatogramTools: need exactly one precursor (given {}) and one or more product ions ({}), skipping conversion of this spectrum to chromatogram. If this is a MS1 chromatogram, please force conversion (e.g. with -convert_to_chromatograms).",
                    spectrum.precursors.len(),
                    spectrum.peaks.len()
                ),
            );
        }
    }
    let defaults = ChromatogramConversionLimits::default();
    let spectra = experiment.spectra.len();
    let chromatograms = experiment.chromatograms.len();
    let moved_bytes = spectra
        .saturating_mul(std::mem::size_of::<MSSpectrum>())
        .saturating_add(chromatograms.saturating_mul(std::mem::size_of::<MSChromatogram>()));
    let limits = ChromatogramConversionLimits {
        max_input_records: defaults.max_input_records.max(spectra),
        max_output_records: defaults.max_output_records.saturating_add(chromatograms),
        max_points: defaults.max_points,
        max_work: defaults
            .max_work
            .saturating_add(spectra.saturating_mul(3))
            .saturating_add(chromatograms),
        max_bytes: defaults.max_bytes.saturating_add(moved_bytes),
    };
    ChromatogramTools { limits }.convert_spectra_to_chromatograms(experiment, true, false)?;
    Ok(())
}

/// Source `loadExperiment`'s `rewrite_source_file` branch
/// (`FileHandler.cpp:1014-1051`); `file` is the name the caller passed.
fn rewrite_source_file(
    experiment: &mut MSExperiment,
    file: &str,
    kind: FileType,
    compute_hash: bool,
) -> Result<()> {
    let warn = |message: std::fmt::Arguments<'_>| {
        let _ = crate::concept::log_stream::log_message(
            crate::concept::log_stream::LogLevel::Warn,
            message,
        );
    };
    let sources = &experiment.settings.source_files;
    let mut source = match sources.first() {
        None => {
            warn(format_args!("No source file annotated."));
            SourceFile::default()
        }
        Some(first) => {
            if sources.len() > 1 {
                warn(format_args!(
                    "Expecting a single source file in mzML. Found {} will take only first one for rewriting.",
                    sources.len()
                ));
            }
            first.clone()
        }
    };
    let normalized = file.trim_end_matches(['/', '\\']);
    source.name = crate::system::file::basename(normalized).to_owned();
    let absolute = crate::system::file::absolute_path(normalized)?;
    let absolute = absolute
        .to_str()
        .ok_or_else(|| Error::InvalidValue("absolute path of the loaded file is not UTF-8".into()))?;
    // `std::filesystem::path::generic_string` writes `/` on every platform.
    let absolute = if cfg!(windows) {
        absolute.replace('\\', "/")
    } else {
        absolute.to_owned()
    };
    let directory = crate::system::file::path(&absolute);
    source.path = if directory.starts_with('/') {
        format!("file://{directory}")
    } else {
        format!("file:///{directory}")
    };
    source.file_type = kind.mzml_name().to_owned();
    #[cfg(feature = "mzml")]
    if compute_hash && kind != FileType::BrukerTdf {
        source.checksum = FileHandler::compute_file_hash(file);
        source.checksum_type = ChecksumType::Sha1;
    }
    #[cfg(not(feature = "mzml"))]
    let _ = (compute_hash, ChecksumType::Sha1);
    experiment.settings.source_files = vec![source];
    Ok(())
}

/// Name-based type, else bounded content recognition over replayed bytes, so
/// a nonseekable decompression stream is read once.
fn detect(path: &Path) -> Result<(FileType, Box<dyn BufRead>)> {
    let mut kind = type_by_file_name(filename(path)?);
    let mut reader = super::path_io::open(path)?;
    if kind == FileType::Unknown {
        let mut preview = Vec::new();
        reader
            .by_ref()
            .take(SNIFF_BYTES)
            .read_to_end(&mut preview)?;
        kind = type_by_content(&preview);
        reader = Box::new(BufReader::new(std::io::Cursor::new(preview).chain(reader)));
    }
    Ok((kind, reader))
}

fn typed_input(path: &Path, allowed: &[FileType], expected: FileType) -> Result<Box<dyn BufRead>> {
    let (kind, reader) = detect(path)?;
    check_allowed(kind, allowed)?;
    if kind != expected {
        return Err(Error::Unsupported(format!(
            "expected {} input, found {}",
            expected.name(),
            kind.name()
        )));
    }
    Ok(reader)
}

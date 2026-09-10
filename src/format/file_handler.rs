// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared experiment dispatch for TOPP-style pipelines. Format recognition is
//! separate from available adapters, whose documented representation limits apply.

use super::file_types::{consistent_output_type, type_by_file_name};
use super::{FileType, dta, mgf};
use crate::{Error, MSExperiment, Result};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

/// Native reader/writer dispatch; unsupported formats return explicit errors.
pub struct FileHandler;

impl FileHandler {
    pub fn can_read_experiment(kind: FileType) -> bool {
        matches!(
            kind,
            FileType::Dta | FileType::Mgf | FileType::Dta2d | FileType::Ms2
        ) || (kind == FileType::MzMl && cfg!(feature = "mzml"))
    }

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
                experiment.validate()?;
                if experiment.spectra.len() != 1
                    || !experiment.chromatograms.is_empty()
                    || !experiment.metadata.is_empty()
                {
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

    /// Load by extension, falling back to bounded content recognition only for
    /// unknown extensions. An empty allowed list accepts every available adapter.
    /// Gzip and bzip2 containers use the optional `file-compression` feature.
    pub fn load_experiment(path: impl AsRef<Path>, allowed: &[FileType]) -> Result<MSExperiment> {
        let path = path.as_ref();
        let kind = type_by_file_name(path.to_str().unwrap_or(""));
        load_stream(super::path_io::open(path)?, kind, allowed)
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
    pub fn load_feature_map(
        path: impl AsRef<Path>,
        allowed: &[FileType],
    ) -> Result<crate::kernel::FeatureMap> {
        let path = path.as_ref();
        let reader = map_input(path, allowed, FileType::FeatureXml)?;
        let mut map = Self::read_feature_map(reader, FileType::FeatureXml)?;
        map.loaded_file_path = filename(path)?.into();
        map.loaded_file_type = FileType::FeatureXml;
        Ok(map)
    }
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
    pub fn load_consensus_map(
        path: impl AsRef<Path>,
        allowed: &[FileType],
    ) -> Result<crate::kernel::ConsensusMap> {
        let path = path.as_ref();
        let reader = map_input(path, allowed, FileType::ConsensusXml)?;
        let mut map = Self::read_consensus_map(reader, FileType::ConsensusXml)?;
        map.loaded_file_path = filename(path)?.into();
        map.loaded_file_type = FileType::ConsensusXml;
        Ok(map)
    }
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

fn load_stream(
    mut reader: impl BufRead,
    kind: FileType,
    allowed: &[FileType],
) -> Result<MSExperiment> {
    if kind != FileType::Unknown {
        check_allowed(kind, allowed)?;
        return FileHandler::read_experiment(reader, kind);
    }
    let mut preview = Vec::new();
    // The same bytes are replayed to the parser. This also works with nonseekable
    // decompression streams and caps both sniffing memory and decompression work.
    reader.by_ref().take(65_536).read_to_end(&mut preview)?;
    let kind = type_by_content(&preview);
    check_allowed(kind, allowed)?;
    FileHandler::read_experiment(BufReader::new(preview.as_slice().chain(reader)), kind)
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
    FileType::Unknown
}

fn filename(path: &Path) -> Result<&str> {
    path.to_str()
        .ok_or_else(|| Error::InvalidValue("filename must be UTF-8".into()))
}
fn map_input(path: &Path, allowed: &[FileType], expected: FileType) -> Result<Box<dyn BufRead>> {
    let mut kind = type_by_file_name(filename(path)?);
    let mut reader = super::path_io::open(path)?;
    if kind == FileType::Unknown {
        let mut preview = Vec::new();
        reader.by_ref().take(65_536).read_to_end(&mut preview)?;
        kind = type_by_content(&preview);
        reader = Box::new(BufReader::new(std::io::Cursor::new(preview).chain(reader)));
    }
    check_allowed(kind, allowed)?;
    if kind != expected {
        return Err(Error::Unsupported(format!(
            "expected {} map, found {}",
            expected.name(),
            kind.name()
        )));
    }
    Ok(reader)
}

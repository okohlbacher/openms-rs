// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared experiment dispatch for TOPP-style pipelines. Format recognition is
//! separate from available adapters, whose documented representation limits apply.

use super::file_types::{consistent_output_type, type_by_file_name};
use super::{FileType, dta, mgf};
use crate::{Error, MSExperiment, Result};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

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
    /// Gzip containers use the existing `mzml` feature's compression dependency.
    pub fn load_experiment(path: impl AsRef<Path>, allowed: &[FileType]) -> Result<MSExperiment> {
        let path = path.as_ref();
        let kind = type_by_file_name(path.to_str().unwrap_or(""));
        let mut reader = BufReader::new(File::open(path)?);
        let magic = reader.fill_buf()?;
        if magic.starts_with(&[0x1f, 0x8b]) {
            #[cfg(feature = "mzml")]
            {
                return load_stream(
                    BufReader::new(flate2::read::MultiGzDecoder::new(reader)),
                    kind,
                    allowed,
                );
            }
            #[cfg(not(feature = "mzml"))]
            {
                return Err(Error::Unsupported(
                    "gzip input requires the mzml feature".into(),
                ));
            }
        }
        if magic.starts_with(b"BZh") || magic.starts_with(b"PK\x03\x04") {
            return Err(Error::Unsupported("bzip2/ZIP experiment input".into()));
        }
        load_stream(reader, kind, allowed)
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
        let suffix = path.extension().and_then(|s| s.to_str()).unwrap_or("");
        if ["bz2", "zip"]
            .iter()
            .any(|s| suffix.eq_ignore_ascii_case(s))
        {
            return Err(Error::Unsupported("bzip2/ZIP experiment output".into()));
        }
        let gzip = suffix.eq_ignore_ascii_case("gz");
        if gzip && !cfg!(feature = "mzml") {
            return Err(Error::Unsupported(
                "gzip output requires the mzml feature".into(),
            ));
        }
        let (temporary, mut file) = TemporaryFile::create(path)?;
        {
            let mut writer = BufWriter::new(&mut file);
            if gzip {
                #[cfg(feature = "mzml")]
                {
                    let mut encoder =
                        flate2::write::GzEncoder::new(&mut writer, flate2::Compression::default());
                    Self::write_experiment(&mut encoder, experiment, kind)?;
                    encoder.finish()?;
                }
            } else {
                Self::write_experiment(&mut writer, experiment, kind)?;
            }
            writer.flush()?;
        }
        file.sync_all()?;
        drop(file);
        std::fs::rename(&temporary.0, path)?;
        Ok(())
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

struct TemporaryFile(PathBuf);
impl TemporaryFile {
    fn create(destination: &Path) -> Result<(Self, File)> {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self::create_with_counter(destination, &NEXT)
    }

    fn create_with_counter(destination: &Path, next: &AtomicU64) -> Result<(Self, File)> {
        let parent = destination.parent().unwrap_or_else(|| Path::new("."));
        for _ in 0..32 {
            let id = next.fetch_add(1, Ordering::Relaxed);
            let name = format!(".openms-{}-{id}.tmp", std::process::id());
            // Renaming a temporary file onto itself would succeed, after which
            // its cleanup guard would delete the output. Compare basenames so
            // relative paths and case-insensitive filesystems are safe too.
            if destination
                .file_name()
                .and_then(|s| s.to_str())
                .is_some_and(|s| s.eq_ignore_ascii_case(&name))
            {
                continue;
            }
            let path = parent.join(name);
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(file) => return Ok((Self(path), file)),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.into()),
            }
        }
        Err(Error::Io(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "cannot allocate temporary output",
        )))
    }
}
impl Drop for TemporaryFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_named_like_temporary_survives_cleanup() {
        let directory = std::env::temp_dir().join(format!(
            "openms-handler-temp-collision-{}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        for uppercase in [false, true] {
            let mut name = format!(".openms-{}-0.tmp", std::process::id());
            if uppercase {
                name.make_ascii_uppercase();
            }
            let destination = directory.join(".").join(name);
            let next = AtomicU64::new(0);
            let (temporary, mut file) =
                TemporaryFile::create_with_counter(&destination, &next).unwrap();
            assert_eq!(next.load(Ordering::Relaxed), 2);
            file.write_all(b"complete output").unwrap();
            drop(file);
            std::fs::rename(&temporary.0, &destination).unwrap();
            drop(temporary);
            assert_eq!(std::fs::read(&destination).unwrap(), b"complete output");
            std::fs::remove_file(destination).unwrap();
        }
        assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 0);
        std::fs::remove_dir(directory).unwrap();
    }
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Complete pinned SDK file-type registry and lexical filename helpers.
//! Registry properties describe the C++ SDK. Use `FileHandler` capability queries
//! to determine which native Rust adapters are available in this build.

use crate::{Error, Result};

/// Capabilities advertised by the pinned C++ SDK, independently of Rust support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FileProperty {
    Readable,
    Writeable,
    ProvidesSpectrum,
    ProvidesExperiment,
    ProvidesFeatures,
    ProvidesConsensusFeatures,
    ProvidesIdentifications,
    ProvidesTransitions,
    ProvidesQuantifications,
    ProvidesTransformations,
    ProvidesQc,
}

// One table keeps the enum, names, descriptions and property lists in agreement.
macro_rules! file_types {
    ($($variant:ident => ($name:literal, $description:literal, [$($property:ident),*])),* $(,)?) => {
        /// Recognized SDK formats. Recognition does not imply an available reader.
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum FileType { $($variant),* }

        impl FileType {
            /// All formats in source annotation order (XML deliberately last).
            pub const ALL: &'static [Self] = &[$(Self::$variant),*];
            /// Preferred source extension, without a leading dot.
            pub fn name(self) -> &'static str { match self { $(Self::$variant => $name),* } }
            /// Human-readable source format description.
            pub fn description(self) -> &'static str { match self { $(Self::$variant => $description),* } }
            /// Source capabilities, not native reader/writer availability.
            pub fn source_properties(self) -> &'static [FileProperty] {
                match self { $(Self::$variant => &[$(FileProperty::$property),*]),* }
            }
        }
    }
}

// Source table: Core SDK 6bfc0e4, FORMAT/FileTypes.cpp.
file_types! {
    Unknown => ("unknown", "unknown file extension", []),
    Dta => ("dta", "dta raw data file", [ProvidesExperiment, ProvidesSpectrum, Readable, Writeable]),
    Dta2d => ("dta2d", "dta2d raw data file", [ProvidesExperiment, Readable, Writeable]),
    MzData => ("mzData", "mzData raw data file", [ProvidesExperiment, Readable, Writeable]),
    MzXml => ("mzXML", "mzXML raw data file", [ProvidesExperiment, Readable, Writeable]),
    FeatureXml => ("featureXML", "OpenMS feature map", [ProvidesFeatures, Readable, Writeable]),
    IdXml => ("idXML", "OpenMS peptide identification file", [ProvidesIdentifications, Readable, Writeable]),
    ConsensusXml => ("consensusXML", "OpenMS consensus feature map", [ProvidesConsensusFeatures, Readable, Writeable]),
    Mgf => ("mgf", "mascot generic format file", [ProvidesExperiment, Readable, Writeable]),
    Ini => ("ini", "OpenMS parameter file", [Readable]),
    Toppas => ("toppas", "OpenMS TOPPAS pipeline", [Readable]),
    TransformationXml => ("trafoXML", "RT transformation file", [ProvidesTransformations, Readable, Writeable]),
    MzMl => ("mzML", "mzML raw data file", [ProvidesExperiment, Readable, Writeable]),
    CachedMzMl => ("cachedMzML", "cachedMzML raw data file", [Readable, Writeable]),
    Ms2 => ("ms2", "ms2 file", [ProvidesExperiment, Readable]),
    PepXml => ("pepXML", "pepXML file", [Readable, Writeable]),
    ProtXml => ("protXML", "protXML file", [ProvidesIdentifications, Readable]),
    MzIdentMl => ("mzid", "mzIdentML file", [ProvidesIdentifications, Readable, Writeable]),
    QcMl => ("qcml", "quality control file", [ProvidesQc, Writeable]),
    MzQc => ("mzqc", "quality control file in json format", [ProvidesQc, Writeable]),
    GelMl => ("gelML", "gelML file", []),
    TraMl => ("traML", "transition file", [ProvidesTransitions, Readable, Writeable]),
    Msp => ("msp", "NIST spectra library file format", [ProvidesExperiment, Readable, Writeable]),
    OmssaXml => ("omssaXML", "omssaXML file", [ProvidesIdentifications, Readable]),
    MascotXml => ("mascotXML", "mascotXML file", []),
    Png => ("png", "portable network graphics file", []),
    Xmass => ("fid", "XMass analysis file", [ProvidesExperiment, ProvidesSpectrum, Readable, Writeable]),
    Tsv => ("tsv", "tab-separated file", [ProvidesFeatures, Readable, Writeable]),
    MzTab => ("mzTab", "mzTab file", []),
    Peplist => ("peplist", "SpecArray file", [ProvidesFeatures, Readable, Writeable]),
    Hardkloer => ("hardkloer", "hardkloer file", []),
    Kroenik => ("kroenik", "kroenik file", [ProvidesFeatures, Readable, Writeable]),
    Fasta => ("fasta", "FASTA file", [Readable, Writeable]),
    Peff => ("peff", "PEFF protein file", [Readable, Writeable]),
    Edta => ("edta", "enhanced dta file", [ProvidesFeatures, ProvidesConsensusFeatures, Readable, Writeable]),
    Csv => ("csv", "comma-separated values file", [Readable, Writeable]),
    Txt => ("txt", "generic text file", []),
    Obo => ("obo", "controlled vocabulary file", []),
    Html => ("html", "any HTML file", []),
    AnalysisXml => ("analysisXML", "analysisXML file", []),
    Xsd => ("xsd", "XSD schema format", []),
    Psq => ("psq", "NCBI binary blast db", []),
    Mrm => ("mrm", "SpectraST MRM list", [Readable]),
    SqMass => ("sqMass", "SQLite format for mass and chromatograms", [Readable, Writeable]),
    Pqp => ("pqp", "pqp file", [Readable, Writeable]),
    Oswpq => ("oswpq", "OpenSwath Parquet bundle", [Readable, Writeable]),
    Ms => ("ms", "SIRIUS file", []),
    Osw => ("osw", "OpenSwath output files", [Readable, Writeable]),
    ChromParquet => ("xic", "OpenSwath Parquet chromatogram output", [Readable, Writeable]),
    MobilParquet => ("xim", "OpenSwath Parquet mobilogram output", [Readable, Writeable]),
    PeakMapParquet => ("xipm", "OpenSwath Parquet peak-map output", [Readable, Writeable]),
    Psms => ("psms", "Percolator tab-delimited output (PSM level)", [Readable]),
    Pin => ("pin", "Percolator tab-delimited input (PSM level)", []),
    ParamXml => ("paramXML", "OpenMS internal XML file", []),
    Splib => ("splib", "SpectraST binary spectral library file", []),
    Novor => ("novor", "Novor custom parameter file", []),
    XquestXml => ("xquest.xml", "xquest.xml file", [ProvidesIdentifications, Readable, Writeable]),
    SpecXml => ("spec.xml", "spec.xml file", []),
    Json => ("json", "JavaScript Object Notation file", [Readable, Writeable]),
    Raw => ("raw", "(Thermo) Raw data file", [ProvidesExperiment, Readable]),
    Oms => ("oms", "OpenMS SQLite file", [ProvidesIdentifications, ProvidesFeatures, ProvidesConsensusFeatures]),
    Exe => ("exe", "Windows executable", []),
    Bz2 => ("bz2", "bzip2 compressed file", [Readable]),
    Gz => ("gz", "gzip compressed file", [Readable]),
    Zip => ("zip", "ZIP compressed file", [Readable]),
    Parquet => ("parquet", "Apache Parquet file", [Readable, Writeable]),
    IdParquet => ("idparquet", "OpenMS identification parquet bundle (directory)", [ProvidesIdentifications, Readable, Writeable]),
    FeatureParquet => ("featureparquet", "OpenMS feature map parquet bundle (directory)", [ProvidesFeatures, ProvidesIdentifications, Readable, Writeable]),
    ConsensusParquet => ("consensusparquet", "OpenMS consensus map parquet bundle (directory)", [ProvidesConsensusFeatures, ProvidesIdentifications, Readable, Writeable]),
    BrukerTdf => ("d", "Bruker TDF", [ProvidesExperiment, Readable]),
    ImzMl => ("imzML", "imzML mass spectrometry imaging file", [ProvidesExperiment, Readable, Writeable]),
    Yaml => ("yaml", "YAML file", [Writeable]),
    Xml => ("xml", "any XML file", [Readable]),
}

impl FileType {
    /// Case-insensitive source extension lookup, including the `pqt` alias.
    pub fn from_name(name: &str) -> Self {
        if name.eq_ignore_ascii_case("pqt") {
            return Self::Parquet;
        }
        Self::ALL
            .iter()
            .copied()
            .find(|t| t.name().eq_ignore_ascii_case(name))
            .unwrap_or(Self::Unknown)
    }

    /// Whether the format represents a directory rather than an ordinary file.
    pub fn is_directory(self) -> bool {
        matches!(
            self,
            Self::BrukerTdf | Self::IdParquet | Self::FeatureParquet | Self::ConsensusParquet
        )
    }

    /// mzML source-file CV name used by the pinned SDK, or an empty string.
    pub fn mzml_name(self) -> &'static str {
        match self {
            Self::Dta | Self::Dta2d => "DTA file",
            Self::MzMl | Self::ImzMl => "mzML file",
            Self::MzData => "PSI mzData file",
            Self::MzXml => "ISB mzXML file",
            Self::Mgf => "Mascot MGF file",
            Self::Xmass => "Bruker FID file",
            Self::BrukerTdf => "Bruker TDF format",
            Self::Raw => "Thermo RAW format",
            _ => "",
        }
    }
}

/// Layout of a source-compatible file selection filter string.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilterLayout {
    Compact,
    OneByOne,
    Both,
}

/// Ordered file types; repeated entries and empty lists retain source behavior.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FileTypeList(pub Vec<FileType>);

impl FileTypeList {
    pub fn contains(&self, value: FileType) -> bool {
        self.0.contains(&value)
    }

    /// Select source formats having every requested property, in source order.
    pub fn with_source_properties(properties: &[FileProperty]) -> Self {
        Self(
            FileType::ALL
                .iter()
                .copied()
                .filter(|t| properties.iter().all(|p| t.source_properties().contains(p)))
                .collect(),
        )
    }

    fn filter_elements(&self, layout: FilterLayout, add_all: bool) -> Vec<(String, FileType)> {
        let mut result = Vec::new();
        if matches!(layout, FilterLayout::Compact | FilterLayout::Both) {
            let extensions = self
                .0
                .iter()
                .map(|t| format!("*.{}", t.name()))
                .collect::<Vec<_>>()
                .join(" ");
            result.push((
                format!("all readable files ({extensions})"),
                FileType::Unknown,
            ));
        }
        if matches!(layout, FilterLayout::OneByOne | FilterLayout::Both) {
            result.extend(
                self.0
                    .iter()
                    .map(|t| (format!("{} (*.{})", t.description(), t.name()), *t)),
            );
        }
        if add_all {
            result.push(("all files (*)".into(), FileType::Unknown));
        }
        result
    }

    pub fn to_file_dialog_filter(&self, layout: FilterLayout, add_all: bool) -> String {
        self.filter_elements(layout, add_all)
            .into_iter()
            .map(|(s, _)| s)
            .collect::<Vec<_>>()
            .join(";;")
    }

    /// Resolve an exact filter item; use `fallback` for aggregate/all-file items.
    pub fn from_file_dialog_filter(&self, filter: &str, fallback: FileType) -> Result<FileType> {
        self.filter_elements(FilterLayout::Both, true)
            .into_iter()
            .find(|(s, _)| s == filter)
            .map(|(_, t)| if t == FileType::Unknown { fallback } else { t })
            .ok_or_else(|| Error::InvalidValue(format!("unknown file dialog filter: {filter}")))
    }
}

/// Lexical filename lookup, without opening input or output paths.
/// Compression suffixes are peeled iteratively. Both path separators are valid.
pub fn type_by_file_name(filename: &str) -> FileType {
    let mut base = filename.rsplit(['/', '\\']).next().unwrap_or("");
    loop {
        for (suffix, kind) in [
            (".pep.xml", FileType::PepXml),
            (".prot.xml", FileType::ProtXml),
            (".xquest.xml", FileType::XquestXml),
            (".spec.xml", FileType::SpecXml),
        ] {
            if base.ends_with(suffix) {
                return kind;
            }
        }
        let Some((stem, extension)) = base.rsplit_once('.') else {
            return if base == "fid" {
                FileType::Xmass
            } else {
                FileType::Unknown
            };
        };
        if ["bz2", "gz", "zip"]
            .iter()
            .any(|s| extension.eq_ignore_ascii_case(s))
        {
            base = stem;
        } else {
            return FileType::from_name(extension);
        }
    }
}

/// An unknown extension is permitted by the source's output-extension check.
pub fn has_valid_extension(filename: &str, expected: FileType) -> bool {
    let actual = type_by_file_name(filename);
    actual == expected || actual == FileType::Unknown
}

/// Source extension stripping, including its lexical handling of alias suffixes.
/// For example, `sample.pep.xml` becomes `sample.pep`, as in the source tests.
pub fn strip_extension(filename: &str) -> &str {
    let Some(last_dot) = filename.rfind('.') else {
        return filename;
    };
    if let Some(pos) = filename
        .to_ascii_lowercase()
        .rfind(&type_by_file_name(filename).name().to_ascii_lowercase())
    {
        // Source underflows if the name starts with the type without a dot.
        if let Some(end) = pos.checked_sub(1).filter(|&p| filename.is_char_boundary(p)) {
            return &filename[..end];
        }
        return filename;
    }
    if filename.rfind(['/', '\\']).is_some_and(|p| p > last_dot) {
        filename
    } else {
        &filename[..last_dot]
    }
}

pub fn swap_extension(filename: &str, kind: FileType) -> String {
    format!("{}.{}", strip_extension(filename), kind.name())
}

/// Resolve an output suffix and explicit type using the source conflict rules.
pub fn consistent_output_type(filename: &str, requested_type: &str) -> FileType {
    let file = type_by_file_name(filename);
    let requested = FileType::from_name(requested_type);
    if file == FileType::Unknown {
        requested
    } else if requested == FileType::Unknown || file == requested {
        file
    } else {
        FileType::Unknown
    }
}

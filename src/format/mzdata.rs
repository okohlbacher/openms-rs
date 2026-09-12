// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded mzData 1.05 I/O: `FORMAT/MzDataFile.h` and its
//! `FORMAT/HANDLERS/MzDataHandler.h`. See `docs/MZDATA_SUPPORT.md` for the
//! member-by-member mapping and the metadata mzData cannot represent.
//!
//! mzData is the PSI peak format that preceded mzXML and mzML. A document is a
//! `<description>` block of experiment metadata followed by a `<spectrumList>`
//! of `<spectrum>` elements, each carrying a `<spectrumDesc>` and one
//! `<mzArrayBinary>` / `<intenArrayBinary>` pair plus any number of
//! `<supDataArrayBinary>` annotation arrays. Metadata is expressed as
//! `<cvParam accession="PSI:1000xxx" value="…"/>` against a small fixed
//! vocabulary that the handler hard-codes rather than loading from an OBO file,
//! and as free `<userParam name="…" value="…"/>` pairs.
//!
//! Unlike mzML, each binary array declares its own `precision` (`"32"` or
//! `"64"`), `endian` (`"little"` or `"big"`) and `length` as XML attributes
//! instead of CV terms, and the source honours all three on reading. Both byte
//! orders are supported here for every array, including the m/z and intensity
//! arrays; see
//! [`Endian`](crate::format::mzdata::Endian) and
//! [`Precision`](crate::format::mzdata::Precision). The writer emits little
//! endian only, as `MzDataHandler::writeBinary_` does.
//!
//! Entry points are
//! [`MzDataFile`](crate::format::mzdata::MzDataFile) for the source file
//! adapter, [`load`](crate::format::mzdata::load) /
//! [`store`](crate::format::mzdata::store) for paths and
//! [`read`](crate::format::mzdata::read) /
//! [`write`](crate::format::mzdata::write) for streams. Every declared array
//! length is attacker-controlled, so reads are preflighted against the
//! explicit ceilings in
//! [`ReadLimits`](crate::format::mzdata::ReadLimits) before anything is
//! allocated, and the parsed experiment is committed only once the whole
//! document has been read.
//!
//! The source parallelises nothing here, and neither does this port; the
//! `ProgressLogger` the source handler takes by reference is replaced by the
//! counters in [`LoadReport`](crate::format::mzdata::LoadReport).

use crate::format::peak_options::PeakFileOptions;
use crate::kernel::{
    DataArray, MSExperiment, MSSpectrum, NumericRange, Peak1D, Precursor, SpectrumType,
};
use crate::metadata::{
    Acquisition, ActivationMethod, AnalyzerType, ContactPerson, DataProcessing,
    DetectorAcquisitionMode, DetectorType, InletType, IonizationMethod, MetaInfo, MetaValue,
    MetaValueData, Polarity, ProcessingAction, ReflectronState, ResolutionMethod, ResolutionType,
    SampleState, ScanDirection, ScanLaw, ScanMode, ScanWindow, SourceFile,
};
use crate::{Error, Result};
use base64::{Engine, engine::general_purpose::STANDARD};
use quick_xml::{
    Reader,
    events::{BytesEnd, BytesStart, Event},
};
use std::collections::BTreeMap;
use std::io::{BufRead, Read, Write};
use std::path::Path;
use std::sync::Arc;

/// Schema version the source `MzDataFile` constructor pins and the writer
/// emits as the root element's `version` attribute (`MzDataFile.cpp:18`).
pub const SCHEMA_VERSION: &str = "1.05";

/// Schema resource the source `MzDataFile` constructor names
/// (`MzDataFile.cpp:18`). No schema validation is implemented here; see
/// [`MzDataFile::is_semantically_valid`].
pub const SCHEMA_FILE: &str = "/SCHEMAS/mzData_1_05.xsd";

/// `xsi:noNamespaceSchemaLocation` the writer emits, verbatim from
/// `MzDataHandler.cpp:584`. mzData has no XML namespace, so nothing in a
/// document is namespace-qualified except this attribute.
pub const SCHEMA_LOCATION: &str = "http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd";

// ---------------------------------------------------------------------------
// The hard-coded controlled vocabulary
// ---------------------------------------------------------------------------
//
// `MzDataHandler::init_` (MzDataHandler.cpp:42-83) fills 19 `cv_terms_` tables
// by splitting semicolon-separated literals. The leading empty entry in most
// of them is the enum's zero value, so the index of a term in its table is the
// numeric value of the corresponding OpenMS enum; `cvStringToEnum_`
// (XMLHandler.cpp:131-145) returns that index and falls back to 0 with a
// warning. Tables 4, 12, 15, 16 and 17 are commented "no longer used" and are
// empty; table 18 is the one with no leading empty entry, so an unknown
// activation method resolves to `CID` rather than to an unknown value.

/// `cv_terms_[0]`, sample state.
const SAMPLE_STATE: &[&str] = &[
    "",
    "Solid",
    "Liquid",
    "Gas",
    "Solution",
    "Emulsion",
    "Suspension",
];
/// `cv_terms_[1]`, ionization mode, which OpenMS stores as a polarity.
const IONIZATION_MODE: &[&str] = &["", "PositiveIonMode", "NegativeIonMode"];
/// `cv_terms_[2]`, resolution method.
const RESOLUTION_METHOD: &[&str] = &["", "FWHM", "TenPercentValley", "Baseline"];
/// `cv_terms_[3]`, resolution type.
const RESOLUTION_TYPE: &[&str] = &["", "Constant", "Proportional"];
/// `cv_terms_[5]`, scan direction.
const SCAN_DIRECTION: &[&str] = &["", "Up", "Down"];
/// `cv_terms_[6]`, scan law.
const SCAN_LAW: &[&str] = &["", "Exponential", "Linear", "Quadratic"];
/// `cv_terms_[8]`, reflectron state.
const REFLECTRON_STATE: &[&str] = &["", "On", "Off", "None"];
/// `cv_terms_[9]`, detector acquisition mode.
const ACQUISITION_MODE: &[&str] = &["", "PulseCounting", "ADC", "TDC", "TransientRecorder"];
/// `cv_terms_[10]`, ionization type.
const IONIZATION_TYPE: &[&str] = &[
    "", "ESI", "EI", "CI", "FAB", "TSP", "LD", "FD", "FI", "PD", "SI", "TI", "API", "ISI", "CID",
    "CAD", "HN", "APCI", "APPI", "ICP",
];
/// `cv_terms_[11]`, inlet type.
const INLET_TYPE: &[&str] = &[
    "",
    "Direct",
    "Batch",
    "Chromatography",
    "ParticleBeam",
    "MembraneSeparator",
    "OpenSplit",
    "JetSeparator",
    "Septum",
    "Reservoir",
    "MovingBelt",
    "MovingWire",
    "FlowInjectionAnalysis",
    "ElectrosprayInlet",
    "ThermosprayInlet",
    "Infusion",
    "ContinuousFlowFastAtomBombardment",
    "InductivelyCoupledPlasma",
];
/// `cv_terms_[13]`, detector type.
const DETECTOR_TYPE: &[&str] = &[
    "",
    "EM",
    "Photomultiplier",
    "FocalPlaneArray",
    "FaradayCup",
    "ConversionDynodeElectronMultiplier",
    "ConversionDynodePhotomultiplier",
    "Multi-Collector",
    "ChannelElectronMultiplier",
];
/// `cv_terms_[14]`, mass analyzer type.
const ANALYZER_TYPE: &[&str] = &[
    "",
    "Quadrupole",
    "PaulIonTrap",
    "RadialEjectionLinearIonTrap",
    "AxialEjectionLinearIonTrap",
    "TOF",
    "Sector",
    "FourierTransform",
    "IonStorage",
];
/// `cv_terms_[18]`, activation method. The only table without a leading empty
/// entry, so index 0 is a real value.
const ACTIVATION_METHOD: &[&str] = &["CID", "PSD", "PD", "SID"];

// ---------------------------------------------------------------------------
// Binary array transport
// ---------------------------------------------------------------------------

/// Byte order a binary array declares in its `endian` attribute.
///
/// mzData states the byte order per array rather than fixing it for the format,
/// and `MzDataHandler::fillData_` (`MzDataHandler.cpp:504-520`) dispatches on
/// it for every array independently, so one spectrum may mix both orders.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Endian {
    /// `endian="little"`, and every value the source does not spell `"big"`.
    #[default]
    Little,
    /// `endian="big"`.
    Big,
}

impl Endian {
    /// The attribute spelling, as the writer would emit it.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Little => "little",
            Self::Big => "big",
        }
    }

    /// Classify an `endian` attribute exactly as the source does.
    ///
    /// The source compares against `"big"` and treats every other spelling —
    /// including a misspelling and the empty string — as little endian
    /// (`MzDataHandler.cpp:506`, `:523`). `Ok(false)` in the second tuple
    /// position marks a spelling that is neither `"big"` nor `"little"`, which
    /// the caller records as a warning; the source is silent about it.
    pub fn parse(value: &str) -> (Self, bool) {
        match value {
            "big" => (Self::Big, true),
            "little" => (Self::Little, true),
            _ => (Self::Little, false),
        }
    }
}

/// Element width a binary array declares in its `precision` attribute.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Precision {
    /// `precision="32"`: IEEE-754 binary32 elements.
    Bits32,
    /// `precision="64"`: IEEE-754 binary64 elements, and every value the
    /// source does not spell `"32"`.
    #[default]
    Bits64,
}

impl Precision {
    /// The attribute spelling, as the writer would emit it.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Bits32 => "32",
            Self::Bits64 => "64",
        }
    }

    /// Bytes per element.
    pub const fn width(self) -> usize {
        match self {
            Self::Bits32 => 4,
            Self::Bits64 => 8,
        }
    }

    /// Classify a `precision` attribute exactly as the source does.
    ///
    /// The source tests only for `"32"` and treats every other spelling as
    /// 64-bit (`MzDataHandler.cpp:503`, `:522`, `:527`). `Ok(false)` in the
    /// second tuple position marks a spelling that is neither `"32"` nor
    /// `"64"`, which the caller records as a warning.
    pub fn parse(value: &str) -> (Self, bool) {
        match value {
            "32" => (Self::Bits32, true),
            "64" => (Self::Bits64, true),
            _ => (Self::Bits64, false),
        }
    }
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

/// Resource ceilings applied while loading an mzData document.
///
/// Every one of these bounds something an attacker controls through the file:
/// the `length` attribute of a binary array, the `count` attribute of
/// `<spectrumList>`, the number of elements a base64 payload decodes to and
/// the total text the document carries. Each is checked *before* the
/// allocation it bounds, so a refused document never commits the memory it
/// asked for.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReadLimits {
    /// Maximum input bytes, counted on the document as it lies on disk and
    /// before any ISO-8859-1 transcoding.
    pub max_xml_bytes: u64,
    /// Maximum nesting depth of open elements.
    pub max_xml_depth: usize,
    /// Maximum retained spectra, and the ceiling the `<spectrumList count>`
    /// attribute is checked against.
    pub max_spectra: usize,
    /// Maximum `<data>` elements in one spectrum, counting the m/z and
    /// intensity arrays.
    pub max_arrays_per_spectrum: usize,
    /// Maximum decoded bytes for any one binary array.
    pub max_array_bytes: usize,
    /// Maximum elements in any one binary array, and the ceiling a declared
    /// `length` attribute is checked against.
    pub max_array_elements: usize,
    /// Maximum retained peaks across the whole document, counted after the
    /// `PeakFileOptions` filters.
    pub max_total_peaks: usize,
    /// Cumulative character-data bytes, including every base64 payload.
    pub max_text_bytes: usize,
    /// Maximum metadata entries the document may create across all records.
    pub max_metadata_entries: usize,
    /// Maximum warnings retained in [`LoadReport::warnings`];
    /// [`LoadReport::warning_count`] keeps counting past it.
    pub max_warnings: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_xml_bytes: 512 * 1024 * 1024,
            max_xml_depth: 64,
            max_spectra: 1_000_000,
            max_arrays_per_spectrum: 4096,
            max_array_bytes: 64 * 1024 * 1024,
            max_array_elements: 20_000_000,
            max_total_peaks: 50_000_000,
            max_text_bytes: 512 * 1024 * 1024,
            max_metadata_entries: 1_000_000,
            max_warnings: 256,
        }
    }
}

/// What a load could not represent faithfully, and what the filters removed.
///
/// The source reports all of this through `XMLHandler::warning`, which
/// `#ifdef`s itself down to a debug-level log in a release build
/// (`XMLHandler.cpp:87-107`), so in practice none of it reaches the caller.
/// Returning it is what lets a caller notice that a document declared a
/// `length` that disagreed with its payload, or that a `cvParam` named a term
/// the handler's fixed vocabulary does not contain.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LoadReport {
    /// Warnings in document order, at most [`ReadLimits::max_warnings`] of them.
    pub warnings: Vec<String>,
    /// How many warnings were raised in total.
    pub warning_count: usize,
    /// Spectra dropped whole by the MS-level, retention-time or precursor-m/z
    /// filters, the source `skip_spectrum_`.
    pub spectra_skipped: usize,
    /// Peaks dropped by the m/z and intensity range filters.
    pub peaks_filtered_out: usize,
}

impl LoadReport {
    /// Whether the document loaded with nothing dropped and nothing warned
    /// about.
    ///
    /// Native convenience; the source has no report to summarise.
    pub fn is_clean(&self) -> bool {
        self.warning_count == 0 && self.spectra_skipped == 0 && self.peaks_filtered_out == 0
    }
}

/// An experiment together with what its load dropped or warned about.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Loaded {
    /// The parsed experiment.
    pub experiment: MSExperiment,
    /// Warnings and filter counts.
    pub report: LoadReport,
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// Choices the mzData writer makes, and whether it may discard.
///
/// mzData stores far less than an [`MSExperiment`] holds: one source file, one
/// ion source, one ion detector, one scan window per scan, one data-processing
/// record for the whole document, no chromatograms and no integer or string
/// data arrays. `MzDataHandler::writeTo` discards all of it, warning about
/// four cases and silent about the rest. The default here refuses instead;
/// [`WriteOptions::source`] selects the source behaviour, following the
/// established `dta::WriteOptions::source` pattern.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WriteOptions {
    /// Write `<supDesc>` and `<supDataArrayBinary>` for each float data array,
    /// the source `PeakFileOptions::getWriteSupplementalData()`
    /// (`MzDataHandler.cpp:983`, `:1023`).
    pub write_supplemental_data: bool,
    /// Write the m/z array as `precision="32"`, as `writeBinary_` hardcodes
    /// (`MzDataHandler.cpp:1500-1506`).
    ///
    /// Default `false`: an f64 m/z is written as `precision="64"`, which mzData
    /// declares per array and which both this reader and the source reader
    /// honour. The source narrows every coordinate through a
    /// `std::vector<float>` and so loses about nine significant digits of every
    /// mass.
    pub mz_32_bit: bool,
    /// Discard what mzData cannot represent instead of refusing to write it.
    pub discard_unrepresentable: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            write_supplemental_data: true,
            mz_32_bit: false,
            discard_unrepresentable: false,
        }
    }
}

impl WriteOptions {
    /// The `MzDataHandler::writeTo` behaviour: 32-bit coordinates, and discard
    /// every record mzData has no element for.
    pub fn source() -> Self {
        Self {
            write_supplemental_data: true,
            mz_32_bit: true,
            discard_unrepresentable: true,
        }
    }
}

/// What a store discarded, and what it had to renumber or substitute.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct StoreReport {
    /// Warnings in document order. Unbounded work is impossible here: the
    /// writer raises at most a fixed number of warnings per spectrum.
    pub warnings: Vec<String>,
    /// How many warnings were raised in total.
    pub warning_count: usize,
    /// True when the native IDs could not all be expressed as mzData's numeric
    /// `id` attribute and the spectra were renumbered from 1
    /// (`MzDataHandler.cpp:790-793`).
    pub renumbered: bool,
}

// ---------------------------------------------------------------------------
// The file adapter
// ---------------------------------------------------------------------------

/// mzData file adapter, the source `MzDataFile`.
///
/// The source class is `Internal::XMLFile` plus `ProgressLogger` plus a
/// `PeakFileOptions` member, and its three real methods are `load`, `store`
/// and `isSemanticallyValid`. This struct keeps the options member and adds the
/// native read ceilings and the write policy, because neither has a source
/// counterpart to store them.
///
/// ```
/// # #[cfg(feature = "mzml")] {
/// use openms::format::mzdata::MzDataFile;
/// let mut file = MzDataFile::new();
/// assert!(!file.options().has_ms_levels());
/// file.options_mut().add_ms_level(1).unwrap();
/// assert!(file.options().has_ms_levels());
/// assert_eq!(file.version(), "1.05");
/// # }
/// ```
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MzDataFile {
    options: PeakFileOptions,
    limits: ReadLimits,
    discard_unrepresentable: bool,
}

impl MzDataFile {
    /// A file adapter with default options, default ceilings and a writer that
    /// refuses to discard.
    ///
    /// The source default constructor also pins the schema version and file,
    /// which are the constants [`SCHEMA_VERSION`] and [`SCHEMA_FILE`] here.
    pub fn new() -> Self {
        Self::default()
    }

    /// Non-mutable access to the loading and storing options, the source
    /// `getOptions() const`.
    pub fn options(&self) -> &PeakFileOptions {
        &self.options
    }

    /// Mutable access to the loading and storing options, the source
    /// `getOptions()`.
    pub fn options_mut(&mut self) -> &mut PeakFileOptions {
        &mut self.options
    }

    /// Replace the loading and storing options, the source `setOptions`.
    pub fn set_options(&mut self, options: PeakFileOptions) {
        self.options = options;
    }

    /// The resource ceilings loading applies. Native addition; the source has
    /// no ceilings at all.
    pub fn limits(&self) -> &ReadLimits {
        &self.limits
    }

    /// Replace the resource ceilings loading applies.
    pub fn set_limits(&mut self, limits: ReadLimits) {
        self.limits = limits;
    }

    /// Whether [`MzDataFile::store`] discards what mzData cannot represent
    /// instead of refusing. `false` by default; see [`WriteOptions`].
    pub fn discards_unrepresentable(&self) -> bool {
        self.discard_unrepresentable
    }

    /// Choose whether [`MzDataFile::store`] discards or refuses.
    pub fn set_discard_unrepresentable(&mut self, discard: bool) {
        self.discard_unrepresentable = discard;
    }

    /// The schema version this adapter pins, the source `XMLFile::getVersion`.
    pub fn version(&self) -> &'static str {
        SCHEMA_VERSION
    }

    /// The write options [`MzDataFile::store`] derives from this adapter's
    /// state.
    pub fn write_options(&self) -> WriteOptions {
        WriteOptions {
            write_supplemental_data: self.options.write_supplemental_data,
            mz_32_bit: self.options.mz_32_bit,
            discard_unrepresentable: self.discard_unrepresentable,
        }
    }

    /// Load an mzData document, the source `MzDataFile::load`.
    ///
    /// The source signature is a template whose documentation says "`map` has
    /// to be a MSExperiment or have the same interface"; only `PeakMap` is
    /// ever instantiated, so this port is monomorphic on [`MSExperiment`].
    ///
    /// The returned experiment records the document path and the
    /// content-detected file type, as the source's `setLoadedFilePath` /
    /// `setLoadedFileType` pair does before parsing. Plain, gzip and bzip2
    /// input is accepted by content, independently of the suffix.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened — the source
    /// `Exception::FileNotFound` — and [`Error::Parse`],
    /// [`Error::InvalidValue`] or [`Error::Unsupported`] for a malformed,
    /// oversized or unrepresentable document, all of which the source reports
    /// as `Exception::ParseError`. The destination of
    /// [`MzDataFile::load_into`] is left unchanged on every error.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MSExperiment> {
        Ok(self.load_report(path)?.experiment)
    }

    /// Load an mzData document and keep what the load dropped or warned about.
    ///
    /// # Errors
    ///
    /// As [`MzDataFile::load`].
    pub fn load_report(&self, path: impl AsRef<Path>) -> Result<Loaded> {
        load_with_options(path, &self.options, &self.limits)
    }

    /// Replace `destination` only after the whole document parsed.
    ///
    /// # Errors
    ///
    /// As [`MzDataFile::load`].
    pub fn load_into(&self, path: impl AsRef<Path>, destination: &mut MSExperiment) -> Result<()> {
        *destination = self.load(path)?;
        Ok(())
    }

    /// Store an experiment as mzData, the source `MzDataFile::store`.
    ///
    /// As with [`MzDataFile::load`], the source is a template over anything
    /// with the `MSExperiment` interface and this port is monomorphic.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the output cannot be created — the source
    /// `Exception::UnableToCreateFile` — and [`Error::Unsupported`] when the
    /// experiment holds metadata mzData cannot represent and
    /// [`MzDataFile::discards_unrepresentable`] is false. The destination file
    /// is published only after the whole document has been serialised.
    pub fn store(&self, path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
        store_with_options(path, experiment, &self.write_options())
    }

    /// Store an experiment as mzData and keep what the store discarded.
    ///
    /// # Errors
    ///
    /// As [`MzDataFile::store`].
    pub fn store_report(
        &self,
        path: impl AsRef<Path>,
        experiment: &MSExperiment,
    ) -> Result<StoreReport> {
        store_report_with_options(path, experiment, &self.write_options())
    }

    /// Not ported: the source `isSemanticallyValid` validates a document
    /// against `/MAPPING/mzdata-mapping.xml` and `/CV/psi-mzdata.obo`.
    ///
    /// `MzDataFile_test.cpp:849-854` marks its own section `NOT_TESTABLE` with
    /// the comment that this "is not officially supported - the mapping file
    /// was hand-crafted", and neither resource ships with this crate. The
    /// generic machinery exists in
    /// [`crate::format::semantic_validator`] behind the
    /// `semantic-validation` feature and can be pointed at an mzData mapping
    /// by a caller that has one. Its two `@param[out]` lists, `errors` and
    /// `warnings`, would be returned rather than filled in; the source's
    /// `@exception Exception::FileNotFound` would be [`Error::Io`].
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`].
    pub fn is_semantically_valid(&self, _path: impl AsRef<Path>) -> Result<bool> {
        Err(Error::Unsupported(
            "mzData semantic validation needs the unshipped mzdata-mapping.xml and psi-mzdata.obo"
                .into(),
        ))
    }
}

// ---------------------------------------------------------------------------
// Free entry points
// ---------------------------------------------------------------------------

/// Load an mzData document with default options and ceilings.
///
/// # Errors
///
/// As [`MzDataFile::load`].
pub fn load(path: impl AsRef<Path>) -> Result<MSExperiment> {
    Ok(load_with_options(path, &PeakFileOptions::default(), &ReadLimits::default())?.experiment)
}

/// Load an mzData document under explicit scientific options and ceilings.
///
/// `options` carries the source `PeakFileOptions` filters: `metadata_only`
/// stops the parse at `<spectrumList>`, the MS-level set and the retention-time
/// and precursor-m/z ranges drop whole spectra, and the m/z and intensity
/// ranges drop individual peaks together with the aligned annotation values.
///
/// # Errors
///
/// As [`MzDataFile::load`].
pub fn load_with_options(
    path: impl AsRef<Path>,
    options: &PeakFileOptions,
    limits: &ReadLimits,
) -> Result<Loaded> {
    let path = path.as_ref();
    let text = path
        .to_str()
        .ok_or_else(|| Error::InvalidValue("mzData filename is not UTF-8".into()))?;
    let mut document = crate::metadata::DocumentIdentifier::new();
    document.set_loaded_file_path(text)?;
    document.set_loaded_file_type(path)?;
    let mut loaded = read_with_options(crate::format::path_io::open(path)?, options, limits)?;
    loaded.experiment.settings.document.loaded_file_path = document.loaded_file_path;
    loaded.experiment.settings.document.loaded_file_type = document.loaded_file_type;
    Ok(loaded)
}

/// Replace `destination` with a loaded mzData document, only on success.
///
/// # Errors
///
/// As [`MzDataFile::load`].
pub fn load_into(path: impl AsRef<Path>, destination: &mut MSExperiment) -> Result<()> {
    *destination = load(path)?;
    Ok(())
}

/// Read an mzData document from a stream with default options and ceilings.
///
/// Stream reads leave the document path and type unset, as the crate's other
/// stream readers do; [`load`] records both.
///
/// # Errors
///
/// As [`MzDataFile::load`], minus the open failure.
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    Ok(read_with_options(reader, &PeakFileOptions::default(), &ReadLimits::default())?.experiment)
}

/// Read an mzData document from a stream under explicit options and ceilings.
///
/// # Errors
///
/// As [`load_with_options`], minus the open failure.
pub fn read_with_options(
    reader: impl BufRead,
    options: &PeakFileOptions,
    limits: &ReadLimits,
) -> Result<Loaded> {
    let bytes = slurp(reader, limits)?;
    let decoded = transcode(&bytes)?;
    let mut parser = Parser::new(options, limits);
    parser.run(&decoded)?;
    Ok(parser.finish())
}

/// Store an experiment as mzData at `path`, refusing to discard.
///
/// # Errors
///
/// As [`MzDataFile::store`].
pub fn store(path: impl AsRef<Path>, experiment: &MSExperiment) -> Result<()> {
    store_with_options(path, experiment, &WriteOptions::default())
}

/// Store an experiment as mzData at `path` under explicit options.
///
/// `.gz` and `.bz2` suffixes select outer compression, as every other writer in
/// this crate does; the source `XMLFile::save_` does the same. The destination
/// is replaced only after the whole document has been serialised and flushed.
///
/// # Errors
///
/// As [`MzDataFile::store`].
pub fn store_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    store_report_with_options(path, experiment, options).map(|_| ())
}

/// Store an experiment as mzData and keep what the store discarded.
///
/// # Errors
///
/// As [`MzDataFile::store`].
pub fn store_report_with_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<StoreReport> {
    // The whole document is serialised into memory first so that a refusal
    // reaches the caller before any file is created, and so that the temporary
    // file the publication uses is written in one pass.
    let mut buffer = Vec::new();
    let report = write_with_options(&mut buffer, experiment, options)?;
    crate::format::path_io::write(path.as_ref(), |writer| {
        writer.write_all(&buffer)?;
        Ok(())
    })?;
    Ok(report)
}

/// Write an experiment as mzData to a stream, refusing to discard.
///
/// # Errors
///
/// As [`MzDataFile::store`].
pub fn write(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_with_options(writer, experiment, &WriteOptions::default()).map(|_| ())
}

/// Write an experiment as mzData to a stream under explicit options.
///
/// The layout is `MzDataHandler::writeTo` verbatim — tab indentation, `\n` line
/// endings, the `ISO-8859-1` declaration and the element order — with the two
/// documented corrections: text is XML-escaped, and any character outside
/// ISO-8859-1 is written as a numeric character reference so the declared
/// encoding stays truthful.
///
/// # Errors
///
/// As [`MzDataFile::store`].
pub fn write_with_options(
    writer: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<StoreReport> {
    let mut report = StoreReport::default();
    preflight_store(experiment, options, &mut report)?;
    let mut sink = writer;
    write_document(&mut sink, experiment, options, &mut report)?;
    Ok(report)
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn limit(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// `DRange<1>::encloses` (`DRange.h:152-161`): inclusive at the minimum,
/// exclusive at the maximum, and false for a NaN bound or value.
fn encloses(range: NumericRange, value: f64) -> bool {
    !(value < range.min || value >= range.max)
}

/// The text a C++ `std::ostream` with default classic-locale settings produces
/// for a double, which is what every numeric `writeCVS_` and attribute in
/// `writeTo` goes through.
fn stream_float(value: f64) -> Result<String> {
    crate::param::ParamValue::Float(value).to_stream_text()
}

/// `DataValue::operator<<`, the stringification `writeUserParam_` applies to
/// every metadata value.
fn meta_text(value: &MetaValue) -> Result<String> {
    fn bracket(parts: Vec<String>) -> String {
        format!("[{}]", parts.join(", "))
    }
    Ok(match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => stream_float(*number)?,
        MetaValueData::StringList(values) => bracket(values.clone()),
        MetaValueData::IntegerList(values) => {
            bracket(values.iter().map(|v| v.to_string()).collect())
        }
        MetaValueData::FloatList(values) => bracket(
            values
                .iter()
                .map(|v| stream_float(*v))
                .collect::<Result<Vec<_>>>()?,
        ),
    })
}

/// Position of `value` in its enum's source-ordered `ALL` array, which is the
/// index `cv_terms_` tables are built to agree with.
fn enum_index<T: PartialEq>(all: &[T], value: &T) -> Option<usize> {
    all.iter().position(|item| item == value)
}

/// The `cv_terms_` term for an enum value, or `None` when the table has no
/// entry for it or the entry is the empty placeholder that `writeCVS_` skips.
fn cv_term<'a, T: PartialEq>(table: &'a [&'a str], all: &[T], value: &T) -> Option<&'a str> {
    let index = enum_index(all, value)?;
    table.get(index).copied().filter(|term| !term.is_empty())
}

/// The enum value a `cvParam` term names, following `cvStringToEnum_`: the
/// index of the term in its table, or `fallback` with a warning when the table
/// has no such term.
///
/// `cvStringToEnum_` falls back to index 0, which is the enum's unknown value
/// for every table except `cv_terms_[18]`, where index 0 is `CID`; callers
/// pass that value as `fallback` so the difference is visible at the call site.
fn cv_enum<T: Copy>(table: &[&str], all: &[T], term: &str, fallback: T) -> (T, bool) {
    match table.iter().position(|entry| *entry == term) {
        // Every table is a prefix of its enum's `ALL`, so the index is always
        // in range; the fallback keeps that assumption from becoming a panic.
        Some(index) => (all.get(index).copied().unwrap_or(fallback), true),
        None => (fallback, false),
    }
}

/// Read the whole input, refusing anything above [`ReadLimits::max_xml_bytes`].
///
/// mzData payloads are split across arbitrarily many character-data chunks that
/// have to be concatenated before they can be decoded, so the document is held
/// in memory either way; reading it up front is also what lets the ISO-8859-1
/// declaration be honoured.
fn slurp(reader: impl BufRead, limits: &ReadLimits) -> Result<Vec<u8>> {
    let ceiling = limits
        .max_xml_bytes
        .checked_add(1)
        .ok_or_else(|| limit("mzData byte ceiling must be below u64::MAX"))?;
    let mut bytes = Vec::new();
    reader.take(ceiling).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > limits.max_xml_bytes {
        return Err(limit("mzData document exceeds the configured byte ceiling"));
    }
    Ok(bytes)
}

/// Honour the XML declaration's encoding.
///
/// mzData's own writer emits `encoding="ISO-8859-1"` and every upstream fixture
/// declares it, so a port that assumed UTF-8 would corrupt any accented
/// metadata. ISO-8859-1 maps byte `b` to `U+00b`, so the transcode is exact and
/// is skipped entirely when the document happens to be pure ASCII, which is the
/// common case. `US-ASCII` and `UTF-8` are accepted as they lie; anything else
/// is refused rather than silently misread.
fn transcode(bytes: &[u8]) -> Result<std::borrow::Cow<'_, [u8]>> {
    let mut reader = Reader::from_reader(bytes);
    let mut buffer = Vec::new();
    let declared = match reader.read_event_into(&mut buffer) {
        Ok(Event::Decl(declaration)) => declaration
            .encoding()
            .transpose()
            .map_err(|e| parse_error(e.to_string()))?
            .map(|encoding| encoding.to_vec()),
        _ => None,
    };
    let latin1 = match declared.as_deref() {
        None => false,
        Some(encoding) if encoding.eq_ignore_ascii_case(b"UTF-8") => false,
        Some(encoding) if encoding.eq_ignore_ascii_case(b"US-ASCII") => false,
        Some(encoding)
            if encoding.eq_ignore_ascii_case(b"ISO-8859-1")
                || encoding.eq_ignore_ascii_case(b"latin1") =>
        {
            true
        }
        Some(_) => {
            return Err(Error::Unsupported(
                "mzData declares an encoding other than UTF-8, US-ASCII or ISO-8859-1".into(),
            ));
        }
    };
    if !latin1 || bytes.is_ascii() {
        return Ok(std::borrow::Cow::Borrowed(bytes));
    }
    let mut text = String::new();
    text.try_reserve(bytes.len())
        .map_err(|_| limit("cannot allocate the transcoded mzData document"))?;
    for &byte in bytes {
        text.push(byte as char);
    }
    Ok(std::borrow::Cow::Owned(text.into_bytes()))
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// Which spectrum-level array a `<data>` element belongs to.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ArrayKind {
    Mz,
    Intensity,
    /// The payload of the `<supDataArrayBinary>` at this position among the
    /// spectrum's float data arrays. Carrying the index is what makes the
    /// payload-to-array pairing explicit instead of positional.
    Supplemental(usize),
}

/// One `<data>` element: its declared transport and its accumulated payload.
struct Encoded {
    kind: ArrayKind,
    precision: Precision,
    endian: Endian,
    /// The `length` attribute, which the source header documents as "the number
    /// of peaks in the current spectrum (according to the length attribute --
    /// which should not be trusted)". It is spent as a resource ceiling and
    /// then compared with the payload, which wins.
    declared: Option<usize>,
    payload: String,
}

/// One `<supDesc supDataArrayRef>` and the description it carries.
struct Description {
    reference: String,
    metadata: MetaInfo,
}

struct Parser<'a> {
    options: &'a PeakFileOptions,
    limits: &'a ReadLimits,
    report: LoadReport,
    experiment: MSExperiment,
    /// Open element names, innermost last.
    stack: Vec<String>,
    /// Character data of the innermost open element.
    text: String,
    /// Remaining character-data budget.
    text_budget: usize,
    /// Remaining metadata-entry budget.
    metadata_budget: usize,
    /// Remaining retained-peak budget.
    peak_budget: usize,
    /// The document-wide `<dataProcessing>` record, shared by every spectrum.
    processing: Option<DataProcessing>,
    /// The spectrum being assembled.
    spectrum: MSSpectrum,
    /// `<data>` elements of the open spectrum, in document order.
    arrays: Vec<Encoded>,
    /// `<supDesc>` entries of the open spectrum.
    descriptions: Vec<Description>,
    /// The source `skip_spectrum_`.
    skip: bool,
    /// Whether `<spectrumList>` was reached, which ends a metadata-only load.
    stop: bool,
}

impl<'a> Parser<'a> {
    fn new(options: &'a PeakFileOptions, limits: &'a ReadLimits) -> Self {
        Self {
            options,
            limits,
            report: LoadReport::default(),
            experiment: MSExperiment::new(),
            stack: Vec::new(),
            text: String::new(),
            text_budget: limits.max_text_bytes,
            metadata_budget: limits.max_metadata_entries,
            peak_budget: limits.max_total_peaks,
            processing: None,
            spectrum: MSSpectrum::default(),
            arrays: Vec::new(),
            descriptions: Vec::new(),
            skip: false,
            stop: false,
        }
    }

    fn finish(self) -> Loaded {
        Loaded {
            experiment: self.experiment,
            report: self.report,
        }
    }

    fn warn(&mut self, message: impl Into<String>) {
        self.report.warning_count += 1;
        if self.report.warnings.len() < self.limits.max_warnings {
            self.report.warnings.push(message.into());
        }
    }

    /// The element enclosing the innermost open one, the source `parent_tag`.
    /// Called while the current element is still on the stack.
    fn parent(&self) -> &str {
        self.stack
            .len()
            .checked_sub(2)
            .and_then(|index| self.stack.get(index))
            .map_or("", String::as_str)
    }

    fn run(&mut self, bytes: &[u8]) -> Result<()> {
        let mut reader = Reader::from_reader(bytes);
        reader.config_mut().expand_empty_elements = true;
        reader.config_mut().enable_all_checks(true);
        let mut buffer = Vec::new();
        let mut seen_root = false;
        loop {
            let event = reader
                .read_event_into(&mut buffer)
                .map_err(|e| parse_error(e.to_string()))?;
            match event {
                Event::Start(element) => {
                    let tag = local_name(&element)?;
                    if !seen_root {
                        if tag != "mzData" {
                            return Err(parse_error("mzData root element is not <mzData>"));
                        }
                        seen_root = true;
                    }
                    if self.stack.len() >= self.limits.max_xml_depth {
                        return Err(limit("mzData nesting exceeds the configured depth ceiling"));
                    }
                    self.stack.push(tag.clone());
                    self.text.clear();
                    self.start(&tag, &element)?;
                    if self.stop {
                        return Ok(());
                    }
                }
                Event::End(element) => {
                    let closing = local_name_end(&element)?;
                    let tag = self
                        .stack
                        .pop()
                        .ok_or_else(|| parse_error("unbalanced mzData end element"))?;
                    if tag != closing {
                        return Err(parse_error(format!(
                            "mzData closes <{closing}> while <{tag}> is open"
                        )));
                    }
                    // The element has been popped, so its parent is now the
                    // innermost open element.
                    let parent = self.stack.last().cloned().unwrap_or_default();
                    self.end(&tag, &parent)?;
                    self.text.clear();
                }
                Event::Text(text) => {
                    let decoded = text.decode().map_err(|e| parse_error(e.to_string()))?;
                    self.push_text(&decoded)?;
                }
                // Xerces delivers a CDATA section's content through
                // `characters()` like any other character data, so the source
                // handler cannot tell the two apart and neither does this.
                Event::CData(text) => {
                    let decoded = text.decode().map_err(|e| parse_error(e.to_string()))?;
                    self.push_text(&decoded)?;
                }
                // quick-xml reports every `&…;` in character data as its own
                // event. The five predefined entities and numeric character
                // references are resolved; an external or undeclared entity is
                // refused rather than expanded, so no document can pull in a
                // file this parser was not given.
                Event::GeneralRef(reference) => {
                    let resolved = if let Some(c) = reference
                        .resolve_char_ref()
                        .map_err(|e| parse_error(e.to_string()))?
                    {
                        c
                    } else {
                        let name = reference.decode().map_err(|e| parse_error(e.to_string()))?;
                        match name.as_ref() {
                            "amp" => '&',
                            "lt" => '<',
                            "gt" => '>',
                            "apos" => '\'',
                            "quot" => '"',
                            other => {
                                return Err(Error::Unsupported(format!(
                                    "mzData names the undeclared XML entity '{other}'"
                                )));
                            }
                        }
                    };
                    self.push_text(resolved.encode_utf8(&mut [0u8; 4]))?;
                }
                Event::Eof => break,
                Event::Decl(_) | Event::Comment(_) | Event::PI(_) => {}
                // `expand_empty_elements` turns `<x/>` into a Start/End pair,
                // so this variant cannot be produced; it is an error rather
                // than a panic because the input that would reach it is
                // file-derived.
                Event::Empty(_) => {
                    return Err(parse_error("unexpanded empty mzData element"));
                }
                Event::DocType(_) => {
                    return Err(Error::Unsupported("XML DTDs are not supported".into()));
                }
            }
        }
        if !seen_root {
            return Err(parse_error("input is not an mzData document"));
        }
        if !self.stack.is_empty() {
            return Err(parse_error("mzData document ends with open elements"));
        }
        Ok(())
    }

    /// Accumulate character data for the innermost open element, charging it
    /// against [`ReadLimits::max_text_bytes`] whether or not it is kept.
    fn push_text(&mut self, text: &str) -> Result<()> {
        self.text_budget = self
            .text_budget
            .checked_sub(text.len())
            .ok_or_else(|| limit("mzData character data exceeds its ceiling"))?;
        if !self.skip {
            self.text.push_str(text);
        }
        Ok(())
    }

    fn spend_metadata(&mut self) -> Result<()> {
        self.metadata_budget = self
            .metadata_budget
            .checked_sub(1)
            .ok_or_else(|| limit("mzData metadata exceeds its entry ceiling"))?;
        Ok(())
    }

    // -- start elements -----------------------------------------------------

    fn start(&mut self, tag: &str, element: &BytesStart<'_>) -> Result<()> {
        // `MzDataHandler.cpp:214-217`: once a spectrum is being skipped nothing
        // but the next `<spectrum>` is interpreted.
        if self.skip && tag != "spectrum" {
            return Ok(());
        }
        let parent = self.parent().to_owned();
        match tag {
            "mzData" => {
                self.experiment.settings.document.identifier =
                    required(element, "accessionNumber")?;
            }
            "sourceFile" => self
                .experiment
                .settings
                .source_files
                .push(SourceFile::default()),
            "contact" => self
                .experiment
                .settings
                .contacts
                .push(ContactPerson::default()),
            "source" => {
                self.experiment.settings.instrument.ion_sources =
                    vec![crate::metadata::IonSource::default()];
            }
            "detector" => {
                self.experiment.settings.instrument.ion_detectors =
                    vec![crate::metadata::IonDetector::default()];
            }
            "analyzer" => self
                .experiment
                .settings
                .instrument
                .mass_analyzers
                .push(crate::metadata::MassAnalyzer::default()),
            "software" => {
                let mut processing = DataProcessing::default();
                if let Some(stamp) = optional(element, "completionTime")? {
                    // `asDateTime_` (XMLHandler.h:359-377) trims, keeps the
                    // first 19 characters and logs a non-fatal error on
                    // failure, leaving the unset sentinel behind.
                    let trimmed = stamp.trim();
                    let head: String = trimmed.chars().take(19).collect();
                    match crate::data_structures::DateTime::parse(&head) {
                        Ok(value) => processing.completion_time = Some(value),
                        Err(_) => {
                            self.warn(format!("DateTime conversion error of \"{head}\""));
                        }
                    }
                }
                self.processing = Some(processing);
            }
            "precursor" => self.spectrum.precursors.push(Precursor::default()),
            "cvParam" => {
                let accession = required(element, "accession")?;
                let value = optional(element, "value")?.unwrap_or_default();
                self.cv_param(&parent, &accession, &value)?;
            }
            "supDataDesc" => {
                if let Some(comment) = optional(element, "comment")? {
                    self.spend_metadata()?;
                    let description = self
                        .descriptions
                        .last_mut()
                        .ok_or_else(|| parse_error("<supDataDesc> outside a <supDesc> element"))?;
                    description
                        .metadata
                        .insert("comment".into(), MetaValue::from(comment));
                }
            }
            "userParam" => {
                let name = required(element, "name")?;
                let value = optional(element, "value")?.unwrap_or_default();
                self.user_param(&parent, name, value)?;
            }
            "supDataArrayBinary" => {
                let id = required(element, "id")?;
                let metadata = self
                    .descriptions
                    .iter()
                    .find(|entry| entry.reference == id)
                    .map(|entry| entry.metadata.clone())
                    .unwrap_or_default();
                if self.spectrum.float_data_arrays.len() >= self.limits.max_arrays_per_spectrum {
                    return Err(limit("mzData spectrum exceeds its array ceiling"));
                }
                self.spectrum.float_data_arrays.push(DataArray {
                    metadata,
                    ..DataArray::<f32>::default()
                });
            }
            "spectrum" => {
                // A `<spectrum>` always starts clean, including the peak count
                // and the skip flag; the source leaves `peak_count_` from the
                // previous spectrum in place (`MzDataHandler.cpp:415` is the
                // only assignment).
                self.reset_spectrum();
                self.spectrum.native_id = format!("spectrum={}", required(element, "id")?);
                if let Some(processing) = &self.processing {
                    self.spectrum
                        .data_processing
                        .push(Arc::new(processing.clone()));
                }
            }
            "spectrumList" => {
                if self.options.metadata_only {
                    // `MzDataHandler.cpp:350`: `EndParsingSoftly` keeps the
                    // metadata parsed so far and abandons the rest.
                    self.stop = true;
                    return Ok(());
                }
                // `MzDataHandler.cpp:352-353` reads this straight into
                // `MSExperiment::reserve`. Here it is only a ceiling: a
                // declared count that disagrees with the number of children is
                // advisory on reading, as it is for every list in this crate's
                // mzML reader.
                let count: u64 = required(element, "count")?
                    .parse()
                    .map_err(|_| parse_error("<spectrumList count> is not a whole number"))?;
                if count > self.limits.max_spectra as u64 {
                    return Err(limit(
                        "<spectrumList count> exceeds the configured spectrum ceiling",
                    ));
                }
            }
            "acqSpecification" => {
                let kind = required(element, "spectrumType")?;
                self.spectrum.spectrum_type = match kind.as_str() {
                    "discrete" => SpectrumType::Centroid,
                    "continuous" => SpectrumType::Profile,
                    _ => {
                        self.warn(format!("Invalid spectrum type '{kind}'."));
                        SpectrumType::Unknown
                    }
                };
                self.spectrum.acquisition_info.method_of_combination =
                    required(element, "methodOfCombination")?;
            }
            "acquisition" => {
                if self.spectrum.acquisition_info.acquisitions.len()
                    >= self.limits.max_arrays_per_spectrum
                {
                    return Err(limit("mzData spectrum exceeds its acquisition ceiling"));
                }
                self.spectrum
                    .acquisition_info
                    .acquisitions
                    .push(Acquisition {
                        identifier: required(element, "acqNumber")?,
                        metadata: MetaInfo::new(),
                    });
            }
            "spectrumInstrument" | "acqInstrument" => {
                let level: u32 = required(element, "msLevel")?.parse().map_err(|_| {
                    parse_error("<spectrumInstrument msLevel> is not a whole number")
                })?;
                self.spectrum.ms_level = level;
                let mut window = ScanWindow::default();
                let start = optional(element, "mzRangeStart")?;
                let stop = optional(element, "mzRangeStop")?;
                if let Some(value) = &start {
                    window.begin = finite(value, "mzRangeStart")?;
                }
                if let Some(value) = &stop {
                    window.end = finite(value, "mzRangeStop")?;
                }
                // `MzDataHandler.cpp:392-396`: a window is kept only when at
                // least one bound is nonzero, so an absent `mzRangeStop`
                // leaves a window whose end is 0 and whose begin may exceed it.
                if window.begin != 0.0 || window.end != 0.0 {
                    self.spectrum.instrument_settings.scan_windows.push(window);
                }
                if self.options.has_ms_levels() {
                    let keep = i32::try_from(level)
                        .is_ok_and(|level| self.options.contains_ms_level(level));
                    if !keep {
                        self.skip = true;
                    }
                }
            }
            "supDesc" => {
                if self.descriptions.len() >= self.limits.max_arrays_per_spectrum {
                    return Err(limit("mzData spectrum exceeds its description ceiling"));
                }
                self.descriptions.push(Description {
                    reference: required(element, "supDataArrayRef")?,
                    metadata: MetaInfo::new(),
                });
            }
            "data" => {
                // `fillData_` identifies the m/z and intensity arrays by their
                // *position* among the `<data>` elements — `precisions_[0]` and
                // `precisions_[1]` (`MzDataHandler.cpp:522-527`) — so a
                // document that writes `<intenArrayBinary>` before
                // `<mzArrayBinary>` has its two arrays swapped. This port
                // classifies by parent element instead.
                let kind = match parent.as_str() {
                    "mzArrayBinary" => ArrayKind::Mz,
                    "intenArrayBinary" => ArrayKind::Intensity,
                    "supDataArrayBinary" => ArrayKind::Supplemental(
                        self.spectrum
                            .float_data_arrays
                            .len()
                            .checked_sub(1)
                            .ok_or_else(|| {
                                parse_error("<data> inside an unopened <supDataArrayBinary>")
                            })?,
                    ),
                    other => {
                        return Err(parse_error(format!(
                            "<data> inside unexpected element <{other}>"
                        )));
                    }
                };
                // Both attributes are required by `attributeAsString_`
                // (`MzDataHandler.cpp:413-414`); `length` is required by the
                // schema but only read for the m/z array, so it stays optional
                // here as well.
                let (precision, known_precision) =
                    Precision::parse(&required(element, "precision")?);
                if !known_precision {
                    self.warn("<data precision> is neither \"32\" nor \"64\"; reading it as 64");
                }
                let (endian, known_endian) = Endian::parse(&required(element, "endian")?);
                if !known_endian {
                    self.warn(
                        "<data endian> is neither \"little\" nor \"big\"; reading it as little",
                    );
                }
                let declared = match optional(element, "length")? {
                    Some(text) => {
                        let value: u64 = text
                            .parse()
                            .map_err(|_| parse_error("<data length> is not a whole number"))?;
                        if value > self.limits.max_array_elements as u64 {
                            return Err(limit(
                                "<data length> exceeds the configured element ceiling",
                            ));
                        }
                        Some(value as usize)
                    }
                    None => None,
                };
                if self.arrays.len() >= self.limits.max_arrays_per_spectrum {
                    return Err(limit("mzData spectrum exceeds its array ceiling"));
                }
                self.arrays.push(Encoded {
                    kind,
                    precision,
                    endian,
                    declared,
                    payload: String::new(),
                });
            }
            _ => {}
        }
        Ok(())
    }

    fn reset_spectrum(&mut self) {
        self.spectrum = MSSpectrum::default();
        self.arrays.clear();
        self.descriptions.clear();
        self.skip = false;
    }

    // -- end elements -------------------------------------------------------

    fn end(&mut self, tag: &str, parent: &str) -> Result<()> {
        if tag == "spectrum" {
            let skipped = self.skip;
            if skipped {
                self.report.spectra_skipped += 1;
            } else {
                self.fill_data()?;
                if self.experiment.spectra.len() >= self.limits.max_spectra {
                    return Err(limit("mzData exceeds the configured spectrum ceiling"));
                }
                let spectrum = std::mem::take(&mut self.spectrum);
                self.experiment.spectra.push(spectrum);
            }
            self.reset_spectrum();
            return Ok(());
        }
        if self.skip {
            return Ok(());
        }
        let text = std::mem::take(&mut self.text);
        match tag {
            "sampleName" => self.experiment.settings.sample.name = text,
            "instrumentName" => self.experiment.settings.instrument.name = text,
            "version" => {
                self.processing_mut()?.software.version = text;
            }
            "institution" => self.contact_mut()?.institution = text,
            "contactInfo" => self.contact_mut()?.contact_info = text,
            "name" if parent == "contact" => {
                let contact = self.contact_mut()?;
                contact.set_name(&text)?;
            }
            "name" if parent == "software" => self.processing_mut()?.software.name = text,
            "comments" if parent == "software" => {
                self.spend_metadata()?;
                self.processing_mut()?
                    .software
                    .cv_terms
                    .metadata
                    .insert("comment".into(), MetaValue::from(text));
            }
            // `spec_.setComment` (`MzDataHandler.cpp:137`) has no `MSSpectrum`
            // counterpart in this crate — the comment lives on
            // `SpectrumSettings`, which `MSSpectrum` does not embed — so it is
            // kept as a `comment` metadata entry. The source writer never
            // emits this element, so the value is read-only either way.
            "comments" if parent == "spectrumDesc" => {
                self.spend_metadata()?;
                self.spectrum
                    .metadata
                    .insert("comment".into(), MetaValue::from(text));
            }
            "data" => {
                let array = self
                    .arrays
                    .last_mut()
                    .ok_or_else(|| parse_error("</data> without an open array"))?;
                // `MzDataHandler.cpp:492-494`: "line breaks inside the base64
                // data are unfortunately no exception".
                array
                    .payload
                    .extend(text.chars().filter(|c| !c.is_whitespace()));
                if array.payload.len() / 4 * 3 > self.limits.max_array_bytes {
                    return Err(limit(
                        "mzData binary array exceeds the configured byte ceiling",
                    ));
                }
            }
            // The source pushes the decode slot when `<arrayName>` opens and
            // names the array from its character data
            // (`MzDataHandler.cpp:426-429`, `:144-147`), so a
            // `<supDataArrayBinary>` without an `<arrayName>` child shifts
            // every following payload by one slot. Here the array is created
            // by its own element and only named here, so a missing
            // `<arrayName>` leaves the name empty and nothing shifts.
            "arrayName" if parent == "supDataArrayBinary" => {
                if let Some(array) = self.spectrum.float_data_arrays.last_mut() {
                    array.name = text;
                }
            }
            "nameOfFile" if parent == "sourceFile" => self.source_file_mut()?.name = text,
            "pathToFile" if parent == "sourceFile" => self.source_file_mut()?.path = text,
            "fileType" if parent == "sourceFile" => self.source_file_mut()?.file_type = text,
            // `MzDataHandler.cpp:152-171`: the `supSourceFile` children are
            // parsed and deliberately dropped.
            "nameOfFile" | "pathToFile" | "fileType" if parent == "supSourceFile" => {}
            _ => {
                if !text.trim().is_empty() {
                    self.warn(format!(
                        "Unhandled character content in tag '{tag}': {}",
                        text.trim()
                    ));
                }
            }
        }
        Ok(())
    }

    fn contact_mut(&mut self) -> Result<&mut ContactPerson> {
        self.experiment
            .settings
            .contacts
            .last_mut()
            .ok_or_else(|| parse_error("contact detail outside a <contact> element"))
    }

    fn source_file_mut(&mut self) -> Result<&mut SourceFile> {
        self.experiment
            .settings
            .source_files
            .last_mut()
            .ok_or_else(|| parse_error("source file detail outside a <sourceFile> element"))
    }

    fn processing_mut(&mut self) -> Result<&mut DataProcessing> {
        self.processing
            .as_mut()
            .ok_or_else(|| parse_error("software detail outside a <software> element"))
    }

    // -- parameters ---------------------------------------------------------

    fn user_param(&mut self, parent: &str, name: String, value: String) -> Result<()> {
        self.spend_metadata()?;
        let value = MetaValue::from(value);
        let target = match parent {
            "spectrumInstrument" => &mut self.spectrum.instrument_settings.metadata,
            "acquisition" => {
                &mut self
                    .spectrum
                    .acquisition_info
                    .acquisitions
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside an <acquisition> element"))?
                    .metadata
            }
            "ionSelection" | "activation" => {
                // Both write onto the precursor, and the source stores them in
                // the precursor's MetaInfo rather than its CV term list.
                let precursor = self
                    .spectrum
                    .precursors
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside a <precursor> element"))?;
                &mut precursor.cv_terms.metadata
            }
            "supDataDesc" => {
                &mut self
                    .descriptions
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside a <supDesc> element"))?
                    .metadata
            }
            "detector" => {
                &mut self
                    .experiment
                    .settings
                    .instrument
                    .ion_detectors
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside a <detector> element"))?
                    .metadata
            }
            "source" => {
                &mut self
                    .experiment
                    .settings
                    .instrument
                    .ion_sources
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside a <source> element"))?
                    .metadata
            }
            "sampleDescription" => &mut self.experiment.settings.sample.metadata,
            "analyzer" => {
                &mut self
                    .experiment
                    .settings
                    .instrument
                    .mass_analyzers
                    .last_mut()
                    .ok_or_else(|| parse_error("<userParam> outside an <analyzer> element"))?
                    .metadata
            }
            "additional" => &mut self.experiment.settings.instrument.metadata,
            "processingMethod" => &mut self.processing_mut()?.metadata,
            other => {
                self.warn(format!(
                    "Invalid userParam: name=\"{name}, value=\"{}\"",
                    meta_text(&value)?
                ));
                let _ = other;
                return Ok(());
            }
        };
        target.insert(name, value);
        Ok(())
    }

    fn cv_param(&mut self, parent: &str, accession: &str, value: &str) -> Result<()> {
        // The source builds a `error` string naming the section and warns once
        // at the end when a term is not recognised inside a known section;
        // an unknown section warns immediately (`MzDataHandler.cpp:1404-1414`).
        let section = match parent {
            "spectrumInstrument" => self.cv_spectrum_instrument(accession, value)?,
            "ionSelection" => self.cv_ion_selection(accession, value)?,
            "activation" => self.cv_activation(accession, value)?,
            "supDataDesc" => Some("supDataDesc.UserParam"),
            "acquisition" => {
                Some("spectrumDesc.spectrumSettings.acquisitionSpecification.acquisition.UserParam")
            }
            "detector" => self.cv_detector(accession, value)?,
            "source" => self.cv_source(accession, value)?,
            "sampleDescription" => self.cv_sample(accession, value)?,
            "analyzer" => self.cv_analyzer(accession, value)?,
            "additional" => self.cv_additional(accession, value)?,
            "processingMethod" => self.cv_processing(accession)?,
            other => {
                self.warn(format!(
                    "Unexpected cvParam: accession=\"{accession}\" value=\"{value}\" in tag {other}"
                ));
                None
            }
        };
        if let Some(section) = section {
            self.warn(format!(
                "Invalid cvParam: accession=\"{accession}\" value=\"{value}\" in {section}"
            ));
        }
        Ok(())
    }

    fn cv_spectrum_instrument(
        &mut self,
        accession: &str,
        value: &str,
    ) -> Result<Option<&'static str>> {
        match accession {
            "PSI:1000036" => {
                // `EnhancedResolutionScan` and `Zoom` are the same state.
                let zoom = matches!(value, "Zoom" | "EnhancedResolutionScan");
                let mode = match value {
                    "Zoom" | "EnhancedResolutionScan" | "MassScan" => Some(ScanMode::MassSpectrum),
                    "SelectedIonDetection" => Some(ScanMode::SelectedIonMonitoring),
                    "SelectedReactionMonitoring" => Some(ScanMode::SelectedReactionMonitoring),
                    "ConsecutiveReactionMonitoring" => {
                        Some(ScanMode::ConsecutiveReactionMonitoring)
                    }
                    "ConstantNeutralGainScan" => Some(ScanMode::ConstantNeutralGain),
                    "ConstantNeutralLossScan" => Some(ScanMode::ConstantNeutralLoss),
                    "ProductIonScan" => None,
                    "PrecursorIonScan" => Some(ScanMode::Precursor),
                    _ => None,
                };
                if zoom {
                    self.spectrum.instrument_settings.zoom_scan = true;
                }
                match (mode, value) {
                    (Some(mode), _) => self.spectrum.instrument_settings.scan_mode = mode,
                    (None, "ProductIonScan") => {
                        self.spectrum.instrument_settings.scan_mode = ScanMode::MsnSpectrum;
                        // The source overrides the `msLevel` attribute here.
                        self.spectrum.ms_level = 2;
                    }
                    (None, other) => {
                        // `MzDataHandler.cpp:1136-1140` writes `MSNSPECTRUM`
                        // onto `exp_->getSpectra().back()` here — the
                        // *previous* spectrum, and `.back()` on an empty
                        // vector when this is the first one, which is
                        // undefined behaviour. This port sets the mode on the
                        // spectrum the parameter belongs to and always warns;
                        // the source warns only on the MS1 branch.
                        let level = self.spectrum.ms_level;
                        if level >= 2 {
                            self.spectrum.instrument_settings.scan_mode = ScanMode::MsnSpectrum;
                            self.warn(format!(
                                "Unknown scan mode '{other}' on an MS{level} spectrum. Assuming MSn scan"
                            ));
                        } else {
                            self.spectrum.instrument_settings.scan_mode = ScanMode::MassSpectrum;
                            self.warn(format!("Unknown scan mode '{other}'. Assuming full scan"));
                        }
                    }
                }
                Ok(None)
            }
            // PSI:1000038 is minutes, PSI:1000039 seconds; OpenMS stores
            // seconds (`MzDataHandler.cpp:1143`, `:1152`).
            "PSI:1000038" | "PSI:1000039" => {
                let seconds = finite(value, "retention time")?
                    * if accession == "PSI:1000038" {
                        60.0
                    } else {
                        1.0
                    };
                self.spectrum.rt = seconds;
                if self.options.has_rt_range() && !encloses(self.options.rt_range(), seconds) {
                    self.skip = true;
                }
                Ok(None)
            }
            "PSI:1000037" => {
                // The source accepts three spellings each way and warns
                // otherwise (`MzDataHandler.cpp:1159-1174`).
                match value {
                    "Positive" | "positive" | "+" => {
                        self.spectrum.instrument_settings.polarity = Polarity::Positive;
                    }
                    "Negative" | "negative" | "-" => {
                        self.spectrum.instrument_settings.polarity = Polarity::Negative;
                    }
                    other => self.warn(format!(
                        "Invalid scan polarity (PSI:1000037) detected: \"{other}\". Valid are 'Positive' or 'Negative'."
                    )),
                }
                Ok(None)
            }
            _ => Ok(Some(
                "SpectrumDescription.SpectrumSettings.SpectrumInstrument",
            )),
        }
    }

    fn cv_ion_selection(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        match accession {
            "PSI:1000040" => {
                let mz = finite(value, "precursor m/z")?;
                self.precursor_mut()?.mz = mz;
                if self.options.has_precursor_mz_range()
                    && !encloses(self.options.precursor_mz_range(), mz)
                {
                    self.skip = true;
                }
                Ok(None)
            }
            "PSI:1000041" => {
                let charge: i32 = value
                    .parse()
                    .map_err(|_| parse_error(format!("Int conversion error of \"{value}\"")))?;
                let precursor = self.precursor_mut()?;
                if precursor.charge != 0 {
                    precursor.charge = 0;
                    self.warn(format!(
                        "Multiple precursor charges detected, expected only one! Ignoring this charge settings! accession=\"{accession}\", value=\"{value}\""
                    ));
                } else {
                    precursor.charge = charge;
                }
                Ok(None)
            }
            "PSI:1000042" => {
                let magnitude = intensity(finite(value, "precursor intensity")?)?;
                self.precursor_mut()?.intensity = magnitude;
                Ok(None)
            }
            // PSI:1000043 is the intensity unit, deliberately ignored.
            "PSI:1000043" => Ok(None),
            _ => Ok(Some("PrecursorList.Precursor.IonSelection.UserParam")),
        }
    }

    fn cv_activation(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        match accession {
            "PSI:1000044" => {
                let (method, known) = cv_enum(
                    ACTIVATION_METHOD,
                    ActivationMethod::ALL,
                    value,
                    ActivationMethod::Cid,
                );
                if !known {
                    self.warn(format!("Unexpected CV entry 'activation method'='{value}'"));
                }
                self.precursor_mut()?.activation_methods.insert(method);
                Ok(None)
            }
            "PSI:1000045" => {
                let energy = finite(value, "activation energy")?;
                self.precursor_mut()?.activation_energy = energy;
                Ok(None)
            }
            // PSI:1000046 is the energy unit; the source assumes electronvolt.
            "PSI:1000046" => Ok(None),
            _ => Ok(Some("PrecursorList.Precursor.Activation.UserParam")),
        }
    }

    fn precursor_mut(&mut self) -> Result<&mut Precursor> {
        self.spectrum
            .precursors
            .last_mut()
            .ok_or_else(|| parse_error("precursor parameter outside a <precursor> element"))
    }

    fn cv_detector(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        let mut unexpected = None;
        {
            let detector = self
                .experiment
                .settings
                .instrument
                .ion_detectors
                .last_mut()
                .ok_or_else(|| parse_error("detector parameter outside a <detector> element"))?;
            match accession {
                "PSI:1000026" => {
                    let (kind, known) = cv_enum(
                        DETECTOR_TYPE,
                        DetectorType::ALL,
                        value,
                        DetectorType::Unknown,
                    );
                    detector.detector_type = kind;
                    if !known {
                        unexpected = Some("detector type");
                    }
                }
                "PSI:1000027" => {
                    let (mode, known) = cv_enum(
                        ACQUISITION_MODE,
                        DetectorAcquisitionMode::ALL,
                        value,
                        DetectorAcquisitionMode::Unknown,
                    );
                    detector.acquisition_mode = mode;
                    if !known {
                        unexpected = Some("acquisition mode");
                    }
                }
                "PSI:1000028" => detector.resolution = finite(value, "detector resolution")?,
                "PSI:1000029" => {
                    detector.adc_sampling_frequency = finite(value, "sampling frequency")?;
                }
                _ => return Ok(Some("Description.Instrument.Detector.UserParam")),
            }
        }
        if let Some(what) = unexpected {
            self.warn(format!("Unexpected CV entry '{what}'='{value}'"));
        }
        Ok(None)
    }

    fn cv_source(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        let mut unexpected = None;
        {
            let source = self
                .experiment
                .settings
                .instrument
                .ion_sources
                .last_mut()
                .ok_or_else(|| parse_error("source parameter outside a <source> element"))?;
            match accession {
                "PSI:1000007" => {
                    let (inlet, known) =
                        cv_enum(INLET_TYPE, InletType::ALL, value, InletType::Unknown);
                    source.inlet_type = inlet;
                    if !known {
                        unexpected = Some("inlet type");
                    }
                }
                "PSI:1000008" => {
                    let (method, known) = cv_enum(
                        IONIZATION_TYPE,
                        IonizationMethod::ALL,
                        value,
                        IonizationMethod::Unknown,
                    );
                    source.ionization_method = method;
                    if !known {
                        unexpected = Some("ion source");
                    }
                }
                "PSI:1000009" => {
                    let (polarity, known) =
                        cv_enum(IONIZATION_MODE, Polarity::ALL, value, Polarity::Unknown);
                    source.polarity = polarity;
                    if !known {
                        unexpected = Some("polarity");
                    }
                }
                _ => return Ok(Some("Description.Instrument.Source.UserParam")),
            }
        }
        if let Some(what) = unexpected {
            self.warn(format!("Unexpected CV entry '{what}'='{value}'"));
        }
        Ok(None)
    }

    fn cv_sample(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        match accession {
            "PSI:1000001" => self.experiment.settings.sample.number = value.to_owned(),
            "PSI:1000003" => {
                let (state, known) =
                    cv_enum(SAMPLE_STATE, SampleState::ALL, value, SampleState::Unknown);
                self.experiment.settings.sample.state = state;
                if !known {
                    self.warn(format!("Unexpected CV entry 'sample state'='{value}'"));
                }
            }
            "PSI:1000004" => {
                self.experiment.settings.sample.mass = finite(value, "sample mass")?;
            }
            "PSI:1000005" => {
                self.experiment.settings.sample.volume = finite(value, "sample volume")?;
            }
            "PSI:1000006" => {
                self.experiment.settings.sample.concentration =
                    finite(value, "sample concentration")?;
            }
            _ => return Ok(Some("Description.Admin.SampleDescription.UserParam")),
        }
        Ok(None)
    }

    fn cv_analyzer(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        let mut unexpected = None;
        {
            let analyzer = self
                .experiment
                .settings
                .instrument
                .mass_analyzers
                .last_mut()
                .ok_or_else(|| parse_error("analyzer parameter outside an <analyzer> element"))?;
            match accession {
                "PSI:1000010" => {
                    let (kind, known) = cv_enum(
                        ANALYZER_TYPE,
                        AnalyzerType::ALL,
                        value,
                        AnalyzerType::Unknown,
                    );
                    analyzer.analyzer_type = kind;
                    if !known {
                        unexpected = Some("analyzer type");
                    }
                }
                "PSI:1000011" => analyzer.resolution = finite(value, "resolution")?,
                "PSI:1000012" => {
                    let (method, known) = cv_enum(
                        RESOLUTION_METHOD,
                        ResolutionMethod::ALL,
                        value,
                        ResolutionMethod::Unknown,
                    );
                    analyzer.resolution_method = method;
                    if !known {
                        unexpected = Some("resolution method");
                    }
                }
                "PSI:1000013" => {
                    let (kind, known) = cv_enum(
                        RESOLUTION_TYPE,
                        ResolutionType::ALL,
                        value,
                        ResolutionType::Unknown,
                    );
                    analyzer.resolution_type = kind;
                    if !known {
                        unexpected = Some("resolution type");
                    }
                }
                "PSI:1000014" => analyzer.accuracy = finite(value, "accuracy")?,
                "PSI:1000015" => analyzer.scan_rate = finite(value, "scan rate")?,
                "PSI:1000016" => analyzer.scan_time = finite(value, "scan time")?,
                "PSI:1000018" => {
                    let (direction, known) = cv_enum(
                        SCAN_DIRECTION,
                        ScanDirection::ALL,
                        value,
                        ScanDirection::Unknown,
                    );
                    analyzer.scan_direction = direction;
                    if !known {
                        unexpected = Some("scan direction");
                    }
                }
                "PSI:1000019" => {
                    let (law, known) = cv_enum(SCAN_LAW, ScanLaw::ALL, value, ScanLaw::Unknown);
                    analyzer.scan_law = law;
                    if !known {
                        unexpected = Some("scan law");
                    }
                }
                "PSI:1000021" => {
                    let (state, known) = cv_enum(
                        REFLECTRON_STATE,
                        ReflectronState::ALL,
                        value,
                        ReflectronState::Unknown,
                    );
                    analyzer.reflectron_state = state;
                    if !known {
                        unexpected = Some("reflectron state");
                    }
                }
                "PSI:1000022" => {
                    analyzer.tof_total_path_length = finite(value, "TOF total path length")?;
                }
                "PSI:1000023" => analyzer.isolation_width = finite(value, "isolation width")?,
                "PSI:1000024" => {
                    analyzer.final_ms_exponent = value
                        .parse()
                        .map_err(|_| parse_error(format!("Int conversion error of \"{value}\"")))?;
                }
                "PSI:1000025" => {
                    analyzer.magnetic_field_strength = finite(value, "magnetic field strength")?;
                }
                // PSI:1000017 ScanFunction and PSI:1000020
                // TandemScanningMethod are read and deliberately dropped;
                // `cv_terms_[4]` and `cv_terms_[12]` are the two tables
                // `init_` leaves empty.
                "PSI:1000017" | "PSI:1000020" => {}
                _ => return Ok(Some("AnalyzerList.Analyzer.UserParam")),
            }
        }
        if let Some(what) = unexpected {
            self.warn(format!("Unexpected CV entry '{what}'='{value}'"));
        }
        Ok(None)
    }

    fn cv_additional(&mut self, accession: &str, value: &str) -> Result<Option<&'static str>> {
        let instrument = &mut self.experiment.settings.instrument;
        match accession {
            "PSI:1000030" => instrument.vendor = value.to_owned(),
            "PSI:1000031" => instrument.model = value.to_owned(),
            "PSI:1000032" => instrument.customizations = value.to_owned(),
            _ => return Ok(Some("Description.Instrument.Additional")),
        }
        Ok(None)
    }

    fn cv_processing(&mut self, accession: &str) -> Result<Option<&'static str>> {
        let processing = self.processing_mut()?;
        match accession {
            "PSI:1000033" => {
                processing.actions.insert(ProcessingAction::Deisotoping);
            }
            "PSI:1000034" => {
                processing
                    .actions
                    .insert(ProcessingAction::ChargeDeconvolution);
            }
            "PSI:1000127" => {
                processing.actions.insert(ProcessingAction::PeakPicking);
            }
            // PSI:1000035 PeakProcessing is read and deliberately dropped.
            "PSI:1000035" => {}
            _ => return Ok(Some("DataProcessing.DataProcessing.UserParam")),
        }
        Ok(None)
    }

    // -- binary arrays ------------------------------------------------------

    /// `MzDataHandler::fillData_` (`MzDataHandler.cpp:487-577`).
    ///
    /// The source reads `precisions_[0]` and `precisions_[1]` before its
    /// `data_to_decode_.size() < 2` guard, so a spectrum with fewer than two
    /// binary arrays reads past the end of a `std::vector`; it then reports an
    /// m/z-versus-intensity length disagreement with the non-fatal
    /// `error(LOAD, …)` and indexes every other array at every position below
    /// the m/z length regardless. Both are refused here.
    fn fill_data(&mut self) -> Result<()> {
        let arrays = std::mem::take(&mut self.arrays);
        if arrays.is_empty() {
            // A `<spectrum>` with no binary array at all carries no peaks; the
            // source reads `precisions_[0]` past the end of an empty vector
            // here. An annotation array declared without any `<data>` at all
            // is still an error, because the source would then read past the
            // end of its decoded lists instead.
            if !self.spectrum.float_data_arrays.is_empty() {
                return Err(parse_error(
                    "mzData spectrum declares <supDataArrayBinary> but no <data> at all",
                ));
            }
            return Ok(());
        }
        let mut mz: Option<&Encoded> = None;
        let mut intensity_array: Option<&Encoded> = None;
        let mut supplemental = Vec::new();
        for array in &arrays {
            match array.kind {
                ArrayKind::Mz if mz.is_none() => mz = Some(array),
                ArrayKind::Intensity if intensity_array.is_none() => {
                    intensity_array = Some(array);
                }
                ArrayKind::Mz | ArrayKind::Intensity => {
                    return Err(parse_error(
                        "mzData spectrum declares a second m/z or intensity array",
                    ));
                }
                // Each `<supDataArrayBinary>` must hold exactly one `<data>`,
                // so the payloads and the float data arrays pair up one to one
                // and in order. The source pairs them by position alone.
                ArrayKind::Supplemental(slot) => {
                    if slot != supplemental.len() {
                        return Err(parse_error(
                            "mzData <supDataArrayBinary> does not hold exactly one <data> child",
                        ));
                    }
                    supplemental.push(array);
                }
            }
        }
        let (Some(mz), Some(intensities)) = (mz, intensity_array) else {
            return Err(parse_error(
                "mzData spectrum has an m/z array without an intensity array or the reverse",
            ));
        };
        let positions = self.decode(mz, "m/z array")?;
        let magnitudes = self.decode(intensities, "intensity array")?;
        if positions.len() != magnitudes.len() {
            return Err(parse_error(format!(
                "Length of data array for m/z differs from length of intensity data: {} vs. {} .",
                positions.len(),
                magnitudes.len()
            )));
        }
        if let Some(declared) = mz.declared {
            if declared != positions.len() {
                self.warn(format!(
                    "Length of data arrays (m/z and int) differs from value in attribute 'length': {} vs. {}.",
                    positions.len(),
                    declared
                ));
            }
        }
        if supplemental.len() != self.spectrum.float_data_arrays.len() {
            return Err(parse_error(
                "mzData spectrum has a <supDataArrayBinary> without a <data> child or the reverse",
            ));
        }
        let mut annotations = Vec::new();
        annotations
            .try_reserve_exact(supplemental.len())
            .map_err(|_| limit("cannot allocate mzData annotation arrays"))?;
        for array in &supplemental {
            let values = self.decode(array, "supplemental array")?;
            if values.len() != positions.len() {
                return Err(parse_error(format!(
                    "Length of meta data array differs from spectrum length. meta data array: {} / spectrum: {} .",
                    values.len(),
                    positions.len()
                )));
            }
            annotations.push(values);
        }
        // Selection happens before anything is appended to the spectrum, so a
        // refusal leaves it untouched.
        let mz_range = self.options.has_mz_range().then(|| self.options.mz_range());
        let intensity_range = self
            .options
            .has_intensity_range()
            .then(|| self.options.intensity_range());
        let mut keep = Vec::new();
        keep.try_reserve_exact(positions.len())
            .map_err(|_| limit("cannot allocate the mzData peak selection"))?;
        for (index, (&position, &magnitude)) in positions.iter().zip(&magnitudes).enumerate() {
            if mz_range.is_none_or(|range| encloses(range, position))
                && intensity_range.is_none_or(|range| encloses(range, magnitude))
            {
                keep.push(index);
            }
        }
        self.report.peaks_filtered_out += positions.len() - keep.len();
        self.peak_budget = self
            .peak_budget
            .checked_sub(keep.len())
            .ok_or_else(|| limit("mzData exceeds the configured total peak ceiling"))?;
        let mut peaks = Vec::new();
        peaks
            .try_reserve_exact(keep.len())
            .map_err(|_| limit("cannot allocate mzData peaks"))?;
        for &index in &keep {
            let mz = positions[index];
            if !mz.is_finite() {
                // The source stores whatever the bytes decode to; a nonfinite
                // coordinate makes every later sort and binary search
                // undefined, so it is refused here as an attribute value
                // would be.
                return Err(parse_error("mzData m/z array holds a nonfinite value"));
            }
            peaks.push(Peak1D::new(mz, intensity(magnitudes[index])?));
        }
        let mut selected = Vec::new();
        selected
            .try_reserve_exact(annotations.len())
            .map_err(|_| limit("cannot allocate mzData annotation arrays"))?;
        for values in &annotations {
            let mut column = Vec::new();
            column
                .try_reserve_exact(keep.len())
                .map_err(|_| limit("cannot allocate an mzData annotation array"))?;
            for &index in &keep {
                column.push(intensity(values[index])?);
            }
            selected.push(column);
        }
        self.spectrum.peaks = peaks;
        for (array, values) in self.spectrum.float_data_arrays.iter_mut().zip(selected) {
            array.data = values;
        }
        Ok(())
    }

    /// Decode one array's payload as `count` f64 values, honouring its declared
    /// precision and byte order.
    ///
    /// `Base64::decodeUncompressed_` (`Base64.h:309-345`) returns nothing for a
    /// payload shorter than four characters, refuses a length that is not a
    /// multiple of four, and then assigns `s.size() / element_size` elements —
    /// so a payload holding a trailing partial element silently loses it. The
    /// upstream fixtures `MzDataFile_3_minimal.mzData` and
    /// `MzDataFile_4_64bit.mzData` both rely on that: their 36-character m/z
    /// payload decodes to 25 bytes and yields three doubles, with one stray
    /// byte dropped. This port therefore truncates to whole elements too, and
    /// warns.
    fn decode(&mut self, array: &Encoded, what: &str) -> Result<Vec<f64>> {
        let width = array.precision.width();
        if array.payload.len() < 4 {
            return Ok(Vec::new());
        }
        if array.payload.len() % 4 != 0 {
            return Err(parse_error(format!(
                "Malformed base64 {what}, length is not a multiple of 4."
            )));
        }
        if array.payload.len() / 4 * 3 > self.limits.max_array_bytes {
            return Err(limit(format!(
                "mzData {what} exceeds the configured byte ceiling"
            )));
        }
        let raw = STANDARD
            .decode(array.payload.as_bytes())
            .map_err(|e| parse_error(format!("invalid base64 in mzData {what}: {e}")))?;
        let count = raw.len() / width;
        if count > self.limits.max_array_elements {
            return Err(limit(format!(
                "mzData {what} holds {count} elements, above the configured ceiling"
            )));
        }
        if raw.len() % width != 0 {
            self.warn(format!(
                "mzData {what} holds {} bytes, which is not a whole number of {width}-byte elements; the trailing bytes are dropped",
                raw.len()
            ));
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| limit(format!("cannot allocate the mzData {what}")))?;
        for chunk in raw.chunks_exact(width) {
            values.push(match (array.precision, array.endian) {
                (Precision::Bits32, Endian::Little) => f64::from(f32::from_le_bytes(four(chunk))),
                (Precision::Bits32, Endian::Big) => f64::from(f32::from_be_bytes(four(chunk))),
                (Precision::Bits64, Endian::Little) => f64::from_le_bytes(eight(chunk)),
                (Precision::Bits64, Endian::Big) => f64::from_be_bytes(eight(chunk)),
            });
        }
        Ok(values)
    }
}

fn local_name(element: &BytesStart<'_>) -> Result<String> {
    // mzData has no XML namespace at all, so the qualified name is the name.
    std::str::from_utf8(element.name().as_ref())
        .map_err(|e| parse_error(e.to_string()))
        .map(str::to_owned)
}

fn local_name_end(element: &BytesEnd<'_>) -> Result<String> {
    std::str::from_utf8(element.name().as_ref())
        .map_err(|e| parse_error(e.to_string()))
        .map(str::to_owned)
}

/// The value of one attribute, with XML character references resolved.
///
/// quick-xml's duplicate-attribute check is left on: a repeated attribute is
/// not well-formed XML and Xerces refuses it, so the source parser never sees
/// one either.
fn attribute(element: &BytesStart<'_>, name: &str) -> Result<Option<String>> {
    for attribute in element.attributes() {
        let attribute = attribute.map_err(|e| parse_error(e.to_string()))?;
        if attribute.key.as_ref() == name.as_bytes() {
            let value = attribute
                .unescape_value()
                .map_err(|e| parse_error(e.to_string()))?;
            return Ok(Some(value.into_owned()));
        }
    }
    Ok(None)
}

/// `attributeAsString_` (`XMLHandler.h:385-390`), which calls `fatalError` and
/// so throws `Exception::ParseError` when the attribute is absent.
fn required(element: &BytesStart<'_>, name: &str) -> Result<String> {
    attribute(element, name)?
        .ok_or_else(|| parse_error(format!("Required attribute '{name}' not present!")))
}

/// `optionalAttributeAsString_`: absent is not an error.
fn optional(element: &BytesStart<'_>, name: &str) -> Result<Option<String>> {
    attribute(element, name)
}

/// The source `asDouble_` logs a non-fatal error and substitutes `0.0` for an
/// unparsable value (`XMLHandler.h:305-317`). Substituting zero for a mass or
/// an intensity is worse than refusing the document, so this returns an error;
/// nonfinite values are refused for the same reason.
fn finite(value: &str, what: &str) -> Result<f64> {
    value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|number| number.is_finite())
        .ok_or_else(|| parse_error(format!("Double conversion error of \"{value}\" ({what})")))
}

fn intensity(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(parse_error("mzData intensity overflows f32"))
    }
}

fn four(chunk: &[u8]) -> [u8; 4] {
    [chunk[0], chunk[1], chunk[2], chunk[3]]
}

fn eight(chunk: &[u8]) -> [u8; 8] {
    [
        chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
    ]
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Scan modes that survive a store followed by a load.
///
/// `writeTo` maps `MS1SPECTRUM` and `MSNSPECTRUM` onto the same `MassScan`
/// value as `MASSSPECTRUM` (`MzDataHandler.cpp:868-881`), and emits
/// `PhotodiodeArrayDetector`, `EnhancedMultiplyChargedScan` and
/// `TimeDelayedFragmentationScan` for three more (`:899-915`) — none of which
/// its own `cvParam_` recognises, so each reads back as `MASSSPECTRUM` with a
/// warning. Everything not in this list is therefore lossy.
const ROUND_TRIP_SCAN_MODES: &[ScanMode] = &[
    ScanMode::Unknown,
    ScanMode::MassSpectrum,
    ScanMode::SelectedIonMonitoring,
    ScanMode::SelectedReactionMonitoring,
    ScanMode::ConsecutiveReactionMonitoring,
    ScanMode::ConstantNeutralGain,
    ScanMode::ConstantNeutralLoss,
    ScanMode::Precursor,
];

/// Processing actions mzData has a `cvParam` for (`MzDataHandler.cpp:739-751`).
const WRITABLE_ACTIONS: &[ProcessingAction] = &[
    ProcessingAction::Deisotoping,
    ProcessingAction::ChargeDeconvolution,
    ProcessingAction::PeakPicking,
];

fn unsupported(what: &str) -> Error {
    Error::Unsupported(format!(
        "mzData cannot store {what}; set WriteOptions::discard_unrepresentable to drop it"
    ))
}

/// Refuse every record mzData has no element for, before anything is written.
fn preflight_store(
    experiment: &MSExperiment,
    options: &WriteOptions,
    report: &mut StoreReport,
) -> Result<()> {
    let settings = &experiment.settings;
    if !options.discard_unrepresentable {
        if !experiment.chromatograms.is_empty() {
            return Err(unsupported("chromatograms"));
        }
        if experiment.sql_run_id != 0 {
            return Err(unsupported("an SQL run identifier"));
        }
        if settings.source_files.len() > 1 {
            return Err(unsupported("more than one source file"));
        }
        for file in &settings.source_files {
            if !file.checksum.is_empty()
                || file.size_mb != 0.0
                || !file.native_id_type.is_empty()
                || !file.native_id_type_accession.is_empty()
                || !file.cv_terms.terms().is_empty()
                || !file.cv_terms.metadata.is_empty()
            {
                return Err(unsupported(
                    "a source file checksum, size, native-ID type or CV terms",
                ));
            }
        }
        for contact in &settings.contacts {
            if !contact.email.is_empty()
                || !contact.url.is_empty()
                || !contact.address.is_empty()
                || !contact.metadata.is_empty()
            {
                return Err(unsupported(
                    "a contact email, URL, address or contact metadata",
                ));
            }
            if contact.first_name.contains(' ') || contact.last_name.contains(' ') {
                return Err(unsupported(
                    "a contact name whose parts contain spaces, because <name> is a single field",
                ));
            }
        }
        let sample = &settings.sample;
        if !sample.organism.is_empty()
            || !sample.comment.is_empty()
            || !sample.subsamples.is_empty()
        {
            return Err(unsupported("a sample organism, comment or subsample"));
        }
        let instrument = &settings.instrument;
        if instrument.ion_sources.len() > 1 {
            return Err(unsupported("more than one ion source"));
        }
        if instrument.ion_detectors.len() > 1 {
            return Err(unsupported("more than one ion detector"));
        }
        if instrument.ion_optics != crate::metadata::IonOpticsType::Unknown {
            return Err(unsupported("an instrument ion-optics type"));
        }
        if instrument.software != crate::metadata::Software::default() {
            return Err(unsupported("instrument software"));
        }
        for source in &instrument.ion_sources {
            if source.order != 0 {
                return Err(unsupported("an ion source order"));
            }
        }
        for detector in &instrument.ion_detectors {
            if detector.order != 0 {
                return Err(unsupported("an ion detector order"));
            }
        }
        for analyzer in &instrument.mass_analyzers {
            if analyzer.order != 0 {
                return Err(unsupported("a mass analyzer order"));
            }
        }
        if settings.hplc != crate::metadata::HPLC::default() {
            return Err(unsupported("HPLC settings"));
        }
        if settings.date_time != crate::data_structures::DateTime::default() {
            return Err(unsupported("an experiment date"));
        }
        if !settings.comment.is_empty() {
            return Err(unsupported("an experiment comment"));
        }
        if !settings.fraction_identifier.is_empty() {
            return Err(unsupported("a fraction identifier"));
        }
        if !settings.instrument_configurations.is_empty() {
            return Err(unsupported("named instrument configurations"));
        }
        if !settings.metadata.is_empty() {
            return Err(unsupported("experiment-level metadata"));
        }
        if !options.write_supplemental_data
            && experiment
                .spectra
                .iter()
                .any(|spectrum| !spectrum.float_data_arrays.is_empty())
        {
            return Err(unsupported(
                "float data arrays while write_supplemental_data is false",
            ));
        }
        let shared = experiment
            .spectra
            .first()
            .and_then(|spectrum| spectrum.data_processing.first());
        for spectrum in &experiment.spectra {
            preflight_spectrum(spectrum, shared)?;
        }
    }
    // The renumbering verdict is part of the report whether or not discarding
    // is allowed, because refusing it is not an option: mzData's `id` is an
    // integer and a native ID that is not one has nowhere to go.
    let (all_numbers, all_empty, _) = native_id_shape(experiment);
    report.renumbered = !all_numbers;
    if !all_numbers && !all_empty {
        if !options.discard_unrepresentable {
            return Err(unsupported(
                "native IDs that are neither numbers nor 'spectrum=' followed by a number",
            ));
        }
        report.warning_count += 1;
        report.warnings.push(
            "Not all spectrum native IDs are numbers or correctly prefixed with 'spectrum='. The spectra are renumbered and the native IDs are lost!"
                .into(),
        );
    }
    Ok(())
}

fn preflight_spectrum(spectrum: &MSSpectrum, shared: Option<&Arc<DataProcessing>>) -> Result<()> {
    if !spectrum.name.is_empty() {
        return Err(unsupported("a spectrum name"));
    }
    if spectrum.source_file != SourceFile::default() {
        return Err(unsupported("a per-spectrum source file"));
    }
    if !spectrum.products.is_empty() {
        return Err(unsupported("spectrum products"));
    }
    if !spectrum.peptide_identifications.is_empty() {
        return Err(unsupported("peptide identifications"));
    }
    if !spectrum.metadata.is_empty() {
        return Err(unsupported("spectrum-level metadata"));
    }
    if spectrum.drift_time != -1.0 {
        return Err(unsupported("a spectrum drift time"));
    }
    if !spectrum.integer_data_arrays.is_empty() || !spectrum.string_data_arrays.is_empty() {
        return Err(unsupported(
            "integer or string data arrays; only <supDataArrayBinary> float arrays exist",
        ));
    }
    if spectrum.instrument_settings.scan_windows.len() > 1 {
        return Err(unsupported("more than one scan window per scan"));
    }
    let settings = &spectrum.instrument_settings;
    if !ROUND_TRIP_SCAN_MODES.contains(&settings.scan_mode) {
        return Err(unsupported(
            "this scan mode, which mzData cannot round-trip",
        ));
    }
    if settings.zoom_scan && settings.scan_mode != ScanMode::MassSpectrum {
        return Err(unsupported(
            "a zoom scan outside the MassSpectrum scan mode, because the flag is written as the scan mode itself",
        ));
    }
    if spectrum.spectrum_type == SpectrumType::Unknown
        && !spectrum.acquisition_info.acquisitions.is_empty()
    {
        return Err(unsupported(
            "an unknown spectrum type together with acquisitions, which is written as 'discrete'",
        ));
    }
    if spectrum.spectrum_type != SpectrumType::Unknown
        && spectrum.acquisition_info.acquisitions.is_empty()
    {
        return Err(unsupported(
            "a spectrum type without acquisitions, because <acqSpecification> carries it",
        ));
    }
    if !spectrum.acquisition_info.metadata.is_empty() {
        return Err(unsupported("acquisition-info metadata"));
    }
    for acquisition in &spectrum.acquisition_info.acquisitions {
        if !acquisition.identifier.is_empty() && acquisition.identifier.parse::<i32>().is_err() {
            return Err(unsupported(
                "an acquisition identifier that is not a 32-bit integer",
            ));
        }
    }
    for precursor in &spectrum.precursors {
        if precursor.activation_methods.len() > 1 {
            return Err(unsupported("more than one activation method"));
        }
        if precursor.isolation_window_lower_offset != 0.0
            || precursor.isolation_window_upper_offset != 0.0
            || precursor.isolation_target_mz.is_some()
            || precursor.drift_time.is_some()
            || precursor.drift_window_lower_offset != 0.0
            || precursor.drift_window_upper_offset != 0.0
            || !precursor.possible_charge_states.is_empty()
            || !precursor.cv_terms.terms().is_empty()
            || precursor.spectrum_reference.is_some()
        {
            return Err(unsupported(
                "a precursor isolation window, charge-state list, CV term, drift time or spectrum reference",
            ));
        }
    }
    for array in &spectrum.float_data_arrays {
        if !array.data_processing.is_empty() {
            return Err(unsupported("data-array processing history"));
        }
    }
    // `writeTo` takes the *first* data processing of the *first* spectrum and
    // applies it to the whole document (`MzDataHandler.cpp:716-755`).
    if spectrum.data_processing.len() > 1 {
        return Err(unsupported("more than one data-processing record"));
    }
    match (shared, spectrum.data_processing.first()) {
        (Some(shared), Some(own)) if **shared != **own => {
            return Err(unsupported(
                "spectra with different data-processing records, because mzData has one <dataProcessing> per document",
            ));
        }
        (Some(_), None) | (None, Some(_)) => {
            return Err(unsupported(
                "some spectra with and some without a data-processing record",
            ));
        }
        _ => {}
    }
    if let Some(processing) = spectrum.data_processing.first() {
        if !processing.software.cv_terms.metadata.is_empty()
            || !processing.software.cv_terms.terms().is_empty()
        {
            return Err(unsupported(
                "software metadata, because <software> has only <name> and <version>",
            ));
        }
        if processing
            .actions
            .iter()
            .any(|action| !WRITABLE_ACTIONS.contains(action))
        {
            return Err(unsupported(
                "a processing action outside Deisotoping, ChargeDeconvolution and PeakPicking",
            ));
        }
    }
    Ok(())
}

/// `MzDataHandler.cpp:764-800`: whether every native ID is a number, whether
/// every one that is not is empty, and whether all carry the `spectrum=`
/// prefix.
fn native_id_shape(experiment: &MSExperiment) -> (bool, bool, bool) {
    let mut all_numbers = true;
    let mut all_empty = true;
    let mut all_prefixed = true;
    for spectrum in &experiment.spectra {
        let mut native_id = spectrum.native_id.as_str();
        match native_id.strip_prefix("spectrum=") {
            Some(rest) => native_id = rest,
            None => all_prefixed = false,
        }
        if native_id.parse::<i32>().is_err() {
            all_numbers = false;
            all_prefixed = false;
            if !native_id.is_empty() {
                all_empty = false;
            }
        }
    }
    (all_numbers, all_empty, all_prefixed)
}

/// XML-escape and keep the declared ISO-8859-1 encoding truthful.
///
/// `writeTo` streams every `std::string` straight out, so an `&` or a `<` in
/// any metadata value produces a document that is not well-formed XML, and a
/// UTF-8 payload produces bytes the declared encoding cannot describe. Here
/// the five XML delimiters are escaped and any character above `U+00FF` is
/// written as a numeric character reference, which is valid in an ISO-8859-1
/// document and reads back identically.
fn escape(value: &str) -> String {
    let escaped = quick_xml::escape::escape(value);
    if escaped.chars().all(|c| (c as u32) <= 0xff) {
        return escaped.into_owned();
    }
    let mut out = String::with_capacity(escaped.len());
    for c in escaped.chars() {
        if (c as u32) <= 0xff {
            out.push(c);
        } else {
            out.push_str(&format!("&#x{:X};", c as u32));
        }
    }
    out
}

/// ISO-8859-1 bytes for an already-escaped string: every character is at most
/// `U+00FF` after [`escape`], so this is the exact inverse of the reader's
/// transcode.
fn latin1(value: &str) -> Vec<u8> {
    value.chars().map(|c| c as u8).collect()
}

struct Sink<'a> {
    writer: &'a mut dyn Write,
}

impl Sink<'_> {
    fn raw(&mut self, text: &str) -> Result<()> {
        self.writer.write_all(text.as_bytes())?;
        Ok(())
    }

    /// Write escaped, ISO-8859-1 encoded text.
    fn text(&mut self, value: &str) -> Result<()> {
        let escaped = escape(value);
        self.writer.write_all(&latin1(&escaped))?;
        Ok(())
    }

    fn tabs(&mut self, depth: usize) -> Result<()> {
        for _ in 0..depth {
            self.writer.write_all(b"\t")?;
        }
        Ok(())
    }

    /// `writeCVS_` for a string value: nothing is written for an empty value.
    fn cv_string(&mut self, depth: usize, accession: &str, name: &str, value: &str) -> Result<()> {
        if value.is_empty() {
            return Ok(());
        }
        self.tabs(depth)?;
        self.raw("<cvParam cvLabel=\"psi\" accession=\"PSI:")?;
        self.raw(accession)?;
        self.raw("\" name=\"")?;
        self.raw(name)?;
        self.raw("\" value=\"")?;
        self.text(value)?;
        self.raw("\"/>\n")
    }

    /// `writeCVS_` for a numeric value: nothing is written for exactly zero.
    fn cv_number(&mut self, depth: usize, accession: &str, name: &str, value: f64) -> Result<()> {
        if value == 0.0 {
            return Ok(());
        }
        self.cv_string(depth, accession, name, &stream_float(value)?)
    }

    /// `writeCVS_` for an enumeration: nothing is written when the table has no
    /// term for the value, which includes every enum's zero value.
    fn cv_enum<T: PartialEq>(
        &mut self,
        depth: usize,
        table: &[&str],
        all: &[T],
        value: &T,
        accession: &str,
        name: &str,
    ) -> Result<()> {
        match cv_term(table, all, value) {
            Some(term) => self.cv_string(depth, accession, name, term),
            None => Ok(()),
        }
    }

    /// `writeUserParam_`: every metadata key whose first character is not `#`.
    fn user_params(&mut self, depth: usize, metadata: &MetaInfo) -> Result<()> {
        for (name, value) in metadata {
            if name.starts_with('#') {
                continue;
            }
            self.tabs(depth)?;
            self.raw("<userParam name=\"")?;
            self.text(name)?;
            self.raw("\" value=\"")?;
            self.text(&meta_text(value)?)?;
            self.raw("\"/>\n")?;
        }
        Ok(())
    }

    /// `writeBinary_`: one `<mzArrayBinary>`, `<intenArrayBinary>` or
    /// `<supDataArrayBinary>` element with its base64 payload.
    ///
    /// `name` carries the array name and its one-based id, and is `None` for
    /// the two peak arrays — the source states the same condition as "the
    /// `name` and `id` are only used if the `tag` is *supDataArrayBinary* or
    /// *supDataArray*", and only emits the `id` attribute and the
    /// `<arrayName>` child for those two tags. `supDataArray` itself is a tag
    /// `writeTo` never passes.
    fn binary(
        &mut self,
        tag: &str,
        precision: Precision,
        values: &[f64],
        name: Option<(&str, usize)>,
    ) -> Result<()> {
        self.raw("\t\t\t<")?;
        self.raw(tag)?;
        if let Some((_, id)) = name {
            self.raw(&format!(" id=\"{id}\""))?;
        }
        self.raw(">\n")?;
        if let Some((array_name, _)) = name {
            self.raw("\t\t\t\t<arrayName>")?;
            self.text(array_name)?;
            self.raw("</arrayName>\n")?;
        }
        let mut raw = Vec::new();
        raw.try_reserve_exact(values.len().saturating_mul(precision.width()))
            .map_err(|_| limit("cannot allocate the mzData binary array"))?;
        for &value in values {
            match precision {
                Precision::Bits32 => raw.extend_from_slice(&(value as f32).to_le_bytes()),
                Precision::Bits64 => raw.extend_from_slice(&value.to_le_bytes()),
            }
        }
        self.raw(&format!(
            "\t\t\t\t<data precision=\"{}\" endian=\"little\" length=\"{}\">",
            precision.name(),
            values.len()
        ))?;
        self.raw(&STANDARD.encode(&raw))?;
        self.raw("</data>\n\t\t\t</")?;
        self.raw(tag)?;
        self.raw(">\n")
    }
}

fn write_document(
    writer: &mut dyn Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
    report: &mut StoreReport,
) -> Result<()> {
    let mut out = Sink { writer };
    let settings = &experiment.settings;
    out.raw("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n")?;
    out.raw("<mzData version=\"1.05\" accessionNumber=\"")?;
    out.text(&settings.document.identifier)?;
    out.raw(&format!(
        "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"{SCHEMA_LOCATION}\">\n"
    ))?;
    write_description(&mut out, experiment, report)?;
    write_spectra(&mut out, experiment, options, report)?;
    Ok(())
}

fn write_description(
    out: &mut Sink<'_>,
    experiment: &MSExperiment,
    report: &mut StoreReport,
) -> Result<()> {
    let settings = &experiment.settings;
    let sample = &settings.sample;
    out.raw("\t<description>\n\t\t<admin>\n\t\t\t<sampleName>")?;
    out.text(&sample.name)?;
    out.raw("</sampleName>\n")?;
    if !sample.number.is_empty()
        || sample.state != SampleState::Unknown
        || sample.mass != 0.0
        || sample.volume != 0.0
        || sample.concentration != 0.0
        || !sample.metadata.is_empty()
    {
        out.raw("\t\t\t<sampleDescription>\n")?;
        out.cv_string(4, "1000001", "SampleNumber", &sample.number)?;
        out.cv_enum(
            4,
            SAMPLE_STATE,
            SampleState::ALL,
            &sample.state,
            "1000003",
            "SampleState",
        )?;
        out.cv_number(4, "1000004", "SampleMass", sample.mass)?;
        out.cv_number(4, "1000005", "SampleVolume", sample.volume)?;
        out.cv_number(4, "1000006", "SampleConcentration", sample.concentration)?;
        out.user_params(4, &sample.metadata)?;
        out.raw("\t\t\t</sampleDescription>\n")?;
    }
    if let Some(file) = settings.source_files.first() {
        out.raw("\t\t\t<sourceFile>\n\t\t\t\t<nameOfFile>")?;
        out.text(&file.name)?;
        out.raw("</nameOfFile>\n\t\t\t\t<pathToFile>")?;
        out.text(&file.path)?;
        out.raw("</pathToFile>\n")?;
        if !file.file_type.is_empty() {
            out.raw("\t\t\t\t<fileType>")?;
            out.text(&file.file_type)?;
            out.raw("</fileType>\n")?;
        }
        out.raw("\t\t\t</sourceFile>\n")?;
    }
    if settings.source_files.len() > 1 {
        warn(
            report,
            "The MzData format can store only one source file. Only the first one is stored!",
        );
    }
    for contact in &settings.contacts {
        out.raw("\t\t\t<contact>\n\t\t\t\t<name>")?;
        out.text(&contact.first_name)?;
        out.raw(" ")?;
        out.text(&contact.last_name)?;
        out.raw("</name>\n\t\t\t\t<institution>")?;
        out.text(&contact.institution)?;
        out.raw("</institution>\n")?;
        if !contact.contact_info.is_empty() {
            out.raw("\t\t\t\t<contactInfo>")?;
            out.text(&contact.contact_info)?;
            out.raw("</contactInfo>\n")?;
        }
        out.raw("\t\t\t</contact>\n")?;
    }
    if settings.contacts.is_empty() {
        // mzData requires a contact entry, so an empty one is emitted.
        out.raw("\t\t\t<contact>\n\t\t\t\t<name></name>\n\t\t\t\t<institution></institution>\n\t\t\t</contact>\n")?;
    }
    out.raw("\t\t</admin>\n")?;
    let instrument = &settings.instrument;
    out.raw("\t\t<instrument>\n\t\t\t<instrumentName>")?;
    out.text(&instrument.name)?;
    out.raw("</instrumentName>\n\t\t\t<source>\n")?;
    if let Some(source) = instrument.ion_sources.first() {
        out.cv_enum(
            4,
            INLET_TYPE,
            InletType::ALL,
            &source.inlet_type,
            "1000007",
            "InletType",
        )?;
        out.cv_enum(
            4,
            IONIZATION_TYPE,
            IonizationMethod::ALL,
            &source.ionization_method,
            "1000008",
            "IonizationType",
        )?;
        out.cv_enum(
            4,
            IONIZATION_MODE,
            Polarity::ALL,
            &source.polarity,
            "1000009",
            "IonizationMode",
        )?;
        out.user_params(4, &source.metadata)?;
    }
    if instrument.ion_sources.len() > 1 {
        warn(
            report,
            "The MzData format can store only one ion source. Only the first one is stored!",
        );
    }
    out.raw("\t\t\t</source>\n")?;
    if instrument.mass_analyzers.is_empty() {
        out.raw("\t\t\t<analyzerList count=\"1\">\n\t\t\t\t<analyzer>\n\t\t\t\t</analyzer>\n")?;
    } else {
        out.raw(&format!(
            "\t\t\t<analyzerList count=\"{}\">\n",
            instrument.mass_analyzers.len()
        ))?;
        for analyzer in &instrument.mass_analyzers {
            out.raw("\t\t\t\t<analyzer>\n")?;
            out.cv_enum(
                5,
                ANALYZER_TYPE,
                AnalyzerType::ALL,
                &analyzer.analyzer_type,
                "1000010",
                "AnalyzerType",
            )?;
            out.cv_number(5, "1000011", "MassResolution", analyzer.resolution)?;
            out.cv_enum(
                5,
                RESOLUTION_METHOD,
                ResolutionMethod::ALL,
                &analyzer.resolution_method,
                "1000012",
                "ResolutionMethod",
            )?;
            out.cv_enum(
                5,
                RESOLUTION_TYPE,
                ResolutionType::ALL,
                &analyzer.resolution_type,
                "1000013",
                "ResolutionType",
            )?;
            out.cv_number(5, "1000014", "Accuracy", analyzer.accuracy)?;
            out.cv_number(5, "1000015", "ScanRate", analyzer.scan_rate)?;
            out.cv_number(5, "1000016", "ScanTime", analyzer.scan_time)?;
            out.cv_enum(
                5,
                SCAN_DIRECTION,
                ScanDirection::ALL,
                &analyzer.scan_direction,
                "1000018",
                "ScanDirection",
            )?;
            out.cv_enum(
                5,
                SCAN_LAW,
                ScanLaw::ALL,
                &analyzer.scan_law,
                "1000019",
                "ScanLaw",
            )?;
            out.cv_enum(
                5,
                REFLECTRON_STATE,
                ReflectronState::ALL,
                &analyzer.reflectron_state,
                "1000021",
                "ReflectronState",
            )?;
            out.cv_number(
                5,
                "1000022",
                "TOFTotalPathLength",
                analyzer.tof_total_path_length,
            )?;
            out.cv_number(5, "1000023", "IsolationWidth", analyzer.isolation_width)?;
            out.cv_number(
                5,
                "1000024",
                "FinalMSExponent",
                f64::from(analyzer.final_ms_exponent),
            )?;
            out.cv_number(
                5,
                "1000025",
                "MagneticFieldStrength",
                analyzer.magnetic_field_strength,
            )?;
            out.user_params(5, &analyzer.metadata)?;
            out.raw("\t\t\t\t</analyzer>\n")?;
        }
    }
    out.raw("\t\t\t</analyzerList>\n\t\t\t<detector>\n")?;
    if let Some(detector) = instrument.ion_detectors.first() {
        out.cv_enum(
            4,
            DETECTOR_TYPE,
            DetectorType::ALL,
            &detector.detector_type,
            "1000026",
            "DetectorType",
        )?;
        out.cv_enum(
            4,
            ACQUISITION_MODE,
            DetectorAcquisitionMode::ALL,
            &detector.acquisition_mode,
            "1000027",
            "DetectorAcquisitionMode",
        )?;
        out.cv_number(4, "1000028", "DetectorResolution", detector.resolution)?;
        out.cv_number(
            4,
            "1000029",
            "SamplingFrequency",
            detector.adc_sampling_frequency,
        )?;
        out.user_params(4, &detector.metadata)?;
    }
    if instrument.ion_detectors.len() > 1 {
        warn(
            report,
            "The MzData format can store only one ion detector. Only the first one is stored!",
        );
    }
    out.raw("\t\t\t</detector>\n")?;
    if !instrument.vendor.is_empty()
        || !instrument.model.is_empty()
        || !instrument.customizations.is_empty()
    {
        out.raw("\t\t\t<additional>\n")?;
        out.cv_string(4, "1000030", "Vendor", &instrument.vendor)?;
        out.cv_string(4, "1000031", "Model", &instrument.model)?;
        out.cv_string(4, "1000032", "Customization", &instrument.customizations)?;
        out.user_params(4, &instrument.metadata)?;
        out.raw("\t\t\t</additional>\n")?;
    }
    out.raw("\t\t</instrument>\n")?;
    let processing = experiment
        .spectra
        .first()
        .and_then(|spectrum| spectrum.data_processing.first());
    match processing {
        None => out.raw(
            "\t\t<dataProcessing>\n\t\t\t<software>\n\t\t\t\t<name></name>\n\t\t\t\t<version></version>\n\t\t\t</software>\n\t\t</dataProcessing>\n",
        )?,
        Some(processing) => {
            out.raw("\t\t<dataProcessing>\n\t\t\t<software")?;
            if let Some(stamp) = processing.completion_time {
                out.raw(" completionTime=\"")?;
                out.text(&stamp.get().replace(' ', "T"))?;
                out.raw("\"")?;
            }
            out.raw(">\n\t\t\t\t<name>")?;
            out.text(&processing.software.name)?;
            out.raw("</name>\n\t\t\t\t<version>")?;
            out.text(&processing.software.version)?;
            out.raw("</version>\n\t\t\t</software>\n\t\t\t<processingMethod>\n")?;
            if processing.actions.contains(&ProcessingAction::Deisotoping) {
                out.raw(
                    "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"Deisotoping\" accession=\"PSI:1000033\" />\n",
                )?;
            }
            if processing
                .actions
                .contains(&ProcessingAction::ChargeDeconvolution)
            {
                out.raw(
                    "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"ChargeDeconvolution\" accession=\"PSI:1000034\" />\n",
                )?;
            }
            if processing.actions.contains(&ProcessingAction::PeakPicking) {
                out.raw(
                    "\t\t\t\t<cvParam cvLabel=\"psi\" name=\"Centroid Mass Spectrum\" accession=\"PSI:1000127\"/>\n",
                )?;
            }
            out.user_params(4, &processing.metadata)?;
            out.raw("\t\t\t</processingMethod>\n\t\t</dataProcessing>\n")?;
        }
    }
    out.raw("\t</description>\n")
}

fn write_spectra(
    out: &mut Sink<'_>,
    experiment: &MSExperiment,
    options: &WriteOptions,
    report: &mut StoreReport,
) -> Result<()> {
    if experiment.spectra.is_empty() {
        // `MzDataHandler.cpp:1057-1072`: an empty experiment still needs one
        // schema-valid spectrum, so a zero-length placeholder is emitted.
        out.raw("\t<spectrumList count=\"1\">\n\t\t<spectrum id=\"1\">\n")?;
        out.raw("\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n")?;
        out.raw("\t\t\t\t\t<spectrumInstrument msLevel=\"1\"/>\n")?;
        out.raw("\t\t\t\t</spectrumSettings>\n\t\t\t</spectrumDesc>\n")?;
        out.raw("\t\t\t<mzArrayBinary>\n")?;
        out.raw("\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n")?;
        out.raw("\t\t\t</mzArrayBinary>\n\t\t\t<intenArrayBinary>\n")?;
        out.raw("\t\t\t\t<data length=\"0\" endian=\"little\" precision=\"32\"></data>\n")?;
        out.raw("\t\t\t</intenArrayBinary>\n\t\t</spectrum>\n")?;
        return out.raw("\t</spectrumList>\n</mzData>\n");
    }
    let (all_numbers, _, all_prefixed) = native_id_shape(experiment);
    let mut level_id = BTreeMap::<u32, i64>::new();
    out.raw(&format!(
        "\t<spectrumList count=\"{}\">\n",
        experiment.spectra.len()
    ))?;
    for (index, spectrum) in experiment.spectra.iter().enumerate() {
        let position = i64::try_from(index)
            .ok()
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| limit("mzData spectrum index overflows"))?;
        let id = if all_prefixed {
            spectrum
                .native_id
                .strip_prefix("spectrum=")
                .and_then(|rest| rest.parse::<i32>().ok())
                .map_or(position, i64::from)
        } else if all_numbers {
            spectrum
                .native_id
                .parse::<i32>()
                .map_or(position, i64::from)
        } else {
            position
        };
        out.raw(&format!(
            "\t\t<spectrum id=\"{id}\">\n\t\t\t<spectrumDesc>\n\t\t\t\t<spectrumSettings>\n"
        ))?;
        if !spectrum.acquisition_info.acquisitions.is_empty() {
            out.raw("\t\t\t\t\t<acqSpecification spectrumType=\"")?;
            match spectrum.spectrum_type {
                SpectrumType::Centroid => out.raw("discrete")?,
                SpectrumType::Profile => out.raw("continuous")?,
                SpectrumType::Unknown => {
                    warn(report, "Spectrum type is unknown, assuming 'discrete'");
                    out.raw("discrete")?;
                }
            }
            out.raw("\" methodOfCombination=\"")?;
            out.text(&spectrum.acquisition_info.method_of_combination)?;
            out.raw(&format!(
                "\" count=\"{}\">\n",
                spectrum.acquisition_info.acquisitions.len()
            ))?;
            for acquisition in &spectrum.acquisition_info.acquisitions {
                let number = if acquisition.identifier.is_empty() {
                    0
                } else {
                    match acquisition.identifier.parse::<i32>() {
                        Ok(number) => number,
                        Err(_) => {
                            warn(
                                report,
                                &format!(
                                    "Could not convert acquisition identifier '{}' to an integer. Using '0' instead!",
                                    acquisition.identifier
                                ),
                            );
                            0
                        }
                    }
                };
                out.raw(&format!(
                    "\t\t\t\t\t\t<acquisition acqNumber=\"{number}\">\n"
                ))?;
                out.user_params(7, &acquisition.metadata)?;
                out.raw("\t\t\t\t\t\t</acquisition>\n")?;
            }
            out.raw("\t\t\t\t\t</acqSpecification>\n")?;
        }
        let iset = &spectrum.instrument_settings;
        out.raw(&format!(
            "\t\t\t\t\t<spectrumInstrument msLevel=\"{}\"",
            spectrum.ms_level
        ))?;
        level_id.insert(spectrum.ms_level, id);
        if let Some(window) = iset.scan_windows.first() {
            out.raw(&format!(
                " mzRangeStart=\"{}\" mzRangeStop=\"{}\"",
                stream_float(window.begin)?,
                stream_float(window.end)?
            ))?;
        }
        if iset.scan_windows.len() > 1 {
            warn(
                report,
                "The MzData format can store only one scan window for each scan. Only the first one is stored!",
            );
        }
        out.raw(">\n")?;
        let scan_mode = match iset.scan_mode {
            ScanMode::Unknown => None,
            ScanMode::MassSpectrum => Some(if iset.zoom_scan { "Zoom" } else { "MassScan" }),
            ScanMode::SelectedIonMonitoring => Some("SelectedIonDetection"),
            ScanMode::SelectedReactionMonitoring => Some("SelectedReactionMonitoring"),
            ScanMode::ConsecutiveReactionMonitoring => Some("ConsecutiveReactionMonitoring"),
            ScanMode::ConstantNeutralGain => Some("ConstantNeutralGainScan"),
            ScanMode::ConstantNeutralLoss => Some("ConstantNeutralLossScan"),
            ScanMode::Precursor => Some("PrecursorIonScan"),
            // Reached only with `discard_unrepresentable`; the preflight
            // refuses each of these otherwise.
            ScanMode::Ms1Spectrum | ScanMode::MsnSpectrum => Some("MassScan"),
            ScanMode::Absorption => Some("PhotodiodeArrayDetector"),
            ScanMode::EnhancedMultiplyCharged => Some("EnhancedMultiplyChargedScan"),
            ScanMode::TimeDelayedFragmentation => Some("TimeDelayedFragmentationScan"),
            other => {
                warn(
                    report,
                    &format!(
                        "Scan mode '{other}' not supported by mzData. Using 'MassScan' scan mode!"
                    ),
                );
                Some("MassScan")
            }
        };
        if let Some(mode) = scan_mode {
            out.raw(&format!(
                "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000036\" name=\"ScanMode\" value=\"{mode}\"/>\n"
            ))?;
        }
        match iset.polarity {
            Polarity::Positive => out.raw(
                "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000037\" name=\"Polarity\" value=\"Positive\"/>\n",
            )?,
            Polarity::Negative => out.raw(
                "\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000037\" name=\"Polarity\" value=\"Negative\"/>\n",
            )?,
            Polarity::Unknown => {}
        }
        out.cv_number(6, "1000039", "TimeInSeconds", spectrum.rt)?;
        out.user_params(6, &iset.metadata)?;
        out.raw("\t\t\t\t\t</spectrumInstrument>\n\t\t\t\t</spectrumSettings>\n")?;
        if !spectrum.precursors.is_empty() {
            let precursor_level = i64::from(spectrum.ms_level).saturating_sub(1);
            let reference = spectrum
                .ms_level
                .checked_sub(1)
                .and_then(|level| level_id.get(&level).copied())
                .unwrap_or(-1);
            out.raw(&format!(
                "\t\t\t\t<precursorList count=\"{}\">\n",
                spectrum.precursors.len()
            ))?;
            for precursor in &spectrum.precursors {
                out.raw(&format!(
                    "\t\t\t\t\t<precursor msLevel=\"{precursor_level}\" spectrumRef=\"{reference}\">\n"
                ))?;
                out.raw("\t\t\t\t\t\t<ionSelection>\n")?;
                let empty = *precursor == Precursor::default();
                if !empty {
                    out.cv_number(7, "1000040", "MassToChargeRatio", precursor.mz)?;
                    out.cv_number(7, "1000041", "ChargeState", f64::from(precursor.charge))?;
                    out.cv_number(7, "1000042", "Intensity", f64::from(precursor.intensity))?;
                    out.raw(
                        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000043\" name=\"IntensityUnit\" value=\"NumberOfCounts\"/>\n",
                    )?;
                    out.user_params(7, &precursor.cv_terms.metadata)?;
                }
                out.raw("\t\t\t\t\t\t</ionSelection>\n\t\t\t\t\t\t<activation>\n")?;
                if !empty {
                    if let Some(method) = precursor.activation_methods.iter().next() {
                        out.cv_enum(
                            7,
                            ACTIVATION_METHOD,
                            ActivationMethod::ALL,
                            method,
                            "1000044",
                            "ActivationMethod",
                        )?;
                    }
                    out.cv_number(7, "1000045", "CollisionEnergy", precursor.activation_energy)?;
                    out.raw(
                        "\t\t\t\t\t\t\t<cvParam cvLabel=\"psi\" accession=\"PSI:1000046\" name=\"EnergyUnit\" value=\"eV\"/>\n",
                    )?;
                }
                out.raw("\t\t\t\t\t\t</activation>\n\t\t\t\t\t</precursor>\n")?;
            }
            out.raw("\t\t\t\t</precursorList>\n")?;
        }
        out.raw("\t\t\t</spectrumDesc>\n")?;
        if options.write_supplemental_data {
            for (position, array) in spectrum.float_data_arrays.iter().enumerate() {
                out.raw(&format!(
                    "\t\t\t<supDesc supDataArrayRef=\"{}\">\n",
                    position + 1
                ))?;
                if !array.metadata.is_empty() {
                    out.raw("\t\t\t\t<supDataDesc>\n")?;
                    out.user_params(5, &array.metadata)?;
                    out.raw("\t\t\t\t</supDataDesc>\n")?;
                }
                out.raw("\t\t\t</supDesc>\n")?;
            }
        }
        let mz_precision = if options.mz_32_bit {
            Precision::Bits32
        } else {
            Precision::Bits64
        };
        let positions: Vec<f64> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
        out.binary("mzArrayBinary", mz_precision, &positions, None)?;
        let magnitudes: Vec<f64> = spectrum
            .peaks
            .iter()
            .map(|peak| f64::from(peak.intensity))
            .collect();
        out.binary("intenArrayBinary", Precision::Bits32, &magnitudes, None)?;
        if options.write_supplemental_data {
            for (position, array) in spectrum.float_data_arrays.iter().enumerate() {
                if array.data.len() != spectrum.peaks.len() {
                    warn(
                        report,
                        &format!(
                            "Length of meta data array (index:'{position}' name:'{}') differs from spectrum length. meta data array: {} / spectrum: {} .",
                            array.name,
                            array.data.len(),
                            spectrum.peaks.len()
                        ),
                    );
                }
                let values: Vec<f64> = array.data.iter().map(|&v| f64::from(v)).collect();
                out.binary(
                    "supDataArrayBinary",
                    Precision::Bits32,
                    &values,
                    Some((array.name.as_str(), position + 1)),
                )?;
            }
        }
        out.raw("\t\t</spectrum>\n")?;
    }
    out.raw("\t</spectrumList>\n</mzData>\n")
}

fn warn(report: &mut StoreReport, message: &str) {
    report.warning_count += 1;
    report.warnings.push(message.to_owned());
}

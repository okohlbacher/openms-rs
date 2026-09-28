// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PEFF, the PSI Extended FASTA Format: protein sequences with structured
//! annotations for modifications, variants, processed regions and
//! proteoforms (<https://github.com/HUPO-PSI/PEFF>).
//!
//! Ports `FORMAT/PEFFFile.h` and its records `PEFFModification`,
//! `PEFFVariantSimple`, `PEFFVariantComplex`, `PEFFProcessedRegion`,
//! `PEFFDisulfideBond`, `PEFFEntry`, `PEFFCustomKeyDef` and
//! `PEFFDatabaseMetadata`. A file is a header of `#` lines — one file
//! description block and one database description block per database, each
//! ended by `# //` — followed by FASTA-style entries whose description line
//! carries `\Key=Value` annotations:
//!
//! ```text
//! # PEFF 1.0
//! # //
//! # DbName=TestDB
//! # Prefix=test
//! # //
//! >test:P12345 \PName=Test Protein 1 \ModResUnimod=(10|UNIMOD:35|Oxidation)
//! MSEQVENCE
//! ```
//!
//! [`PEFFFile::load`] and [`PEFFFile::store`] read and write a whole file;
//! [`PEFFFile::read_start`] / [`PEFFFile::read_next`] and
//! [`PEFFFile::write_start`] / [`PEFFFile::write_next`] /
//! [`PEFFFile::write_end`] stream one entry at a time. [`PEFFReader`] reads
//! from any [`BufRead`].
//!
//! The source reports specification violations it tolerates through
//! `OPENMS_LOG_WARN`; the port returns those lines, verbatim and in order, as
//! warnings (see [`PEFFFile::warnings`]). See `docs/PEFF_SUPPORT.md` for the
//! member mapping, the native differences and the evidence.

use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::format::fasta::FASTAEntry;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

mod read;
mod sequences;
mod text;
mod write;

pub use read::PEFFReader;
pub use sequences::{
    MAX_COMBINATORIAL_ELEMENTS, PEFFPeptides, PeptideGenerationOptions, SIGNAL_PEPTIDE,
    VariantDigestOptions,
};

/// The annotation identifier of an annotation that has none: source
/// `std::numeric_limits<UInt>::max()`.
pub const ANNOTATION_ID_NOT_SET: u32 = u32::MAX;

/// The vocabulary of a [`PEFFModification`]'s accession
/// (source `PEFFModification::Type`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PEFFModificationType {
    /// PSI-MOD, accessions `MOD:…`; written under `\ModResPsi`.
    PsiMod,
    /// Unimod, accessions `UNIMOD:…`; written under `\ModResUnimod`.
    Unimod,
    /// Anything else; written under `\ModRes`.
    #[default]
    Generic,
}

impl PEFFModificationType {
    /// The type the source constructor derives from an accession prefix.
    pub fn of_accession(accession: &str) -> Self {
        if accession.starts_with("MOD:") {
            Self::PsiMod
        } else if accession.starts_with("UNIMOD:") {
            Self::Unimod
        } else {
            Self::Generic
        }
    }
}

/// A PEFF modification annotation (`\ModResPsi`, `\ModResUnimod`, `\ModRes`).
///
/// Accessions may be PSI-MOD, Unimod or custom. The position is one-based,
/// and 0 means an unknown position (written `?`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFModification {
    /// One-based position; 0 is unknown (`?`).
    pub position: usize,
    /// `MOD:00046`, `UNIMOD:35` or a custom accession.
    pub accession: String,
    /// Human-readable name.
    pub name: String,
    /// The optional tag, the last field of the annotation tuple.
    pub optional_tag: String,
    /// The annotation identifier (under `HasAnnotationIdentifiers=true`), or
    /// [`ANNOTATION_ID_NOT_SET`].
    pub annotation_id: u32,
    /// The vocabulary; source member `type`.
    pub modification_type: PEFFModificationType,
}

impl Default for PEFFModification {
    fn default() -> Self {
        Self {
            position: 0,
            accession: String::new(),
            name: String::new(),
            optional_tag: String::new(),
            annotation_id: ANNOTATION_ID_NOT_SET,
            modification_type: PEFFModificationType::Generic,
        }
    }
}

impl PEFFModification {
    /// A modification with no optional tag and no annotation identifier; the
    /// type follows the accession prefix (`MOD:`, `UNIMOD:`, else generic).
    ///
    /// Source `PEFFModification(pos, acc, n)` with the defaulted tag and id.
    pub fn new(position: usize, accession: impl Into<String>, name: impl Into<String>) -> Self {
        Self::with_details(position, accession, name, "", ANNOTATION_ID_NOT_SET)
    }

    /// Source `PEFFModification(pos, acc, n, tag, aid)`.
    pub fn with_details(
        position: usize,
        accession: impl Into<String>,
        name: impl Into<String>,
        optional_tag: impl Into<String>,
        annotation_id: u32,
    ) -> Self {
        let accession = accession.into();
        Self {
            position,
            modification_type: PEFFModificationType::of_accession(&accession),
            accession,
            name: name.into(),
            optional_tag: optional_tag.into(),
            annotation_id,
        }
    }
}

/// A single amino acid substitution (`\VariantSimple`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFVariantSimple {
    /// One-based position.
    pub position: usize,
    /// The substituted residue; `'\0'` when unset.
    ///
    /// The source stores the first *byte* of the field in a `char`; this is
    /// the first character, which is the same for every ASCII code.
    pub variant_aa: char,
    /// The optional tag, the last field of the annotation tuple.
    pub optional_tag: String,
    /// The annotation identifier, or [`ANNOTATION_ID_NOT_SET`].
    pub annotation_id: u32,
}

impl Default for PEFFVariantSimple {
    fn default() -> Self {
        Self {
            position: 0,
            variant_aa: '\0',
            optional_tag: String::new(),
            annotation_id: ANNOTATION_ID_NOT_SET,
        }
    }
}

impl PEFFVariantSimple {
    /// Source `PEFFVariantSimple(pos, aa, tag, aid)`.
    pub fn new(
        position: usize,
        variant_aa: char,
        optional_tag: impl Into<String>,
        annotation_id: u32,
    ) -> Self {
        Self {
            position,
            variant_aa,
            optional_tag: optional_tag.into(),
            annotation_id,
        }
    }
}

/// An insertion, deletion or multi-residue substitution (`\VariantComplex`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFVariantComplex {
    /// One-based first replaced position.
    pub start_position: usize,
    /// One-based last replaced position.
    pub end_position: usize,
    /// The replacement sequence; empty for a deletion.
    pub replacement: String,
    /// The optional tag, the last field of the annotation tuple.
    pub optional_tag: String,
    /// The annotation identifier, or [`ANNOTATION_ID_NOT_SET`].
    pub annotation_id: u32,
}

impl Default for PEFFVariantComplex {
    fn default() -> Self {
        Self {
            start_position: 0,
            end_position: 0,
            replacement: String::new(),
            optional_tag: String::new(),
            annotation_id: ANNOTATION_ID_NOT_SET,
        }
    }
}

impl PEFFVariantComplex {
    /// Source `PEFFVariantComplex(start, end, repl, tag, aid)`.
    pub fn new(
        start_position: usize,
        end_position: usize,
        replacement: impl Into<String>,
        optional_tag: impl Into<String>,
        annotation_id: u32,
    ) -> Self {
        Self {
            start_position,
            end_position,
            replacement: replacement.into(),
            optional_tag: optional_tag.into(),
            annotation_id,
        }
    }
}

/// A processed region such as a signal or transit peptide (`\Processed`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFProcessedRegion {
    /// One-based start position.
    pub start_position: usize,
    /// One-based end position.
    pub end_position: usize,
    /// PEFF CV accession, for example `PEFF:0001021` (signal peptide).
    pub accession: String,
    /// Optional name, for example `signal peptide`.
    pub name: String,
    /// The optional tag, the last field of the annotation tuple.
    pub optional_tag: String,
    /// The annotation identifier, or [`ANNOTATION_ID_NOT_SET`].
    pub annotation_id: u32,
}

impl Default for PEFFProcessedRegion {
    fn default() -> Self {
        Self {
            start_position: 0,
            end_position: 0,
            accession: String::new(),
            name: String::new(),
            optional_tag: String::new(),
            annotation_id: ANNOTATION_ID_NOT_SET,
        }
    }
}

impl PEFFProcessedRegion {
    /// Source `PEFFProcessedRegion(start, end, acc, n, tag, aid)`.
    pub fn new(
        start_position: usize,
        end_position: usize,
        accession: impl Into<String>,
        name: impl Into<String>,
        optional_tag: impl Into<String>,
        annotation_id: u32,
    ) -> Self {
        Self {
            start_position,
            end_position,
            accession: accession.into(),
            name: name.into(),
            optional_tag: optional_tag.into(),
            annotation_id,
        }
    }
}

/// A disulfide bond between two cysteines (`\DisulfideBond`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFDisulfideBond {
    /// The first cysteine: its annotation identifier (or position).
    pub id1: String,
    /// The second cysteine: its annotation identifier (or position).
    pub id2: String,
    /// The optional tag, for example `between chains`.
    pub optional_tag: String,
    /// The annotation identifier, or [`ANNOTATION_ID_NOT_SET`].
    pub annotation_id: u32,
}

impl Default for PEFFDisulfideBond {
    fn default() -> Self {
        Self {
            id1: String::new(),
            id2: String::new(),
            optional_tag: String::new(),
            annotation_id: ANNOTATION_ID_NOT_SET,
        }
    }
}

impl PEFFDisulfideBond {
    /// Source `PEFFDisulfideBond(i1, i2, tag, aid)`.
    pub fn new(
        id1: impl Into<String>,
        id2: impl Into<String>,
        optional_tag: impl Into<String>,
        annotation_id: u32,
    ) -> Self {
        Self {
            id1: id1.into(),
            id2: id2.into(),
            optional_tag: optional_tag.into(),
            annotation_id,
        }
    }
}

/// One PEFF entry: a description line and its sequence, with every
/// annotation parsed.
///
/// Per the PEFF specification the description line is
/// `>Prefix:DbUniqueId \key=value \key=value …`, where `Prefix` names a
/// database of the header. `identifier` holds the whole `Prefix:DbUniqueId`,
/// and `prefix` only its first part.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PEFFEntry {
    /// The database prefix: `sp` for `>sp:P12345`.
    pub prefix: String,
    /// The whole identifier, `Prefix:DbUniqueId`.
    pub identifier: String,
    /// The sequence, white space removed.
    pub sequence: String,
    /// `\PName`, which may list several names.
    pub protein_names: Vec<String>,
    /// `\GName`.
    pub gene_name: String,
    /// `\NcbiTaxId` or its shorthand `\OX`.
    pub ncbi_tax_id: i32,
    /// `\TaxName`.
    pub taxonomy_name: String,
    /// `\Length`; the sequence length when the entry has none.
    pub sequence_length: usize,
    /// `\SV`.
    pub sequence_version: String,
    /// `\EV`.
    pub entry_version: String,
    /// `\PE`, 1 to 5.
    pub protein_existence: i32,
    /// `\DbUniqueId`.
    pub db_unique_id: String,
    /// `\ID`, for example `NPM_HUMAN`.
    pub entry_id: String,
    /// `\AltAC`: alternative accessions.
    pub alt_accessions: Vec<String>,
    /// `\ModResPsi`, `\ModResUnimod` and `\ModRes`.
    pub modifications: Vec<PEFFModification>,
    /// `\VariantSimple`.
    pub simple_variants: Vec<PEFFVariantSimple>,
    /// `\VariantComplex`.
    pub complex_variants: Vec<PEFFVariantComplex>,
    /// `\Processed`.
    pub processed_regions: Vec<PEFFProcessedRegion>,
    /// `\DisulfideBond`.
    pub disulfide_bonds: Vec<PEFFDisulfideBond>,
    /// `\Proteoform`: ProForma strings.
    pub proteoforms: Vec<String>,
    /// Every other key, with its value exactly as the line holds it
    /// (PEFF escapes and list parentheses included). Source `std::map`.
    pub custom_annotations: BTreeMap<String, String>,
}

impl PEFFEntry {
    /// The entry as FASTA: identifier, first protein name as the description,
    /// and sequence. Every PEFF annotation is dropped.
    ///
    /// Source `toFASTAEntry`.
    pub fn to_fasta_entry(&self) -> FASTAEntry {
        FASTAEntry::new(
            self.identifier.clone(),
            self.protein_names.first().cloned().unwrap_or_default(),
            self.sequence.clone(),
        )
    }

    /// A PEFF entry with the FASTA identifier, sequence and sequence length,
    /// and the description, when it is not empty, as the only protein name.
    ///
    /// Source `fromFASTAEntry`.
    pub fn from_fasta_entry(fasta: &FASTAEntry) -> Self {
        let mut entry = Self {
            identifier: fasta.identifier.clone(),
            sequence: fasta.sequence.clone(),
            sequence_length: fasta.sequence.len(),
            ..Self::default()
        };
        if !fasta.description.is_empty() {
            entry.protein_names.push(fasta.description.clone());
        }
        entry
    }
}

/// A `CustomKeyDef` of a database header: the declaration of a custom
/// annotation key.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PEFFCustomKeyDef {
    /// `KeyName`.
    pub key_name: String,
    /// `Description`, unquoted.
    pub description: String,
    /// `ConceptCURIE`.
    pub concept_curie: String,
    /// `RegExp`, unquoted.
    pub regexp: String,
    /// `FieldNames`.
    pub field_names: Vec<String>,
    /// `FieldTypes`.
    pub field_types: Vec<String>,
}

/// The residue alphabet of a database (source
/// `PEFFDatabaseMetadata::SequenceType`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PEFFSequenceType {
    /// Amino acids.
    #[default]
    Aa,
    /// Nucleic acids.
    Na,
}

/// The header of one database of a PEFF file: its `#` description block.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PEFFDatabaseMetadata {
    /// The PEFF version of the file, `1.0` by default.
    pub version: String,
    /// `DbName`.
    pub db_name: String,
    /// `Prefix`: the identifier prefix of this database's entries.
    pub prefix: String,
    /// `DbDescription`.
    pub db_description: String,
    /// `Decoy`.
    pub is_decoy: bool,
    /// `DbSource`, which may repeat.
    pub db_sources: Vec<String>,
    /// `DbVersion`.
    pub db_version: String,
    /// `NumberOfEntries`.
    pub number_of_entries: usize,
    /// `SequenceType`.
    pub sequence_type: PEFFSequenceType,
    /// `GeneralComment`, which may repeat.
    pub general_comments: Vec<String>,
    /// `Conversion`: conversion notes.
    pub conversion: String,
    /// `SpecificKey`: key name to description.
    pub specific_keys: BTreeMap<String, String>,
    /// `SpecificValue`: key name to type definition.
    pub specific_values: BTreeMap<String, String>,
    /// `OptionalTagDef`.
    pub optional_tag_defs: Vec<String>,
    /// `CustomKeyDef`.
    pub custom_key_defs: Vec<PEFFCustomKeyDef>,
    /// `HasAnnotationIdentifiers`: annotations carry `id:` prefixes.
    pub has_annotation_identifiers: bool,
    /// `ProteoformDb`: a proteoform database.
    pub is_proteoform_db: bool,
    /// Every other header key, kept so a round trip preserves it (for
    /// example UniProt's `DbDate`).
    pub unrecognized_keys: BTreeMap<String, String>,
}

impl Default for PEFFDatabaseMetadata {
    fn default() -> Self {
        Self {
            version: "1.0".to_owned(),
            db_name: String::new(),
            prefix: String::new(),
            db_description: String::new(),
            is_decoy: false,
            db_sources: Vec::new(),
            db_version: String::new(),
            number_of_entries: 0,
            sequence_type: PEFFSequenceType::Aa,
            general_comments: Vec::new(),
            conversion: String::new(),
            specific_keys: BTreeMap::new(),
            specific_values: BTreeMap::new(),
            optional_tag_defs: Vec::new(),
            custom_key_defs: Vec::new(),
            has_annotation_identifiers: false,
            is_proteoform_db: false,
            unrecognized_keys: BTreeMap::new(),
        }
    }
}

/// Resource bounds for reading PEFF. Exceeding one is [`Error::InvalidValue`];
/// entries read before the one that failed have been returned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Bytes of one header line or one description line.
    pub max_line_bytes: usize,
    /// Bytes of one entry's sequence.
    pub max_sequence_bytes: usize,
    /// Header lines.
    pub max_header_lines: usize,
    /// Annotation items of one entry (modifications, variants, regions,
    /// bonds, proteoforms, names and accessions together).
    pub max_items_per_entry: usize,
    /// Entries [`PEFFFile::load`] collects.
    pub max_entries: usize,
    /// Warnings kept; later ones are counted in
    /// [`PEFFReader::suppressed_warnings`].
    pub max_warnings: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_line_bytes: 64 * 1024 * 1024,
            max_sequence_bytes: 256 * 1024 * 1024,
            max_header_lines: 1_000_000,
            max_items_per_entry: 10_000_000,
            max_entries: 1 << 32,
            max_warnings: 100_000,
        }
    }
}

/// The result of [`PEFFFile::load`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PEFFDocument {
    /// The entries, in file order.
    pub entries: Vec<PEFFEntry>,
    /// One header per database block, in file order.
    pub headers: Vec<PEFFDatabaseMetadata>,
    /// The source's `OPENMS_LOG_WARN` lines, without line ends, in order.
    pub warnings: Vec<String>,
}

/// Reads and writes PEFF files (source `PEFFFile`, a `ProgressLogger`).
///
/// Whole files go through [`load`](Self::load) and [`store`](Self::store);
/// single entries through the `read_*` and `write_*` streaming calls, whose
/// input and output sessions are independent of each other.
#[derive(Default)]
pub struct PEFFFile {
    /// Resource bounds for every read.
    pub limits: ReadLimits,
    /// Progress reporting, as the inherited source `ProgressLogger`.
    pub progress: ProgressLogger,
    input: Option<PEFFReader<Box<dyn BufRead>>>,
    output: Option<BufWriter<File>>,
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// Source `File::exists` then `File::readable` then open, for read paths.
fn open(path: &Path) -> Result<BufReader<File>> {
    Ok(BufReader::new(File::open(path)?))
}

/// Source `FileHandler::hasValidExtension(filename, FileTypes::PEFF)`: a
/// known extension other than `.peff` is refused; an unknown one is not.
fn check_extension(path: &Path) -> Result<()> {
    use super::file_types::{FileType, has_valid_extension};
    if has_valid_extension(&path.to_string_lossy(), FileType::Peff) {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!(
            "{}: invalid file extension; expected 'peff'",
            path.display()
        )))
    }
}

impl PEFFFile {
    /// A PEFF file handler with the default [`ReadLimits`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Load every entry and header of a PEFF file.
    ///
    /// Source `load(filename, entries, headers)`, which reports the progress
    /// section "Loading PEFF file" from 0 to 1 and clears both vectors
    /// first; this returns them, with the warnings, only on success.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened (the source's
    /// `FileNotFound` / `FileNotReadable`); [`Error::Parse`] when an entry
    /// cannot be read or a number does not convert (`ParseError` /
    /// `ConversionError`); [`Error::InvalidValue`] when a [`ReadLimits`]
    /// bound is exceeded.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<PEFFDocument> {
        let limits = self.limits;
        let mut progress = ProgressReporter::new(Some(&mut self.progress));
        progress.start(0, 1, "Loading PEFF file")?;
        let document = read_with_limits(open(path.as_ref())?, &limits)?;
        progress.end()?;
        Ok(document)
    }

    /// Store `entries` under one database header.
    ///
    /// Source `store(filename, entries, header)`: a file description block
    /// (`# PEFF`, the header's general comments, `# //`), the database block,
    /// then every entry, with progress "Writing PEFF file" per entry.
    ///
    /// Returns the source's `OPENMS_LOG_WARN` lines for header keys the
    /// specification requires but the header lacks.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the name has a known extension other than
    /// `.peff`, and [`Error::Io`] when the file cannot be created or written
    /// — both the source's `UnableToCreateFile`. An I/O failure part-way
    /// leaves a partial file, as in the source.
    pub fn store(
        &mut self,
        path: impl AsRef<Path>,
        entries: &[PEFFEntry],
        header: &PEFFDatabaseMetadata,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        let text = write::format_header(header, &mut warnings);
        self.store_text(path.as_ref(), text, entries)?;
        Ok(warnings)
    }

    /// Store `entries` under several database headers: one file description
    /// block, one database block per header, then every entry.
    ///
    /// Source `store(filename, entries, headers)`. The PEFF version is the
    /// first non-empty header version (with a warning when the headers
    /// disagree), and the general comments of all headers are merged without
    /// duplicates, in first-seen order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `headers` is empty (the source's
    /// `InvalidParameter`); otherwise as [`store`](Self::store).
    pub fn store_multi(
        &mut self,
        path: impl AsRef<Path>,
        entries: &[PEFFEntry],
        headers: &[PEFFDatabaseMetadata],
    ) -> Result<Vec<String>> {
        if headers.is_empty() {
            return Err(invalid(
                "PEFF: store() requires at least one database header",
            ));
        }
        let mut warnings = Vec::new();
        let text = write::format_headers(headers, &mut warnings)?;
        self.store_text(path.as_ref(), text, entries)?;
        Ok(warnings)
    }

    fn store_text(&mut self, path: &Path, header: String, entries: &[PEFFEntry]) -> Result<()> {
        check_extension(path)?;
        let mut progress = ProgressReporter::new(Some(&mut self.progress));
        progress.start(0, progress_value(entries.len())?, "Writing PEFF file")?;
        let mut output = BufWriter::new(File::create(path)?);
        output.write_all(header.as_bytes())?;
        for entry in entries {
            output.write_all(write::format_entry(entry).as_bytes())?;
            progress.next_progress()?;
        }
        output.flush()?;
        progress.end()
    }

    /// Open a PEFF file for [`read_next`](Self::read_next) and parse its
    /// header, which [`headers`](Self::headers) then returns.
    ///
    /// Source `readStart`. After the header, every database block is checked
    /// for the keys the specification requires, and each missing one is a
    /// warning; a file without any header gets one default header.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened (the source's
    /// `FileNotFound` / `FileNotReadable`); [`Error::Parse`] when a header
    /// number does not convert; [`Error::InvalidValue`] when a bound is
    /// exceeded. On error the previous session is closed.
    pub fn read_start(&mut self, path: impl AsRef<Path>) -> Result<()> {
        self.input = None;
        let reader: Box<dyn BufRead> = Box::new(open(path.as_ref())?);
        self.input = Some(PEFFReader::with_limits(reader, self.limits)?);
        Ok(())
    }

    fn reader(&mut self) -> Result<&mut PEFFReader<Box<dyn BufRead>>> {
        self.input
            .as_mut()
            .ok_or_else(|| invalid("PEFF input is not open"))
    }

    /// Read the next entry into `entry`; `false` at the end of the file.
    ///
    /// Source `readNext`. On `true`, `entry` is replaced; on `false` or an
    /// error it is unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the source's messages when an entry cannot be
    /// read ("The first PEFF entry could not be read!", or "Only N entries
    /// could be read. Parsing next record failed.") or a number does not
    /// convert; [`Error::InvalidValue`] when no file is open or a bound is
    /// exceeded.
    pub fn read_next(&mut self, entry: &mut PEFFEntry) -> Result<bool> {
        match self.reader()?.next_entry()? {
            Some(next) => {
                *entry = next;
                Ok(true)
            }
            None => Ok(false),
        }
    }

    /// The headers parsed by [`read_start`](Self::read_start), or none when
    /// no file is open. Source `getHeaders`.
    pub fn headers(&self) -> &[PEFFDatabaseMetadata] {
        self.input.as_ref().map_or(&[], PEFFReader::headers)
    }

    /// The warnings of the open read session so far (header and entries).
    pub fn warnings(&self) -> &[String] {
        self.input.as_ref().map_or(&[], PEFFReader::warnings)
    }

    /// Whether the open input has no bytes left. Source `atEnd`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no file is open; [`Error::Io`] when the
    /// input fails.
    pub fn at_end(&mut self) -> Result<bool> {
        self.reader()?.at_end()
    }

    /// Create `path` and write the header of one database.
    ///
    /// Source `writeStart(filename, header)`. Returns the header's warnings.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when an output is already open or the name has
    /// a known extension other than `.peff`; [`Error::Io`] when the file
    /// cannot be created — both `UnableToCreateFile` in the source.
    pub fn write_start(
        &mut self,
        path: impl AsRef<Path>,
        header: &PEFFDatabaseMetadata,
    ) -> Result<Vec<String>> {
        let mut warnings = Vec::new();
        let text = write::format_header(header, &mut warnings);
        self.open_output(path.as_ref(), &text)?;
        Ok(warnings)
    }

    /// Create `path` and write the headers of several databases.
    ///
    /// Source `writeStart(filename, headers)`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `headers` is empty (the source's
    /// `InvalidParameter`); otherwise as [`write_start`](Self::write_start).
    pub fn write_start_multi(
        &mut self,
        path: impl AsRef<Path>,
        headers: &[PEFFDatabaseMetadata],
    ) -> Result<Vec<String>> {
        if headers.is_empty() {
            return Err(invalid(
                "PEFF: writeStart() requires at least one database header",
            ));
        }
        let mut warnings = Vec::new();
        let text = write::format_headers(headers, &mut warnings)?;
        self.open_output(path.as_ref(), &text)?;
        Ok(warnings)
    }

    fn open_output(&mut self, path: &Path, header: &str) -> Result<()> {
        if self.output.is_some() {
            return Err(invalid("PEFF output is already open"));
        }
        check_extension(path)?;
        let mut output = BufWriter::new(File::create(path)?);
        output.write_all(header.as_bytes())?;
        self.output = Some(output);
        Ok(())
    }

    /// Append one entry to the open output. Source `writeNext`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no output is open; [`Error::Io`] when the
    /// write fails.
    pub fn write_next(&mut self, entry: &PEFFEntry) -> Result<()> {
        let output = self
            .output
            .as_mut()
            .ok_or_else(|| invalid("PEFF output is not open"))?;
        output.write_all(write::format_entry(entry).as_bytes())?;
        Ok(())
    }

    /// Flush and close the open output; nothing happens when none is open.
    /// Source `writeEnd`, which the source destructor also calls.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the final flush fails.
    pub fn write_end(&mut self) -> Result<()> {
        if let Some(mut output) = self.output.take() {
            output.flush()?;
        }
        Ok(())
    }

    /// Whether `path` looks like PEFF: the first non-blank line is a `#` line
    /// holding `PEFF` and a version digit, or an entry line with a backslash
    /// annotation. Unreadable files are not PEFF. Source `isPEFFFile`.
    pub fn is_peff_file(path: impl AsRef<Path>) -> bool {
        let Ok(file) = File::open(path) else {
            return false;
        };
        read::is_peff(BufReader::new(file))
    }

    /// The entry's sequence in ProForma notation, with its modifications:
    /// each localised one right after its residue, unlocalised ones
    /// (position 0) as a `<?[…]>` prefix, each as `[accession]`, or `[name]`
    /// without an accession. Source `toProForma`.
    pub fn to_proforma(entry: &PEFFEntry) -> String {
        write::to_proforma(entry)
    }
}

impl Drop for PEFFFile {
    fn drop(&mut self) {
        // The source destructor closes the output stream; a flush failure
        // here has nowhere to go, as in the source.
        let _ = self.write_end();
    }
}

/// Read a whole PEFF document from any buffered input.
///
/// # Errors
///
/// As [`PEFFFile::load`].
pub fn read(input: impl BufRead) -> Result<PEFFDocument> {
    read_with_limits(input, &ReadLimits::default())
}

/// [`read`] under explicit resource bounds.
///
/// # Errors
///
/// As [`PEFFFile::load`].
pub fn read_with_limits(input: impl BufRead, limits: &ReadLimits) -> Result<PEFFDocument> {
    let mut reader = PEFFReader::with_limits(input, *limits)?;
    let mut entries = Vec::new();
    while let Some(entry) = reader.next_entry()? {
        if entries.len() == limits.max_entries {
            return Err(invalid("PEFF entry limit exceeded"));
        }
        entries.push(entry);
    }
    let (headers, warnings) = reader.into_parts();
    Ok(PEFFDocument {
        entries,
        headers,
        warnings,
    })
}

/// Write a whole PEFF document (one database header) to any output.
///
/// Returns the header warnings. Source `store` without the file.
///
/// # Errors
///
/// [`Error::Io`] when the output fails.
pub fn write(
    mut output: impl Write,
    entries: &[PEFFEntry],
    header: &PEFFDatabaseMetadata,
) -> Result<Vec<String>> {
    let mut warnings = Vec::new();
    output.write_all(write::format_header(header, &mut warnings).as_bytes())?;
    for entry in entries {
        output.write_all(write::format_entry(entry).as_bytes())?;
    }
    output.flush()?;
    Ok(warnings)
}

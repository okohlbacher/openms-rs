// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! UniProtKB XML protein databases, read into the fields a PEFF header needs.
//!
//! Ports `FORMAT/UniProtXMLFile.h` (with its two records `UniProtEntry` and
//! `UniProtFeature`) and the handler `FORMAT/HANDLERS/UniProtXMLHandler.h`.
//! The source handler follows the scoping rules of the upstream C# UniPEFF
//! tool, and so does this one: isoform sequences and the
//! `<comment type="alternative products">` subtree are skipped, `<gene>` and
//! `<organism>` names are scoped so they never leak into the entry mnemonic,
//! and a `<feature>` is collected whole before it is appended.
//!
//! Input may be plain or gzip/bzip2 compressed; compression is recognised by
//! the leading bytes, as the source's `XMLFile::parse_` does. See
//! `docs/UNIPROT_XML_SUPPORT.md` for the member mapping and evidence.

use super::sax_xml::{SaxAttributes, SaxEvent, SaxLimits, parse};
use crate::{Error, Result};
use std::io::BufRead;
use std::path::Path;

/// One `<feature>` of a UniProtKB entry, after its `<location>` has been
/// collapsed.
///
/// A feature is either a single point (`has_position`, `position`) or a range
/// (`has_range`, `begin`, `end`). A coordinate whose `position` attribute is
/// absent (UniProt's `status="unknown"`) or not a non-negative decimal is 0;
/// the consumer decides whether that is legal for its PEFF key.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UniProtFeature {
    /// `feature/@type`, for example `modified residue`, `disulfide bond` or
    /// `sequence variant`. Source member `type`.
    pub feature_type: String,
    /// `feature/@description`, raw, before any cleanup.
    pub description: String,
    /// A `<position>` element was present.
    pub has_position: bool,
    /// A `<begin>` or `<end>` element was present.
    pub has_range: bool,
    /// One-based position; 0 when unknown or absent.
    pub position: i32,
    /// One-based range start; 0 when unknown or absent.
    pub begin: i32,
    /// One-based range end; 0 when unknown or absent.
    pub end: i32,
    /// `<original>` text of a sequence variant (the last one, when repeated).
    pub original: String,
    /// `<variation>` text of a sequence variant; only the first is kept.
    pub variation: String,
}

/// One `<entry>` of a UniProtKB XML file, in a parser-neutral form.
///
/// Holds only the fields needed to build a PEFF description line. Isoforms
/// (sequences without a `length` attribute) and splice-variant features are
/// deliberately not represented, because UniPEFF's PEFF output ignores them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UniProtEntry {
    /// `entry/@dataset`, for example `Swiss-Prot` or `TrEMBL`.
    pub dataset: String,
    /// The first `<accession>`, the primary identifier.
    pub accession: String,
    /// The second and later `<accession>` elements.
    pub alt_accessions: Vec<String>,
    /// The first `<name>` directly under `<entry>`: the mnemonic, for example
    /// `KSINK_HUMAN`.
    pub name: String,
    /// The first `<fullName>` inside a `<recommendedName>`.
    pub full_name: String,
    /// The first `<gene>/<name type="primary">`.
    pub primary_gene: String,
    /// `<organism>/<dbReference type="NCBI Taxonomy" id="...">`.
    pub ncbi_tax_id: String,
    /// `<organism>/<name type="scientific">`.
    pub tax_name: String,
    /// `<proteinExistence type="...">`, raw, for example
    /// `evidence at protein level`.
    pub protein_existence: String,
    /// `<sequence version="...">` of the canonical sequence.
    pub sequence_version: String,
    /// `<entry version="...">`.
    pub entry_version: String,
    /// The canonical `<sequence>` text with all white space removed; isoform
    /// sequences are skipped.
    pub sequence: String,
    /// `<feature>` elements in document order.
    pub features: Vec<UniProtFeature>,
}

impl UniProtEntry {
    /// Bytes of text this entry holds, for the aggregate bound of [`read`].
    fn text_bytes(&self) -> usize {
        let strings = [
            &self.dataset,
            &self.accession,
            &self.name,
            &self.full_name,
            &self.primary_gene,
            &self.ncbi_tax_id,
            &self.tax_name,
            &self.protein_existence,
            &self.sequence_version,
            &self.entry_version,
            &self.sequence,
        ];
        let features: usize = self
            .features
            .iter()
            .map(|f| {
                f.feature_type.len()
                    + f.description.len()
                    + f.original.len()
                    + f.variation.len()
                    + std::mem::size_of::<UniProtFeature>()
            })
            .sum();
        strings.iter().map(|s| s.len()).sum::<usize>()
            + self.alt_accessions.iter().map(String::len).sum::<usize>()
            + features
    }
}

/// Resource bounds for one read. Exceeding any of them is
/// [`Error::InvalidValue`]; a streaming read has then delivered the entries
/// before the one that failed, and [`read`] returns none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Decompressed input bytes. The default admits a full TrEMBL release,
    /// which a streaming read handles one entry at a time.
    pub max_input_bytes: u64,
    /// Open elements at any one time.
    pub max_depth: usize,
    /// Parser events in the document.
    pub max_events: u64,
    /// Bytes of one run of character data (the longest is a sequence).
    pub max_text_bytes: usize,
    /// Features of one entry.
    pub max_features_per_entry: usize,
    /// Accessions of one entry.
    pub max_accessions_per_entry: usize,
    /// Entries [`read`] collects.
    pub max_entries: usize,
    /// Bytes of text [`read`] collects over all entries.
    pub max_result_bytes: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 1 << 40,
            max_depth: 256,
            max_events: 1 << 40,
            max_text_bytes: 64 * 1024 * 1024,
            max_features_per_entry: 1_000_000,
            max_accessions_per_entry: 1_000_000,
            max_entries: 100_000_000,
            max_result_bytes: 16 * 1024 * 1024 * 1024,
        }
    }
}

impl ReadLimits {
    fn sax(&self) -> SaxLimits {
        SaxLimits {
            max_input_bytes: self.max_input_bytes,
            max_depth: self.max_depth,
            max_events: self.max_events,
            max_attributes: 1024,
            max_text_bytes: self.max_text_bytes,
        }
    }
}

/// Where the buffered character data goes at the next end tag
/// (source `UniProtXMLHandler::CaptureTarget`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Capture {
    None,
    Accession,
    EntryName,
    FullName,
    PrimaryGene,
    TaxName,
    Sequence,
    FeatureOriginal,
    FeatureVariation,
}

/// Source `parsePosition_`: `strtol` in base 10 over the whole attribute,
/// with a negative, partial or empty parse meaning unknown (0).
///
/// `strtol` skips leading white space and accepts one sign. A value beyond
/// `int` is 0 here; the source casts it with `static_cast<int>`, which keeps
/// the low 32 bits (and turns an overflowing `strtol`'s `LONG_MAX` into -1).
fn parse_position(attribute: &str) -> i32 {
    let text = attribute.trim_start_matches([' ', '\t', '\n', '\x0B', '\x0C', '\r']);
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return 0;
    }
    if negative {
        return 0;
    }
    digits.parse::<i32>().unwrap_or(0)
}

/// Source `stripWhitespace`: every byte `std::isspace` accepts in the C
/// locale, which is ASCII space, tab, line feed, vertical tab, form feed and
/// carriage return.
fn strip_whitespace(text: &str) -> String {
    text.chars()
        .filter(|c| !matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r'))
        .collect()
}

/// The state of source `UniProtXMLHandler`, member for member.
struct Handler<'a, F> {
    callback: F,
    limits: &'a ReadLimits,
    entry: UniProtEntry,
    depth: usize,
    entry_depth: usize,
    recommended_name_depth: usize,
    gene_depth: usize,
    organism_depth: usize,
    alt_products_depth: usize,
    feature_depth: usize,
    sequence_depth: usize,
    full_name_captured: bool,
    buffer: String,
    capture: Capture,
    feature: UniProtFeature,
}

fn attribute(attributes: &SaxAttributes, name: &str) -> String {
    attributes.get(name).unwrap_or("").to_owned()
}

impl<F: FnMut(UniProtEntry) -> Result<()>> Handler<'_, F> {
    /// `resetEntry_`.
    fn reset_entry(&mut self) {
        self.entry = UniProtEntry::default();
        self.full_name_captured = false;
        self.recommended_name_depth = 0;
        self.gene_depth = 0;
        self.organism_depth = 0;
        self.sequence_depth = 0;
        self.feature_depth = 0;
    }

    fn event(&mut self, event: SaxEvent<'_>) -> Result<()> {
        match event {
            SaxEvent::Start(tag, attributes) => self.start(tag, attributes),
            SaxEvent::End(tag) => self.end(tag),
            SaxEvent::Characters(chars) => {
                if self.capture != Capture::None {
                    self.buffer.push_str(chars);
                }
                Ok(())
            }
        }
    }

    fn start(&mut self, tag: &str, attributes: &SaxAttributes) -> Result<()> {
        self.depth += 1;
        if self.alt_products_depth != 0 {
            return Ok(());
        }
        self.buffer.clear();
        self.capture = Capture::None;

        if tag == "entry" && self.entry_depth == 0 {
            self.reset_entry();
            self.entry_depth = self.depth;
            self.entry.dataset = attribute(attributes, "dataset");
            self.entry.entry_version = attribute(attributes, "version");
            return Ok(());
        }
        if self.entry_depth == 0 {
            return Ok(());
        }
        if tag == "comment" && attributes.get("type") == Some("alternative products") {
            self.alt_products_depth = self.depth;
            return Ok(());
        }
        if tag == "feature" && self.feature_depth == 0 {
            if self.entry.features.len() == self.limits.max_features_per_entry {
                return Err(Error::InvalidValue(
                    "UniProt XML feature limit exceeded".into(),
                ));
            }
            self.feature = UniProtFeature::default();
            self.feature_depth = self.depth;
            self.feature.feature_type = attribute(attributes, "type");
            self.feature.description = attribute(attributes, "description");
            return Ok(());
        }
        if self.feature_depth != 0 {
            match tag {
                "position" => {
                    self.feature.has_position = true;
                    self.feature.position = parse_position(attributes.get("position").unwrap_or(""));
                }
                "begin" => {
                    self.feature.has_range = true;
                    self.feature.begin = parse_position(attributes.get("position").unwrap_or(""));
                }
                "end" => {
                    self.feature.has_range = true;
                    self.feature.end = parse_position(attributes.get("position").unwrap_or(""));
                }
                "original" => self.capture = Capture::FeatureOriginal,
                // The C# tool keeps only the first <variation>.
                "variation" if self.feature.variation.is_empty() => {
                    self.capture = Capture::FeatureVariation;
                }
                _ => {}
            }
            return Ok(());
        }
        if tag == "accession" {
            self.capture = Capture::Accession;
            return Ok(());
        }
        if tag == "recommendedName" && self.recommended_name_depth == 0 {
            self.recommended_name_depth = self.depth;
            return Ok(());
        }
        if tag == "fullName" && self.recommended_name_depth != 0 && !self.full_name_captured {
            self.capture = Capture::FullName;
            return Ok(());
        }
        if tag == "gene" && self.gene_depth == 0 {
            self.gene_depth = self.depth;
            return Ok(());
        }
        if self.gene_depth != 0
            && tag == "name"
            && attributes.get("type") == Some("primary")
            && self.entry.primary_gene.is_empty()
        {
            self.capture = Capture::PrimaryGene;
            return Ok(());
        }
        if tag == "organism" && self.organism_depth == 0 {
            self.organism_depth = self.depth;
            return Ok(());
        }
        if self.organism_depth != 0
            && tag == "name"
            && attributes.get("type") == Some("scientific")
            && self.entry.tax_name.is_empty()
        {
            self.capture = Capture::TaxName;
            return Ok(());
        }
        if self.organism_depth != 0
            && tag == "dbReference"
            && attributes.get("type") == Some("NCBI Taxonomy")
            && self.entry.ncbi_tax_id.is_empty()
        {
            self.entry.ncbi_tax_id = attribute(attributes, "id");
            return Ok(());
        }
        if tag == "proteinExistence" && self.entry.protein_existence.is_empty() {
            self.entry.protein_existence = attribute(attributes, "type");
            return Ok(());
        }
        // Only the canonical sequence carries a length attribute.
        if tag == "sequence"
            && self.depth == self.entry_depth + 1
            && !attributes.get("length").unwrap_or("").is_empty()
        {
            self.sequence_depth = self.depth;
            if self.entry.sequence_version.is_empty() {
                self.entry.sequence_version = attribute(attributes, "version");
            }
            self.capture = Capture::Sequence;
            return Ok(());
        }
        if tag == "name"
            && self.depth == self.entry_depth + 1
            && self.gene_depth == 0
            && self.organism_depth == 0
            && self.recommended_name_depth == 0
            && self.entry.name.is_empty()
        {
            self.capture = Capture::EntryName;
        }
        Ok(())
    }

    fn end(&mut self, tag: &str) -> Result<()> {
        if self.alt_products_depth != 0 {
            if self.depth == self.alt_products_depth && tag == "comment" {
                self.alt_products_depth = 0;
            }
            self.depth -= 1;
            return Ok(());
        }
        if self.capture != Capture::None {
            let text = std::mem::take(&mut self.buffer);
            match self.capture {
                Capture::Accession => {
                    if self.entry.accession.is_empty() {
                        self.entry.accession = text;
                    } else {
                        if self.entry.alt_accessions.len() + 1
                            >= self.limits.max_accessions_per_entry
                        {
                            return Err(Error::InvalidValue(
                                "UniProt XML accession limit exceeded".into(),
                            ));
                        }
                        self.entry.alt_accessions.push(text);
                    }
                }
                Capture::EntryName => self.entry.name = text,
                Capture::FullName => {
                    self.entry.full_name = text;
                    self.full_name_captured = true;
                }
                Capture::PrimaryGene => self.entry.primary_gene = text,
                Capture::TaxName => self.entry.tax_name = text,
                Capture::Sequence => self.entry.sequence = strip_whitespace(&text),
                Capture::FeatureOriginal => self.feature.original = text,
                Capture::FeatureVariation => self.feature.variation = text,
                Capture::None => {}
            }
            self.capture = Capture::None;
        }
        if self.sequence_depth != 0 && self.depth == self.sequence_depth && tag == "sequence" {
            self.sequence_depth = 0;
        }
        if self.recommended_name_depth != 0
            && self.depth == self.recommended_name_depth
            && tag == "recommendedName"
        {
            self.recommended_name_depth = 0;
        }
        if self.gene_depth != 0 && self.depth == self.gene_depth && tag == "gene" {
            self.gene_depth = 0;
        }
        if self.organism_depth != 0 && self.depth == self.organism_depth && tag == "organism" {
            self.organism_depth = 0;
        }
        if self.feature_depth != 0 && self.depth == self.feature_depth && tag == "feature" {
            let feature = std::mem::take(&mut self.feature);
            self.entry.features.push(feature);
            self.feature_depth = 0;
        }
        if self.entry_depth != 0 && self.depth == self.entry_depth && tag == "entry" {
            let entry = std::mem::take(&mut self.entry);
            (self.callback)(entry)?;
            self.reset_entry();
            self.entry_depth = 0;
        }
        self.depth -= 1;
        Ok(())
    }
}

/// Read every entry of a UniProtKB XML document into memory, with the
/// default [`ReadLimits`].
///
/// Source `UniProtXMLFile::load`, which clears its output vector and then
/// streams into it; this returns the entries, and none if the read fails.
///
/// # Errors
///
/// As [`read_streaming`]; also [`Error::InvalidValue`] when the entries or
/// their text exceed [`ReadLimits::max_entries`] or
/// [`ReadLimits::max_result_bytes`].
pub fn read(input: impl BufRead) -> Result<Vec<UniProtEntry>> {
    read_with_limits(input, &ReadLimits::default())
}

/// [`read`] under explicit resource bounds.
///
/// # Errors
///
/// As [`read`].
pub fn read_with_limits(input: impl BufRead, limits: &ReadLimits) -> Result<Vec<UniProtEntry>> {
    let mut entries = Vec::new();
    let mut bytes = 0usize;
    read_streaming_with_limits(input, limits, |entry| {
        bytes = bytes.saturating_add(entry.text_bytes());
        if entries.len() == limits.max_entries || bytes > limits.max_result_bytes {
            return Err(Error::InvalidValue(
                "UniProt XML result limit exceeded".into(),
            ));
        }
        entries.push(entry);
        Ok(())
    })?;
    Ok(entries)
}

/// Stream the entries of a UniProtKB XML document to `callback`, one at a
/// time as each `</entry>` is reached, with the default [`ReadLimits`].
///
/// Source `UniProtXMLFile::loadStreaming`, the path for databases too large
/// to hold in memory (TrEMBL). The callback owns each entry. It returns a
/// [`Result`] here, and an error stops the read and is returned unchanged;
/// the source's `std::function<void(UniProtEntry&&)>` cannot stop it.
///
/// The field rules of the source handler, all preserved:
///
/// - Only elements inside an `<entry>` are read; an `<entry>` nested in
///   another is ordinary content.
/// - The whole `<comment type="alternative products">` subtree is skipped.
/// - `<accession>`: the first is `accession`, later ones `alt_accessions`.
/// - `name`: the first `<name>` directly under `<entry>` that is not inside
///   `<gene>`, `<organism>` or `<recommendedName>`.
/// - `full_name`: the first `<fullName>` inside any `<recommendedName>`.
/// - `primary_gene` / `tax_name`: the first `<name type="primary">` inside a
///   `<gene>` and the first `<name type="scientific">` inside an
///   `<organism>`, wherever they sit in those elements.
/// - `ncbi_tax_id`: the `id` of the first `<dbReference type="NCBI
///   Taxonomy">` inside an `<organism>`.
/// - `sequence`: the text of a `<sequence>` directly under `<entry>` with a
///   non-empty `length` attribute, white space removed; its `version`, when
///   the entry has none yet, is `sequence_version`.
/// - Inside a `<feature>`, only `<position>`, `<begin>`, `<end>`,
///   `<original>` and the first `<variation>` are read.
/// - An element's text is kept only when no child element starts inside it,
///   because every start tag discards the buffered text; a missing attribute
///   reads as the empty string.
///
/// # Errors
///
/// [`Error::Parse`] when the document is not well-formed XML, the source's
/// `Exception::ParseError`; [`Error::Unsupported`] for a DTD or a non-UTF-8
/// encoding; [`Error::InvalidValue`] when a bound in [`ReadLimits`] is
/// exceeded; [`Error::Io`] when the input fails; or the callback's error.
pub fn read_streaming(
    input: impl BufRead,
    callback: impl FnMut(UniProtEntry) -> Result<()>,
) -> Result<()> {
    read_streaming_with_limits(input, &ReadLimits::default(), callback)
}

/// [`read_streaming`] under explicit resource bounds.
///
/// # Errors
///
/// As [`read_streaming`].
pub fn read_streaming_with_limits(
    input: impl BufRead,
    limits: &ReadLimits,
    callback: impl FnMut(UniProtEntry) -> Result<()>,
) -> Result<()> {
    let mut handler = Handler {
        callback,
        limits,
        entry: UniProtEntry::default(),
        depth: 0,
        entry_depth: 0,
        recommended_name_depth: 0,
        gene_depth: 0,
        organism_depth: 0,
        alt_products_depth: 0,
        feature_depth: 0,
        sequence_depth: 0,
        full_name_captured: false,
        buffer: String::new(),
        capture: Capture::None,
        feature: UniProtFeature::default(),
    };
    parse(input, &limits.sax(), |event| handler.event(event))
}

/// Load every entry of a UniProtKB XML file (`.xml`, `.xml.gz` or
/// `.xml.bz2`) with the default [`ReadLimits`].
///
/// # Errors
///
/// [`Error::Io`] when the file cannot be opened, the source's
/// `Exception::FileNotFound`; otherwise as [`read`].
pub fn load(path: impl AsRef<Path>) -> Result<Vec<UniProtEntry>> {
    load_with_limits(path, &ReadLimits::default())
}

/// [`load`] under explicit resource bounds.
///
/// # Errors
///
/// As [`load`].
pub fn load_with_limits(path: impl AsRef<Path>, limits: &ReadLimits) -> Result<Vec<UniProtEntry>> {
    read_with_limits(super::path_io::open(path.as_ref())?, limits)
}

/// Stream the entries of a UniProtKB XML file to `callback`, with the
/// default [`ReadLimits`]. Source `UniProtXMLFile::loadStreaming`.
///
/// # Errors
///
/// [`Error::Io`] when the file cannot be opened; otherwise as
/// [`read_streaming`].
pub fn load_streaming(
    path: impl AsRef<Path>,
    callback: impl FnMut(UniProtEntry) -> Result<()>,
) -> Result<()> {
    load_streaming_with_limits(path, &ReadLimits::default(), callback)
}

/// [`load_streaming`] under explicit resource bounds.
///
/// # Errors
///
/// As [`load_streaming`].
pub fn load_streaming_with_limits(
    path: impl AsRef<Path>,
    limits: &ReadLimits,
    callback: impl FnMut(UniProtEntry) -> Result<()>,
) -> Result<()> {
    read_streaming_with_limits(super::path_io::open(path.as_ref())?, limits, callback)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positions_follow_strtol_with_a_whole_string_rule() {
        for (text, expected) in [
            ("", 0),
            ("7", 7),
            ("  7", 7),
            ("+7", 7),
            ("-7", 0),
            ("7 ", 0),
            ("0x7", 0),
            ("2147483647", i32::MAX),
            ("2147483648", 0),
        ] {
            assert_eq!(parse_position(text), expected, "{text:?}");
        }
    }

    #[test]
    fn text_with_a_child_element_is_not_captured() {
        let text = "<uniprot><entry><accession>A<b/>B</accession><accession>C</accession>\
                    <name>N</name></entry></uniprot>";
        let entries = read(text.as_bytes()).unwrap();
        // The child <b> discarded "A"; its own end deposited nothing, so the
        // first accession is "C".
        assert_eq!(entries[0].accession, "C");
        assert!(entries[0].alt_accessions.is_empty());
        assert_eq!(entries[0].name, "N");
    }

    #[test]
    fn callback_errors_stop_the_read() {
        let text = "<uniprot><entry><accession>A</accession></entry>\
                    <entry><accession>B</accession></entry></uniprot>";
        let mut seen = Vec::new();
        let result = read_streaming(text.as_bytes(), |entry| {
            seen.push(entry.accession);
            Err(Error::InvalidValue("stop".into()))
        });
        assert!(matches!(result, Err(Error::InvalidValue(_))));
        assert_eq!(seen, ["A"]);
    }
}

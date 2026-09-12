// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! qcML quality-control reports: runs, sets, quality parameters and attachments.
//!
//! Port of `FORMAT/QcMLFile.h` and `FORMAT/QcMLFile.cpp`. A qcML document holds
//! `<runQuality>` and `<setQuality>` entries, each keyed by an `ID` and carrying
//! `<qualityParameter>` records (a CV-identified scalar) and `<attachment>`
//! records (either an opaque base64 `<binary>` payload or an inline `<table>`
//! with its own column types and space-delimited rows). A set additionally
//! names the runs it aggregates.
//!
//! Entry points:
//! [`load`](crate::format::qcml::load) and
//! [`read`](crate::format::qcml::read) parse a document,
//! [`QcMLFile::store`](crate::format::qcml::QcMLFile::store) writes one, and
//! [`QcMLFile`](crate::format::qcml::QcMLFile) is the in-memory document with
//! the source's registration, lookup, removal, merge and export surface.
//!
//! There is no bundled XSD: the source constructs its `XMLFile` base with an
//! empty schema path and version `"0.7"` because the qcML schema is archived,
//! so neither the source nor this port validates against a schema. See
//! `docs/QCML_SUPPORT.md` for the full member-by-member mapping, the source
//! defects this port corrects and the ones it deliberately reproduces.

use super::parse_error;
use crate::{Error, Result};
use quick_xml::{
    Reader,
    events::{BytesStart, Event},
};
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;
use std::path::Path;

/// Format version the source's `XMLFile` base is constructed with.
pub const VERSION: &str = "0.7";

/// Value the source's `exportQP` substitutes for a parameter it cannot find.
pub const NOT_FOUND: &str = "N/A";

/// CV accession whose value names a run (`MS:1000577`, "spectrum file name"),
/// and inside a set names one of its members.
pub const RUN_NAME_ACCESSION: &str = "MS:1000577";

/// CV accession whose value names a set (`QC:0000058`).
pub const SET_NAME_ACCESSION: &str = "QC:0000058";

/// CV accession the writer synthesises for each documented set member
/// (`QC:0000005`, written with the name `set name`).
pub const SET_MEMBER_ACCESSION: &str = "QC:0000005";

/// Attributes accepted on one element; no qcML element declares more than nine,
/// and the duplicate-name scan is quadratic in this count.
const MAX_ATTRIBUTES: usize = 64;

fn limit(message: &str) -> Error {
    Error::InvalidRange(message.into())
}

fn xml_char(c: char) -> bool {
    matches!(c, '\t' | '\n' | '\r')
        || ('\u{20}'..='\u{d7ff}').contains(&c)
        || ('\u{e000}'..='\u{fffd}').contains(&c)
        || c >= '\u{10000}'
}

/// Ceilings for one bounded qcML parse.
///
/// Every field is a hard refusal, checked before the corresponding allocation.
/// The source has no ceiling of any kind on this path: `parse_()` hands the
/// whole file to the XML parser and the handler appends to `std::vector` and
/// `std::map` members until the allocator fails.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    /// Maximum decoded document size in bytes.
    pub max_input_bytes: usize,
    /// Maximum element nesting depth, which an embedded XSL stylesheet uses up.
    pub max_depth: usize,
    /// Maximum number of markup nodes in the document.
    ///
    /// Elements, and also the comments and processing instructions the port
    /// discards: each costs the parser the same walk as an element, so a
    /// document of nothing but comments would otherwise be bounded only by
    /// [`max_input_bytes`](Self::max_input_bytes).
    pub max_elements: usize,
    /// Maximum bytes of character data captured for one `<binary>`,
    /// `<tableColumnTypes>` or `<tableRowValues>` element.
    pub max_text_bytes: usize,
    /// Maximum bytes of an internal DOCTYPE subset, which the source's own
    /// writer emits when it injects a report stylesheet.
    pub max_doctype_bytes: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_depth: 100,
            max_elements: 4_000_000,
            max_text_bytes: 64 * 1024 * 1024,
            max_doctype_bytes: 4096,
        }
    }
}

/// Representation of a quality parameter: one CV-identified scalar.
///
/// All eight fields are public `std::string` members in the source and carry no
/// invariant of their own; `store` and
/// [`to_xml_string`](QualityParameter::to_xml_string) impose the format's
/// requirements. The source's `flag` field is documented as "cv accession of
/// the unit", which is a copy-paste slip: the writer treats it as a boolean
/// marker and the reader parses the `flag` attribute into it verbatim.
///
/// # Comparison
///
/// The source's `operator==`, `operator<` and `operator>` compare `name` alone.
/// `PartialEq` here is full structural equality, because a Rust `==` that
/// ignored seven of eight fields would be a trap; use
/// [`same_name`](QualityParameter::same_name) for the source predicate.
/// `Ord` orders by `name` first and breaks ties on the remaining fields in
/// declaration order, so it agrees with the source's ordering on distinct names
/// and is deterministic where the source's `std::sort` is not.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct QualityParameter {
    /// Name, written as the `name` attribute. Required by the reader.
    pub name: String,
    /// Identifier, written as the `ID` attribute. Required by the reader.
    pub id: String,
    /// Value, written as the `value` attribute and omitted when empty.
    pub value: String,
    /// Controlled vocabulary reference, written as `cvRef`. Required.
    pub cv_ref: String,
    /// Controlled vocabulary accession, written as `accession`. Required.
    pub cv_acc: String,
    /// Unit's controlled vocabulary reference, omitted when empty.
    pub unit_ref: String,
    /// Unit's controlled vocabulary accession, omitted when empty.
    pub unit_acc: String,
    /// Boolean marker, omitted when empty. The source writes the literal
    /// `flag="true"` for any non-empty value and so cannot round-trip it;
    /// see [`WriteOptions::source_flag_literal`].
    pub flag: String,
}

impl QualityParameter {
    /// Maximum bytes accepted in any single field.
    pub const MAX_TEXT_BYTES: usize = 16 * 1024 * 1024;

    /// Maximum indentation level accepted by the XML writers.
    ///
    /// The source builds `std::string indent(indentation_level, '\t')` from an
    /// unchecked `UInt`, so a large level attempts a multi-gigabyte allocation.
    pub const MAX_INDENTATION: u32 = 64;

    /// The source's `operator==`, `operator<` and `operator>` predicate: `name`
    /// equality alone, ignoring identifier, value, vocabulary and unit.
    ///
    /// `merge` applies this through `std::unique`, which is why
    /// [`MergeOptions::source`] collapses distinct parameters that share a name.
    pub fn same_name(&self, other: &Self) -> bool {
        self.name == other.name
    }

    /// Serialise as a self-closing `<qualityParameter/>` element indented with
    /// `indentation_level` tabs, under the native default options.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `indentation_level` exceeds
    /// [`MAX_INDENTATION`](Self::MAX_INDENTATION) or a field exceeds
    /// [`MAX_TEXT_BYTES`](Self::MAX_TEXT_BYTES), [`Error::MissingInformation`]
    /// when `name`, `id`, `cv_ref` or `cv_acc` is empty, and
    /// [`Error::InvalidValue`] when a field holds a character XML 1.0 cannot
    /// represent. The source checks none of these and concatenates the values
    /// verbatim, so an empty required attribute or an unescaped `&` produces a
    /// document its own reader rejects.
    pub fn to_xml_string(&self, indentation_level: u32) -> Result<String> {
        self.to_xml_string_with_options(indentation_level, &WriteOptions::default())
    }

    /// Serialise as a `<qualityParameter/>` element under explicit options.
    ///
    /// # Errors
    ///
    /// As [`to_xml_string`](Self::to_xml_string).
    pub fn to_xml_string_with_options(
        &self,
        indentation_level: u32,
        options: &WriteOptions,
    ) -> Result<String> {
        let indent = indent(indentation_level)?;
        for field in self.fields() {
            check_text(field, Self::MAX_TEXT_BYTES)?;
        }
        for (label, value) in [
            ("name", &self.name),
            ("ID", &self.id),
            ("cvRef", &self.cv_ref),
            ("accession", &self.cv_acc),
        ] {
            if value.is_empty() {
                return Err(Error::MissingInformation(format!(
                    "qualityParameter requires a non-empty {label}"
                )));
            }
        }
        let mut out = indent.clone();
        out.push_str("<qualityParameter");
        attribute(&mut out, "name", &self.name)?;
        attribute(&mut out, "ID", &self.id)?;
        attribute(&mut out, "cvRef", &self.cv_ref)?;
        attribute(&mut out, "accession", &self.cv_acc)?;
        if !self.value.is_empty() {
            attribute(&mut out, "value", &self.value)?;
        }
        if !self.unit_ref.is_empty() {
            attribute(&mut out, options.unit_ref_attribute(), &self.unit_ref)?;
        }
        if !self.unit_acc.is_empty() {
            attribute(&mut out, options.unit_acc_attribute(), &self.unit_acc)?;
        }
        if !self.flag.is_empty() {
            let flag = if options.source_flag_literal {
                "true"
            } else {
                self.flag.as_str()
            };
            attribute(&mut out, "flag", flag)?;
        }
        out.push_str("/>\n");
        Ok(out)
    }

    fn fields(&self) -> [&String; 8] {
        [
            &self.name,
            &self.id,
            &self.value,
            &self.cv_ref,
            &self.cv_acc,
            &self.unit_ref,
            &self.unit_acc,
            &self.flag,
        ]
    }

    fn weight(&self) -> usize {
        self.fields().iter().map(|f| f.len() + 24).sum::<usize>() + 32
    }
}

/// Representation of an attachment: an opaque binary payload or an inline table.
///
/// The source's `id` field is documented `///< Name`, a copy-paste slip; it is
/// the `ID` attribute. `quality_ref` is the `qualityParameterRef` IDREF of the
/// parameter this attachment elaborates, and is empty when the attachment hangs
/// off the run or set itself.
///
/// # Comparison
///
/// As [`QualityParameter`]: the source compares `name` alone, `PartialEq` here
/// is structural and `Ord` is `name`-primary with deterministic ties.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Attachment {
    /// Name, written as the `name` attribute. Required by the reader.
    pub name: String,
    /// Identifier, written as the `ID` attribute. Required by the reader.
    pub id: String,
    /// Value, written as the `value` attribute and omitted when empty.
    pub value: String,
    /// Controlled vocabulary reference, written as `cvRef`. Required.
    pub cv_ref: String,
    /// Controlled vocabulary accession, written as `accession`. Required.
    pub cv_acc: String,
    /// Unit's controlled vocabulary reference, omitted when empty.
    pub unit_ref: String,
    /// Unit's controlled vocabulary accession, omitted when empty.
    pub unit_acc: String,
    /// Opaque `<binary>` content, conventionally base64. Neither the source nor
    /// this port decodes or validates it; it is carried verbatim.
    pub binary: String,
    /// `qualityParameterRef` IDREF, empty when attached to the run or set.
    pub quality_ref: String,
    /// Column types of the inline table, written space-delimited inside
    /// `<tableColumnTypes>`.
    pub col_types: Vec<String>,
    /// Cell values of the inline table, one `<tableRowValues>` element per row,
    /// in the column order of [`col_types`](Self::col_types).
    pub table_rows: Vec<Vec<String>>,
}

impl Attachment {
    /// Maximum column types in one table.
    pub const MAX_COLUMNS: usize = 100_000;
    /// Maximum rows in one table.
    pub const MAX_ROWS: usize = 4_000_000;
    /// Maximum total cells in one table.
    pub const MAX_TABLE_CELLS: usize = 16_000_000;
    /// Maximum bytes accepted in any single cell, scalar field or binary payload.
    pub const MAX_TEXT_BYTES: usize = 64 * 1024 * 1024;

    /// The source's `operator==`, `operator<` and `operator>` predicate: `name`
    /// equality alone.
    pub fn same_name(&self, other: &Self) -> bool {
        self.name == other.name
    }

    /// True when this attachment carries an inline table the writers can emit,
    /// that is a non-empty column list **and** at least one row.
    ///
    /// The source's `toXMLString` and `toCSVString` both require both halves; an
    /// attachment with column types but no rows is treated as having no table.
    pub fn has_table(&self) -> bool {
        !self.col_types.is_empty() && !self.table_rows.is_empty()
    }

    /// Render the inline table as delimiter-separated text: the column types,
    /// then one line per row, each line terminated by `\n`.
    ///
    /// Every occurrence of `separator` inside a column type or cell is replaced
    /// by `_`, or by `$` when `separator` is itself `_`, and each assembled line
    /// is trimmed of leading and trailing space, tab, carriage return and line
    /// feed — the source's `StringUtils::trimmed` character set. An attachment
    /// with no table renders as the empty string, as in the source.
    ///
    /// Rows are emitted with exactly the cells they hold: the source does not
    /// pad or check a row against the column count, so a short or long row
    /// yields a short or long line. `store`'s XML table has the same property.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `separator` is empty, because the source's
    /// `substitute` is a no-op for an empty needle and `concatenate` then runs
    /// all cells together into one unparsable field. [`Error::InvalidRange`]
    /// when the table exceeds [`MAX_TABLE_CELLS`](Self::MAX_TABLE_CELLS),
    /// [`MAX_COLUMNS`](Self::MAX_COLUMNS), [`MAX_ROWS`](Self::MAX_ROWS) or
    /// [`QcMLFile::MAX_OUTPUT_BYTES`].
    pub fn to_csv_string(&self, separator: &str) -> Result<String> {
        if separator.is_empty() {
            return Err(Error::InvalidValue(
                "qcML attachment CSV separator must not be empty".into(),
            ));
        }
        self.preflight_table()?;
        if !self.has_table() {
            return Ok(String::new());
        }
        let replacement = if separator == "_" { "$" } else { "_" };
        let mut out = String::new();
        let mut line = String::new();
        join_substituted(&mut line, &self.col_types, separator, replacement);
        out.push_str(line.trim_matches(trimmed));
        out.push('\n');
        for row in &self.table_rows {
            line.clear();
            join_substituted(&mut line, row, separator, replacement);
            out.push_str(line.trim_matches(trimmed));
            out.push('\n');
        }
        Ok(out)
    }

    /// Serialise as an `<attachment>` element indented with `indentation_level`
    /// tabs, under the native default options.
    ///
    /// A non-empty [`binary`](Self::binary) is written as a `<binary>` child and
    /// wins over any table, otherwise a complete table is written. The source
    /// silently returns the empty string when neither is present, and silently
    /// discards the table when both are; the native default refuses both cases
    /// so an attachment cannot vanish from a report without the caller saying
    /// so. [`WriteOptions::source`] selects the source behaviour.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when a required attribute is empty or the
    /// attachment carries neither a binary payload nor a complete table,
    /// [`Error::InvalidValue`] when both are present, when a column type or
    /// cell is empty or contains XML whitespace under the default options, or
    /// when any field holds a character XML 1.0 cannot represent, and
    /// [`Error::InvalidRange`] for the indentation, cell-count and byte
    /// ceilings.
    pub fn to_xml_string(&self, indentation_level: u32) -> Result<String> {
        self.to_xml_string_with_options(indentation_level, &WriteOptions::default())
    }

    /// Serialise as an `<attachment>` element under explicit options.
    ///
    /// # Errors
    ///
    /// As [`to_xml_string`](Self::to_xml_string). Under
    /// [`WriteOptions::source`] an attachment with neither payload yields
    /// `Ok("")` and one with both yields only its binary.
    pub fn to_xml_string_with_options(
        &self,
        indentation_level: u32,
        options: &WriteOptions,
    ) -> Result<String> {
        let indent = indent(indentation_level)?;
        self.preflight_table()?;
        for field in self.scalar_fields() {
            check_text(field, Self::MAX_TEXT_BYTES)?;
        }
        check_text(&self.binary, Self::MAX_TEXT_BYTES)?;
        let table = self.has_table();
        if !self.binary.is_empty() && table && !options.drop_unrepresentable {
            return Err(Error::InvalidValue(
                "qcML attachment holds both a binary payload and a table; the source writes only the binary".into(),
            ));
        }
        if self.binary.is_empty() && !table {
            if options.drop_unrepresentable {
                return Ok(String::new());
            }
            return Err(Error::MissingInformation(
                "qcML attachment has neither a binary payload nor a complete table".into(),
            ));
        }
        for (label, value) in [
            ("name", &self.name),
            ("ID", &self.id),
            ("cvRef", &self.cv_ref),
            ("accession", &self.cv_acc),
        ] {
            if value.is_empty() {
                return Err(Error::MissingInformation(format!(
                    "attachment requires a non-empty {label}"
                )));
            }
        }
        // Two spaces after the tag name reproduce the source literal
        // "<attachment " followed by " name=\"...\"".
        let mut out = indent.clone();
        out.push_str("<attachment ");
        attribute(&mut out, "name", &self.name)?;
        attribute(&mut out, "ID", &self.id)?;
        attribute(&mut out, "cvRef", &self.cv_ref)?;
        attribute(&mut out, "accession", &self.cv_acc)?;
        if !self.value.is_empty() {
            attribute(&mut out, "value", &self.value)?;
        }
        if !self.unit_ref.is_empty() {
            attribute(&mut out, options.unit_ref_attribute(), &self.unit_ref)?;
        }
        if !self.unit_acc.is_empty() {
            attribute(&mut out, options.unit_acc_attribute(), &self.unit_acc)?;
        }
        if !self.quality_ref.is_empty() {
            attribute(&mut out, "qualityParameterRef", &self.quality_ref)?;
        }
        if !self.binary.is_empty() {
            out.push_str(">\n");
            out.push_str(&indent);
            out.push_str("\t<binary>");
            push_text(&mut out, &self.binary)?;
            out.push_str("</binary>\n");
            out.push_str(&indent);
            out.push_str("</attachment>\n");
            return Ok(out);
        }
        out.push_str(">\n");
        // The source emits "<table>" with neither indentation nor a newline and
        // closes it the same way; the layout is preserved verbatim.
        out.push_str("<table>");
        out.push_str(&indent);
        out.push_str("\t<tableColumnTypes>");
        push_cells(&mut out, &self.col_types, options, "column type", true)?;
        out.push_str("</tableColumnTypes>\n");
        for row in &self.table_rows {
            out.push_str(&indent);
            out.push_str("\t<tableRowValues>");
            push_cells(&mut out, row, options, "table cell", false)?;
            out.push_str("</tableRowValues>\n");
        }
        out.push_str("</table>");
        out.push_str(&indent);
        out.push_str("</attachment>\n");
        Ok(out)
    }

    fn scalar_fields(&self) -> [&String; 8] {
        [
            &self.name,
            &self.id,
            &self.value,
            &self.cv_ref,
            &self.cv_acc,
            &self.unit_ref,
            &self.unit_acc,
            &self.quality_ref,
        ]
    }

    fn preflight_table(&self) -> Result<()> {
        if self.col_types.len() > Self::MAX_COLUMNS {
            return Err(limit("qcML attachment column limit exceeded"));
        }
        if self.table_rows.len() > Self::MAX_ROWS {
            return Err(limit("qcML attachment row limit exceeded"));
        }
        let mut cells = self.col_types.len();
        let mut bytes = self.binary.len();
        for cell in &self.col_types {
            bytes = bytes
                .checked_add(cell.len().saturating_add(8))
                .ok_or_else(|| limit("qcML attachment byte limit exceeded"))?;
        }
        for row in &self.table_rows {
            cells = cells
                .checked_add(row.len())
                .ok_or_else(|| limit("qcML attachment cell limit exceeded"))?;
            if cells > Self::MAX_TABLE_CELLS {
                return Err(limit("qcML attachment cell limit exceeded"));
            }
            bytes = bytes
                .checked_add(32)
                .ok_or_else(|| limit("qcML attachment byte limit exceeded"))?;
            for cell in row {
                bytes = bytes
                    .checked_add(cell.len().saturating_add(8))
                    .ok_or_else(|| limit("qcML attachment byte limit exceeded"))?;
            }
            if bytes > QcMLFile::MAX_OUTPUT_BYTES {
                return Err(limit("qcML attachment byte limit exceeded"));
            }
        }
        if cells > Self::MAX_TABLE_CELLS || bytes > QcMLFile::MAX_OUTPUT_BYTES {
            return Err(limit("qcML attachment cell or byte limit exceeded"));
        }
        Ok(())
    }

    fn weight(&self) -> usize {
        let mut bytes = self
            .scalar_fields()
            .iter()
            .map(|f| f.len() + 24)
            .sum::<usize>()
            + self.binary.len()
            + 64;
        for cell in &self.col_types {
            bytes = bytes.saturating_add(cell.len().saturating_add(8));
        }
        for row in &self.table_rows {
            bytes = bytes.saturating_add(32);
            for cell in row {
                bytes = bytes.saturating_add(cell.len().saturating_add(8));
            }
        }
        bytes
    }
}

/// How `merge` collapses parameters and attachments that compare equal.
///
/// The source appends the addendum's records, `std::sort`s and then applies
/// `std::unique` with comparators that look at `name` alone, so two parameters
/// that share a name but differ in identifier, value or unit are silently
/// reduced to one. The native default deduplicates on full structural equality,
/// which removes only genuine duplicates; [`source`](MergeOptions::source)
/// selects the source's name-only collapse.
#[derive(Clone, Copy, Debug, Default)]
pub struct MergeOptions {
    /// Collapse records that merely share a `name`, as the source does.
    pub collapse_by_name: bool,
}

impl MergeOptions {
    /// The source `QcMLFile::merge` behaviour: collapse by `name` alone.
    pub fn source() -> Self {
        Self {
            collapse_by_name: true,
        }
    }
}

/// A report stylesheet injected into the stored document.
///
/// The source's `store` looks for `XSL/QcML_report_sheet.xsl` in the OpenMS
/// share directory and, when it finds it, writes an `xml-stylesheet` processing
/// instruction, a DOCTYPE declaring the `id` attribute of `xsl:stylesheet`, and
/// the stylesheet body just before `</qcML>`, so that a browser renders the
/// report as HTML. When the file is absent it warns "No qcml stylesheet found,
/// result will not be viewable in a browser!" and omits all three.
///
/// This crate ships no XSL resource, so the default is the source's
/// stylesheet-absent path. A caller that has the stylesheet supplies it here.
#[derive(Clone, Debug)]
pub struct Stylesheet {
    /// `id` of the embedded `xsl:stylesheet`, referenced by the processing
    /// instruction as `href="#<id>"`. The source hard-codes
    /// `openms-qc-stylesheet`.
    pub id: String,
    /// Stylesheet body, written verbatim between `</cvList>` and `</qcML>`. It
    /// must already have its own XML declaration removed;
    /// [`from_file_text`](Self::from_file_text) does that. A body holding
    /// `</qcML>`, an XML 1.0-invalid character, or more than
    /// [`MAX_BYTES`](Self::MAX_BYTES) bytes is rejected at write time;
    /// everything else is the caller's responsibility, because the source does
    /// not parse the stylesheet either.
    pub xslt: String,
}

impl Stylesheet {
    /// Maximum stylesheet body in bytes.
    pub const MAX_BYTES: usize = 16 * 1024 * 1024;

    /// Build a stylesheet from the raw text of an `.xsl` file, dropping
    /// everything up to and including the first line feed.
    ///
    /// This is exactly the source's `xslt.erase(0, xslt.find('\n') + 1)`, whose
    /// purpose is to remove the stylesheet's own `<?xml ... ?>` declaration so
    /// it can be embedded. Text with no line feed is kept whole, because
    /// `std::string::npos + 1` is `0` and the erase is then a no-op.
    pub fn from_file_text(id: &str, text: &str) -> Self {
        let xslt = match text.split_once('\n') {
            Some((_, rest)) => rest,
            None => text,
        };
        Self {
            id: id.to_owned(),
            xslt: xslt.to_owned(),
        }
    }
}

/// How a qcML document is written.
///
/// The native default produces a document this port can read back without loss.
/// [`source`](WriteOptions::source) reproduces the pinned C++ `store` byte for
/// byte, including its three round-trip defects, for a caller that must match
/// OpenMS output.
#[derive(Clone, Debug, Default)]
pub struct WriteOptions {
    /// Write the unit attributes as `unitRef` and `unitAcc`, the spellings the
    /// source's writer emits, instead of the `unitCvRef` and `unitAccession`
    /// spellings its reader parses.
    ///
    /// With the source spellings a stored unit is lost on reload: `store`
    /// writes `unitAcc`/`unitRef` while `onStartElement` reads
    /// `unitAccession`/`unitCvRef`. This port's reader accepts both, so the
    /// setting only decides which spelling is produced.
    pub source_unit_attributes: bool,
    /// Write `flag="true"` for any non-empty [`QualityParameter::flag`] instead
    /// of the value itself, as the source does, discarding the value.
    pub source_flag_literal: bool,
    /// Reproduce the source's table text handling: substitute `' '` with `'_'`
    /// in column types only, and write row cells unchanged.
    ///
    /// The source intends to substitute in rows too — it builds a substituted
    /// copy of each row — but then concatenates the *original* row, so the
    /// substitution is discarded and a cell containing a space silently becomes
    /// two cells on reload. The native default instead refuses a column type or
    /// cell that is empty or contains XML whitespace, because the format is
    /// space-delimited and cannot represent either.
    pub source_table_text: bool,
    /// Silently drop an attachment that has neither a binary payload nor a
    /// complete table, and drop the table of one that has both, as the source's
    /// `toXMLString` does by returning the empty string.
    pub drop_unrepresentable: bool,
    /// Declare `encoding="ISO-8859-1"` as the source does, instead of `UTF-8`.
    ///
    /// The source writes `std::string` bytes unchanged under that declaration,
    /// so any non-ASCII content produces a document whose declaration lies about
    /// its bytes. This port refuses non-ASCII content under this setting rather
    /// than emitting such a document.
    pub source_encoding_declaration: bool,
    /// Report stylesheet to inject, and with it the `xml-stylesheet` processing
    /// instruction and the DOCTYPE the source emits alongside.
    pub stylesheet: Option<Stylesheet>,
}

impl WriteOptions {
    /// Reproduce the pinned C++ `QcMLFile::store` output, defects included.
    ///
    /// Byte-identical to the source for content that needs no XML escaping,
    /// which is every accession, identifier and number a qcML report normally
    /// holds. Values containing `&`, `<`, `>` or `"` differ, because the source
    /// concatenates them verbatim and produces a document no XML parser accepts,
    /// while this port always escapes.
    ///
    /// `stylesheet` stays `None`: the source injects one only when it finds the
    /// XSL file, and no such resource is bundled here.
    pub fn source() -> Self {
        Self {
            source_unit_attributes: true,
            source_flag_literal: true,
            source_table_text: true,
            drop_unrepresentable: true,
            source_encoding_declaration: true,
            stylesheet: None,
        }
    }

    fn unit_ref_attribute(&self) -> &'static str {
        if self.source_unit_attributes {
            "unitRef"
        } else {
            "unitCvRef"
        }
    }

    fn unit_acc_attribute(&self) -> &'static str {
        if self.source_unit_attributes {
            "unitAcc"
        } else {
            "unitAccession"
        }
    }
}

/// A qcML document: quality parameters and attachments per run and per set.
///
/// The source keeps seven `protected` maps and exposes them only through the
/// registration, lookup and export methods; a Rust type has no protected access,
/// so the collections are reachable through the borrowing accessors documented
/// as native additions below. Every map is a [`BTreeMap`], matching the source's
/// `std::map` iteration order, which `store` and `getRunIDs` both depend on.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct QcMLFile {
    run_qps: BTreeMap<String, Vec<QualityParameter>>,
    run_ats: BTreeMap<String, Vec<Attachment>>,
    set_qps: BTreeMap<String, Vec<QualityParameter>>,
    set_ats: BTreeMap<String, Vec<Attachment>>,
    set_member_names: BTreeMap<String, BTreeSet<String>>,
    run_name_ids: BTreeMap<String, String>,
    set_name_ids: BTreeMap<String, String>,
}

impl QcMLFile {
    /// Maximum runs, and separately maximum sets, in one document.
    pub const MAX_ENTRIES: usize = 1_000_000;
    /// Maximum quality parameters held by one run or set.
    pub const MAX_PARAMETERS_PER_ENTRY: usize = 1_000_000;
    /// Maximum attachments held by one run or set.
    pub const MAX_ATTACHMENTS_PER_ENTRY: usize = 1_000_000;
    /// Maximum member names recorded for one set.
    pub const MAX_SET_MEMBERS: usize = 1_000_000;
    /// Maximum bytes any one serialisation or export may produce.
    pub const MAX_OUTPUT_BYTES: usize = 1024 * 1024 * 1024;

    /// An empty document, as the source's default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// True when the document holds no run and no set at all.
    ///
    /// Native addition; the source offers no such predicate.
    pub fn is_empty(&self) -> bool {
        self.run_qps.is_empty()
            && self.run_ats.is_empty()
            && self.set_qps.is_empty()
            && self.set_ats.is_empty()
    }

    /// Register a run under `id` and map `name` to it.
    ///
    /// As in the source this **resets** the run: its parameter and attachment
    /// lists are replaced with empty ones, so registering an existing `id` again
    /// discards what it held. Registering a `name` already mapped to another
    /// `id` repoints the mapping and leaves the previous run reachable only by
    /// its identifier.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `id` or `name` is empty — the source accepts
    /// both and produces a run that `store` writes with `ID=""` and that no
    /// name lookup can distinguish — and [`Error::InvalidRange`] when the
    /// document already holds [`MAX_ENTRIES`](Self::MAX_ENTRIES) runs or either
    /// string exceeds [`QualityParameter::MAX_TEXT_BYTES`].
    pub fn register_run(&mut self, id: &str, name: &str) -> Result<()> {
        self.register(id, name, true)
    }

    /// Register a set under `id`, map `name` to it and record its member names.
    ///
    /// `member_names` are the `MS:1000577` values of the runs the set
    /// aggregates. Resetting and name-repointing behave as
    /// [`register_run`](Self::register_run).
    ///
    /// # Errors
    ///
    /// As [`register_run`](Self::register_run), plus [`Error::InvalidRange`]
    /// when `member_names` exceeds [`MAX_SET_MEMBERS`](Self::MAX_SET_MEMBERS).
    pub fn register_set(
        &mut self,
        id: &str,
        name: &str,
        member_names: &BTreeSet<String>,
    ) -> Result<()> {
        if member_names.len() > Self::MAX_SET_MEMBERS {
            return Err(limit("qcML set member limit exceeded"));
        }
        for member in member_names {
            check_text(member, QualityParameter::MAX_TEXT_BYTES)?;
        }
        self.register(id, name, false)?;
        self.set_member_names
            .insert(id.to_owned(), member_names.clone());
        Ok(())
    }

    fn register(&mut self, id: &str, name: &str, run: bool) -> Result<()> {
        let kind = if run { "run" } else { "set" };
        if id.is_empty() {
            return Err(Error::InvalidValue(format!("qcML {kind} needs an ID")));
        }
        if name.is_empty() {
            return Err(Error::InvalidValue(format!("qcML {kind} needs a name")));
        }
        check_text(id, QualityParameter::MAX_TEXT_BYTES)?;
        check_text(name, QualityParameter::MAX_TEXT_BYTES)?;
        let (qps, ats, names) = if run {
            (&mut self.run_qps, &mut self.run_ats, &mut self.run_name_ids)
        } else {
            (&mut self.set_qps, &mut self.set_ats, &mut self.set_name_ids)
        };
        if !qps.contains_key(id) && qps.len() >= Self::MAX_ENTRIES {
            return Err(limit("qcML entry limit exceeded"));
        }
        qps.insert(id.to_owned(), Vec::new());
        ats.insert(id.to_owned(), Vec::new());
        names.insert(name.to_owned(), id.to_owned());
        Ok(())
    }

    /// Append a quality parameter to the run named or identified by `r`.
    ///
    /// `r` is matched against the run identifiers first and against the
    /// registered run names second, so either addresses the run.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when `r` is neither. The source silently
    /// discards the parameter in that case — its own comment reads "TODO warn
    /// that run has to be registered!" — which loses data with no diagnostic,
    /// so this port refuses. [`Error::InvalidRange`] when the run already holds
    /// [`MAX_PARAMETERS_PER_ENTRY`](Self::MAX_PARAMETERS_PER_ENTRY) parameters
    /// or a field exceeds [`QualityParameter::MAX_TEXT_BYTES`].
    pub fn add_run_quality_parameter(&mut self, r: &str, qp: QualityParameter) -> Result<()> {
        let id = Self::resolve(&self.run_qps, &self.run_name_ids, r)
            .ok_or_else(|| Error::MissingInformation(format!("no qcML run {r:?} is registered")))?
            .to_owned();
        Self::push_parameter(&mut self.run_qps, &id, qp)
    }

    /// Append a quality parameter to the set named or identified by `r`.
    ///
    /// # Errors
    ///
    /// As [`add_run_quality_parameter`](Self::add_run_quality_parameter), for
    /// sets.
    pub fn add_set_quality_parameter(&mut self, r: &str, qp: QualityParameter) -> Result<()> {
        let id = Self::resolve(&self.set_qps, &self.set_name_ids, r)
            .ok_or_else(|| Error::MissingInformation(format!("no qcML set {r:?} is registered")))?
            .to_owned();
        Self::push_parameter(&mut self.set_qps, &id, qp)
    }

    /// Append an attachment to the run identified by `r`, creating the run's
    /// attachment list if it has none.
    ///
    /// Unlike the parameter setters this does **not** require registration and
    /// does not consult the name map: the source's `addRunAttachment` indexes
    /// `runQualityAts_[run_id]` directly, deliberately permitting an attachment
    /// with no quality parameter. A run known only this way is invisible to
    /// [`exists_run`](Self::exists_run) and absent from
    /// [`run_ids`](Self::run_ids), because both read the parameter map, but
    /// `store` still writes it: the stored run set is the union of both maps'
    /// keys.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `r` is empty, and [`Error::InvalidRange`]
    /// for the per-entry attachment ceiling, the table ceilings and the field
    /// byte ceiling.
    pub fn add_run_attachment(&mut self, r: &str, at: Attachment) -> Result<()> {
        Self::push_attachment(&mut self.run_ats, r, at)
    }

    /// Append an attachment to the set identified by `r`, creating the set's
    /// attachment list if it has none.
    ///
    /// # Errors
    ///
    /// As [`add_run_attachment`](Self::add_run_attachment), for sets.
    pub fn add_set_attachment(&mut self, r: &str, at: Attachment) -> Result<()> {
        Self::push_attachment(&mut self.set_ats, r, at)
    }

    fn push_parameter(
        map: &mut BTreeMap<String, Vec<QualityParameter>>,
        id: &str,
        qp: QualityParameter,
    ) -> Result<()> {
        check_parameter(&qp)?;
        let list = map.entry(id.to_owned()).or_default();
        if list.len() >= Self::MAX_PARAMETERS_PER_ENTRY {
            return Err(limit("qcML quality parameter limit exceeded"));
        }
        list.push(qp);
        Ok(())
    }

    fn push_attachment(
        map: &mut BTreeMap<String, Vec<Attachment>>,
        r: &str,
        at: Attachment,
    ) -> Result<()> {
        if r.is_empty() {
            return Err(Error::InvalidValue(
                "qcML attachment needs a run or set ID".into(),
            ));
        }
        check_text(r, QualityParameter::MAX_TEXT_BYTES)?;
        check_attachment(&at)?;
        if map.len() >= Self::MAX_ENTRIES && !map.contains_key(r) {
            return Err(limit("qcML entry limit exceeded"));
        }
        let list = map.entry(r.to_owned()).or_default();
        if list.len() >= Self::MAX_ATTACHMENTS_PER_ENTRY {
            return Err(limit("qcML attachment limit exceeded"));
        }
        list.push(at);
        Ok(())
    }

    fn resolve<'a, T>(
        data: &'a BTreeMap<String, T>,
        names: &'a BTreeMap<String, String>,
        key: &'a str,
    ) -> Option<&'a str> {
        if data.contains_key(key) {
            return Some(key);
        }
        names.get(key).map(String::as_str)
    }

    /// The registered run identifiers, in lexical order.
    ///
    /// The source's `getRunIDs` fills a caller's vector with the keys of the run
    /// **parameter** map, so a run known only through
    /// [`add_run_attachment`](Self::add_run_attachment) is not listed.
    pub fn run_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.run_qps.keys().map(String::as_str)
    }

    /// The registered run names, in lexical order.
    ///
    /// The source's `getRunNames` fills a caller's vector with the keys of the
    /// name-to-identifier map. That map holds one entry per *name*, so two runs
    /// registered under the same name appear once, and a run whose name was
    /// repointed to another identifier still appears.
    pub fn run_names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.run_name_ids.keys().map(String::as_str)
    }

    /// The registered set identifiers, in lexical order. Native addition,
    /// symmetric to [`run_ids`](Self::run_ids).
    pub fn set_ids(&self) -> impl ExactSizeIterator<Item = &str> {
        self.set_qps.keys().map(String::as_str)
    }

    /// The registered set names, in lexical order. Native addition, symmetric
    /// to [`run_names`](Self::run_names).
    pub fn set_names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.set_name_ids.keys().map(String::as_str)
    }

    /// The identifier a run name maps to. Native addition; the source reaches
    /// its `run_Name_ID_map_` only from inside its own methods.
    pub fn run_id_for_name(&self, name: &str) -> Option<&str> {
        self.run_name_ids.get(name).map(String::as_str)
    }

    /// The identifier a set name maps to. Native addition.
    pub fn set_id_for_name(&self, name: &str) -> Option<&str> {
        self.set_name_ids.get(name).map(String::as_str)
    }

    /// True when `id` is a registered run identifier.
    ///
    /// The source's `existsRun(filename, checkname = false)`. Its comment "NO,
    /// do not!: permit AT without a QP" records the deliberate choice to consult
    /// only the parameter map.
    pub fn exists_run(&self, id: &str) -> bool {
        self.run_qps.contains_key(id)
    }

    /// True when `filename` is a registered run identifier **or** a registered
    /// run name, checked in that order.
    ///
    /// The source's `existsRun(filename, true)`; a Rust boolean parameter is
    /// replaced by the second name, per the port's overload rule.
    pub fn exists_run_or_name(&self, filename: &str) -> bool {
        self.exists_run(filename) || self.run_name_ids.contains_key(filename)
    }

    /// True when `id` is a registered set identifier.
    /// The source's `existsSet(filename, checkname = false)`.
    pub fn exists_set(&self, id: &str) -> bool {
        self.set_qps.contains_key(id)
    }

    /// True when `filename` is a registered set identifier **or** a registered
    /// set name. The source's `existsSet(filename, true)`.
    pub fn exists_set_or_name(&self, filename: &str) -> bool {
        self.exists_set(filename) || self.set_name_ids.contains_key(filename)
    }

    /// Identifiers of the run's quality parameters whose **accession** is
    /// `qp_accession`, in insertion order.
    ///
    /// `filename` is a run identifier or, failing that, a run name. The source
    /// clears and fills a caller's vector; an empty result means "not found".
    /// Note that the match is on `cvAcc`, not on `name`, despite the parameter
    /// being called `qpname`.
    pub fn exists_run_quality_parameter(&self, filename: &str, qp_accession: &str) -> Vec<String> {
        Self::parameter_ids(&self.run_qps, &self.run_name_ids, filename, qp_accession)
    }

    /// Identifiers of the set's quality parameters whose **accession** is
    /// `qp_accession`, in insertion order.
    pub fn exists_set_quality_parameter(&self, filename: &str, qp_accession: &str) -> Vec<String> {
        Self::parameter_ids(&self.set_qps, &self.set_name_ids, filename, qp_accession)
    }

    fn parameter_ids(
        data: &BTreeMap<String, Vec<QualityParameter>>,
        names: &BTreeMap<String, String>,
        filename: &str,
        accession: &str,
    ) -> Vec<String> {
        let Some(id) = Self::resolve(data, names, filename) else {
            return Vec::new();
        };
        let Some(list) = data.get(id) else {
            return Vec::new();
        };
        list.iter()
            .filter(|qp| qp.cv_acc == accession)
            .map(|qp| qp.id.clone())
            .collect()
    }

    /// The quality parameters of the run identified by `id`, or an empty slice.
    /// Native addition: the source's map is `protected`.
    pub fn run_quality_parameters(&self, id: &str) -> &[QualityParameter] {
        self.run_qps.get(id).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The attachments of the run identified by `id`, or an empty slice.
    /// Native addition.
    pub fn run_attachments(&self, id: &str) -> &[Attachment] {
        self.run_ats.get(id).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The quality parameters of the set identified by `id`, or an empty slice.
    /// Native addition.
    pub fn set_quality_parameters(&self, id: &str) -> &[QualityParameter] {
        self.set_qps.get(id).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The attachments of the set identified by `id`, or an empty slice.
    /// Native addition.
    pub fn set_attachments(&self, id: &str) -> &[Attachment] {
        self.set_ats.get(id).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The member names recorded for the set identified by `id`, in lexical
    /// order. Native addition.
    pub fn set_members(&self, id: &str) -> impl Iterator<Item = &str> {
        self.set_member_names
            .get(id)
            .into_iter()
            .flat_map(|names| names.iter().map(String::as_str))
    }

    /// Remove attachments of the run **and** set identified by `r` whose
    /// `quality_ref` is one of `ids`.
    ///
    /// `name` restricts the removal to attachments with that `name`; `None`
    /// removes every attachment referencing one of `ids`, which is the source's
    /// default empty `at` argument.
    ///
    /// The source reaches its maps with `operator[]`, which inserts an empty
    /// list for an unknown `r` and so silently registers a phantom run and set
    /// on every call. This port touches only lists that already exist.
    pub fn remove_attachments_by_quality_ref(
        &mut self,
        r: &str,
        ids: &[String],
        name: Option<&str>,
    ) {
        for map in [&mut self.run_ats, &mut self.set_ats] {
            if let Some(list) = map.get_mut(r) {
                list.retain(|at| {
                    !ids.iter().any(|id| &at.quality_ref == id)
                        || name.is_some_and(|n| at.name != n)
                });
            }
        }
    }

    /// Remove attachments whose **accession** is `accession` from the run
    /// identified by `r`, and from the set identified by `r`.
    ///
    /// The run's attachments are touched only when `r` is a registered run
    /// identifier and the set's only when `r` is a registered set identifier,
    /// reproducing the source's `existsRun`/`existsSet` guards. Names are not
    /// consulted.
    pub fn remove_attachments_by_accession(&mut self, r: &str, accession: &str) {
        if self.exists_run(r) {
            if let Some(list) = self.run_ats.get_mut(r) {
                list.retain(|at| at.cv_acc != accession);
            }
        }
        if self.exists_set(r) {
            if let Some(list) = self.set_ats.get_mut(r) {
                list.retain(|at| at.cv_acc != accession);
            }
        }
    }

    /// Remove attachments whose accession is `accession` from every run that has
    /// an attachment list.
    ///
    /// The source's comment claims "from all runs/sets", but the loop iterates
    /// the run attachment map only and calls the per-entry removal above, so a
    /// set's attachments are reached only when its identifier is also a run
    /// identifier with an attachment list. That behaviour is preserved; a caller
    /// that means every set must iterate [`set_ids`](Self::set_ids) itself.
    pub fn remove_all_attachments(&mut self, accession: &str) {
        let ids: Vec<String> = self.run_ats.keys().cloned().collect();
        for id in ids {
            self.remove_attachments_by_accession(&id, accession);
        }
    }

    /// Remove the quality parameters of the run and set identified by `r` whose
    /// identifier is one of `ids`, together with the attachments referencing
    /// them.
    ///
    /// The source removes the attachments first, by delegating to the
    /// quality-reference removal above with no name restriction, and this port
    /// keeps that order. As there, no phantom entry is created for an unknown
    /// `r`.
    pub fn remove_quality_parameters(&mut self, r: &str, ids: &[String]) {
        self.remove_attachments_by_quality_ref(r, ids, None);
        for map in [&mut self.run_qps, &mut self.set_qps] {
            if let Some(list) = map.get_mut(r) {
                list.retain(|qp| !ids.iter().any(|id| &qp.id == id));
            }
        }
    }

    /// Merge `addendum` into this document, optionally recording every merged
    /// run as a member of the set named `set_name`.
    ///
    /// Per run and per set the addendum's parameters and attachments are
    /// appended, the combined list is sorted and duplicates are removed as
    /// `options` directs. Set membership from the addendum is added only for
    /// sets this document does not already know: the source uses
    /// `std::map::insert`, which keeps the existing value, so a set present in
    /// both retains this document's member list.
    ///
    /// Ordering is deterministic here. The source sorts with a `name`-only
    /// comparator through `std::sort`, which is not stable, so the relative order
    /// of equally named records is unspecified; this port's [`Ord`] breaks those
    /// ties on the remaining fields.
    ///
    /// The merge is atomic: it is assembled in a temporary and committed only
    /// once every ceiling has held, so a failure leaves this document untouched.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the merged document would exceed
    /// [`MAX_ENTRIES`](Self::MAX_ENTRIES),
    /// [`MAX_PARAMETERS_PER_ENTRY`](Self::MAX_PARAMETERS_PER_ENTRY),
    /// [`MAX_ATTACHMENTS_PER_ENTRY`](Self::MAX_ATTACHMENTS_PER_ENTRY),
    /// [`MAX_SET_MEMBERS`](Self::MAX_SET_MEMBERS) or
    /// [`MAX_OUTPUT_BYTES`](Self::MAX_OUTPUT_BYTES) of payload, and
    /// [`Error::InvalidValue`] when `set_name` is empty. The source takes an
    /// empty `setname` to mean "do not build a set"; pass `None` for that.
    pub fn merge(
        &mut self,
        addendum: &Self,
        set_name: Option<&str>,
        options: &MergeOptions,
    ) -> Result<()> {
        if set_name.is_some_and(str::is_empty) {
            return Err(Error::InvalidValue(
                "qcML merge set name must not be empty; pass None instead".into(),
            ));
        }
        preflight_merge(self, addendum)?;
        let mut next = self.clone();
        for (id, qps) in &addendum.run_qps {
            let list = next.run_qps.entry(id.clone()).or_default();
            list.extend(qps.iter().cloned());
            collapse(list, options, QualityParameter::same_name);
            if let Some(name) = set_name {
                next.set_member_names
                    .entry(name.to_owned())
                    .or_default()
                    .insert(id.clone());
            }
        }
        for (id, ats) in &addendum.run_ats {
            let list = next.run_ats.entry(id.clone()).or_default();
            list.extend(ats.iter().cloned());
            collapse(list, options, Attachment::same_name);
            if let Some(name) = set_name {
                next.set_member_names
                    .entry(name.to_owned())
                    .or_default()
                    .insert(id.clone());
            }
        }
        for (id, members) in &addendum.set_member_names {
            next.set_member_names
                .entry(id.clone())
                .or_insert_with(|| members.clone());
        }
        for (id, qps) in &addendum.set_qps {
            let list = next.set_qps.entry(id.clone()).or_default();
            list.extend(qps.iter().cloned());
            collapse(list, options, QualityParameter::same_name);
        }
        for (id, ats) in &addendum.set_ats {
            let list = next.set_ats.entry(id.clone()).or_default();
            list.extend(ats.iter().cloned());
            collapse(list, options, Attachment::same_name);
        }
        for map in [&next.run_qps, &next.set_qps] {
            if map.len() > Self::MAX_ENTRIES {
                return Err(limit("qcML entry limit exceeded"));
            }
            for list in map.values() {
                if list.len() > Self::MAX_PARAMETERS_PER_ENTRY {
                    return Err(limit("qcML quality parameter limit exceeded"));
                }
            }
        }
        for map in [&next.run_ats, &next.set_ats] {
            if map.len() > Self::MAX_ENTRIES {
                return Err(limit("qcML entry limit exceeded"));
            }
            for list in map.values() {
                if list.len() > Self::MAX_ATTACHMENTS_PER_ENTRY {
                    return Err(limit("qcML attachment limit exceeded"));
                }
            }
        }
        for members in next.set_member_names.values() {
            if members.len() > Self::MAX_SET_MEMBERS {
                return Err(limit("qcML set member limit exceeded"));
            }
        }
        *self = next;
        Ok(())
    }

    /// Values of the parameters with accession `qp_accession` held by the runs
    /// that are members of the set named `set_name`, in member order.
    ///
    /// The source's `collectSetParameter` appends to a caller's vector without
    /// clearing it and, being non-`const`, indexes its maps with `operator[]`,
    /// so asking about an unknown set or member *creates* an empty set and an
    /// empty run as a side effect. This port takes `&self` and returns the
    /// values, so a lookup mutates nothing.
    ///
    /// Members are matched against run **identifiers**, which is what the source
    /// does: `runQualityQPs_[*it]`. A set whose members were recorded as
    /// `MS:1000577` names — which is how the reader and `merge` record them when
    /// the name differs from the identifier — therefore yields nothing.
    pub fn collect_set_parameter(&self, set_name: &str, qp_accession: &str) -> Vec<String> {
        let mut out = Vec::new();
        let Some(members) = self.set_member_names.get(set_name) else {
            return out;
        };
        for member in members {
            for qp in self.run_quality_parameters(member) {
                if qp.cv_acc == qp_accession {
                    out.push(qp.value.clone());
                }
            }
        }
        out
    }

    /// The first matching attachment of the run or set `filename`, rendered as
    /// tab-separated text, or `None` when there is none.
    ///
    /// `qp_name` is matched against each attachment's `name` **or** `accession`.
    /// Runs are searched first, by identifier and then by name, then sets the
    /// same way. The source returns the empty string when nothing matches, which
    /// is indistinguishable from a matched attachment that has no table; `None`
    /// here means "no match" and `Some("")` means "matched, but no table".
    ///
    /// # Errors
    ///
    /// As [`Attachment::to_csv_string`], for the matched attachment.
    pub fn export_attachment(&self, filename: &str, qp_name: &str) -> Result<Option<String>> {
        for (data, names) in [
            (&self.run_ats, &self.run_name_ids),
            (&self.set_ats, &self.set_name_ids),
        ] {
            let Some(id) = Self::resolve(data, names, filename) else {
                continue;
            };
            let Some(list) = data.get(id) else { continue };
            for at in list {
                if at.name == qp_name || at.cv_acc == qp_name {
                    return at.to_csv_string("\t").map(Some);
                }
            }
        }
        Ok(None)
    }

    /// The value of the first matching quality parameter of the run or set
    /// `filename`, or `None`.
    ///
    /// The source matches **runs on accession** and **sets on name**. That
    /// asymmetry is not documented upstream and looks accidental, but it decides
    /// which lookups succeed, so it is reproduced exactly: a set parameter is
    /// found by its `name`, a run parameter by its `cvAcc`. Runs are searched
    /// first, by identifier then by name, then sets.
    ///
    /// The source returns the literal [`NOT_FOUND`] when nothing matches;
    /// [`export_quality_parameters`](Self::export_quality_parameters) reproduces
    /// that in its joined output.
    pub fn export_quality_parameter(&self, filename: &str, qp_name: &str) -> Option<&str> {
        if let Some(id) = Self::resolve(&self.run_qps, &self.run_name_ids, filename) {
            if let Some(list) = self.run_qps.get(id) {
                for qp in list {
                    if qp.cv_acc == qp_name {
                        return Some(qp.value.as_str());
                    }
                }
            }
        }
        let id = Self::resolve(&self.set_qps, &self.set_name_ids, filename)?;
        let list = self.set_qps.get(id)?;
        for qp in list {
            if qp.name == qp_name {
                return Some(qp.value.as_str());
            }
        }
        None
    }

    /// The values of several quality parameters of `filename`, each followed by
    /// a comma.
    ///
    /// Reproduces the source's `exportQPs` exactly, including the trailing
    /// comma and the [`NOT_FOUND`] placeholder for a parameter that is absent,
    /// so the field count is stable and the output remains parsable.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the joined text would exceed
    /// [`MAX_OUTPUT_BYTES`](Self::MAX_OUTPUT_BYTES).
    pub fn export_quality_parameters(&self, filename: &str, qp_names: &[String]) -> Result<String> {
        let mut out = String::new();
        for name in qp_names {
            let value = self
                .export_quality_parameter(filename, name)
                .unwrap_or(NOT_FOUND);
            if out.len().saturating_add(value.len()) + 1 > Self::MAX_OUTPUT_BYTES {
                return Err(limit("qcML export byte limit exceeded"));
            }
            out.push_str(value);
            out.push(',');
        }
        Ok(out)
    }

    /// Identification statistics of the set `filename` as a tab-separated table,
    /// or `None` when the set holds none.
    ///
    /// Two rows are built: `id` from the parameters with accession `QC:0000043`
    /// through `QC:0000047`, and `ms2` from `QC:0000053` through `QC:0000057`.
    /// Each parameter contributes one column, named by the part of its `name`
    /// before the first space, with its `value` as the cell.
    ///
    /// # Warning
    ///
    /// The result is usually misaligned, and that is the source's behaviour, not
    /// this port's. [`map_to_csv`] takes the column names from the *first* row
    /// alone, and the `id` and `ms2` accessions have different CV names, so the
    /// `ms2` row almost always lacks every column the header declares and is
    /// emitted as its label alone. Read the two groups separately if you need a
    /// well-formed table.
    ///
    /// # Errors
    ///
    /// As [`map_to_csv`].
    pub fn export_id_stats(&self, filename: &str) -> Result<Option<String>> {
        let Some(id) = Self::resolve(&self.set_qps, &self.set_name_ids, filename) else {
            return Ok(None);
        };
        let Some(list) = self.set_qps.get(id) else {
            return Ok(None);
        };
        let mut table: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
        for qp in list {
            let row = match qp.cv_acc.as_str() {
                "QC:0000043" | "QC:0000044" | "QC:0000045" | "QC:0000046" | "QC:0000047" => "id",
                "QC:0000053" | "QC:0000054" | "QC:0000055" | "QC:0000056" | "QC:0000057" => "ms2",
                _ => continue,
            };
            let column = qp.name.split(' ').next().unwrap_or("").to_owned();
            table
                .entry(row.to_owned())
                .or_default()
                .insert(column, qp.value.clone());
        }
        if table.is_empty() {
            return Ok(None);
        }
        map_to_csv(&table, "\t").map(Some)
    }

    /// Serialise the whole document, under the native default options.
    ///
    /// # Errors
    ///
    /// As [`to_xml_string_with_options`](Self::to_xml_string_with_options).
    pub fn to_xml_string(&self) -> Result<String> {
        self.to_xml_string_with_options(&WriteOptions::default())
    }

    /// Serialise the whole document under explicit options.
    ///
    /// Runs are written first, then sets, each in lexical identifier order over
    /// the union of the parameter and attachment map keys — the source builds
    /// that union in a `std::set`. Inside a `<setQuality>` the recorded members
    /// come first, each as a synthesised [`SET_MEMBER_ACCESSION`] parameter whose
    /// `ID` is the member's run identifier and whose value is that run's
    /// [`RUN_NAME_ACCESSION`] value. A fixed `<cvList>` of three vocabularies
    /// closes the document.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when the output would exceed
    /// [`MAX_OUTPUT_BYTES`](Self::MAX_OUTPUT_BYTES),
    /// [`Error::Unsupported`] when
    /// [`WriteOptions::source_encoding_declaration`] is set and the document
    /// holds non-ASCII text, and otherwise as
    /// [`QualityParameter::to_xml_string`] and [`Attachment::to_xml_string`].
    ///
    /// [`Error::MissingInformation`] additionally when a run's or set's
    /// registered name differs from its identifier and no
    /// [`RUN_NAME_ACCESSION`] (run) or [`SET_NAME_ACCESSION`] (set) parameter
    /// carries it. The format stores a name only inside such a parameter, so
    /// the source drops any other name without a word and the reloaded entry is
    /// named after its identifier. [`WriteOptions::drop_unrepresentable`],
    /// which [`WriteOptions::source`] sets, selects that behaviour.
    pub fn to_xml_string_with_options(&self, options: &WriteOptions) -> Result<String> {
        self.preflight_output(options)?;
        if !options.drop_unrepresentable {
            self.check_names_are_persisted()?;
        }
        let mut out = String::new();
        if options.source_encoding_declaration {
            out.push_str("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        } else {
            out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        }
        if let Some(sheet) = &options.stylesheet {
            check_text(&sheet.id, QualityParameter::MAX_TEXT_BYTES)?;
            if sheet.id.is_empty() || sheet.id.contains(|c: char| c.is_whitespace() || c == '"') {
                return Err(Error::InvalidValue(
                    "qcML stylesheet id must be a non-empty quote-free token".into(),
                ));
            }
            check_text(&sheet.xslt, Stylesheet::MAX_BYTES)?;
            // The body is written verbatim, so the one thing that must not be
            // in it is the document's own closing tag.
            if sheet.xslt.contains("</qcML>") {
                return Err(Error::InvalidValue(
                    "qcML stylesheet body must not close the qcML element".into(),
                ));
            }
            out.push_str("<?xml-stylesheet type=\"text/xml\" href=\"#");
            out.push_str(&sheet.id);
            out.push_str("\"?>\n");
            out.push_str(
                "<!DOCTYPE catelog [\n  <!ATTLIST xsl:stylesheet\n  id  ID  #REQUIRED>\n  ]>\n",
            );
        }
        out.push_str("<qcML xmlns=\"https://github.com/qcML/qcml\" >\n");
        for id in union_keys(&self.run_qps, &self.run_ats) {
            out.push_str("\t<runQuality ID=\"");
            push_attribute_value(&mut out, id)?;
            out.push_str("\">\n");
            for qp in self.run_quality_parameters(id) {
                out.push_str(&qp.to_xml_string_with_options(4, options)?);
            }
            for at in self.run_attachments(id) {
                out.push_str(&at.to_xml_string_with_options(4, options)?);
            }
            out.push_str("\t</runQuality>\n");
        }
        for id in union_keys(&self.set_qps, &self.set_ats) {
            out.push_str("\t<setQuality ID=\"");
            push_attribute_value(&mut out, id)?;
            out.push_str("\">\n");
            if let Some(members) = self.set_member_names.get(id) {
                for member in members {
                    // The source resolves a member only against the run
                    // identifiers; a member recorded as a run *name* - which is
                    // what its own reader records from MS:1000577 - finds no run
                    // and is silently skipped ("TODO warn - no mzML file
                    // registered for this run"). The native default also tries
                    // the name map and refuses a member it cannot resolve.
                    let resolved = if options.drop_unrepresentable {
                        self.run_qps.get(member).map(|_| member.as_str())
                    } else {
                        Self::resolve(&self.run_qps, &self.run_name_ids, member)
                    };
                    let Some(member) = resolved else {
                        if options.drop_unrepresentable {
                            continue;
                        }
                        return Err(Error::MissingInformation(format!(
                            "qcML set member {member:?} matches no registered run"
                        )));
                    };
                    let Some(run) = self.run_qps.get(member) else {
                        continue;
                    };
                    let mut qp = QualityParameter {
                        id: member.to_owned(),
                        name: "set name".into(),
                        cv_ref: "QC".into(),
                        cv_acc: SET_MEMBER_ACCESSION.into(),
                        ..QualityParameter::default()
                    };
                    for candidate in run {
                        if candidate.cv_acc == RUN_NAME_ACCESSION {
                            qp.value = candidate.value.clone();
                        }
                    }
                    out.push_str(&qp.to_xml_string_with_options(4, options)?);
                }
            }
            for qp in self.set_quality_parameters(id) {
                out.push_str(&qp.to_xml_string_with_options(4, options)?);
            }
            for at in self.set_attachments(id) {
                out.push_str(&at.to_xml_string_with_options(4, options)?);
            }
            out.push_str("\t</setQuality>\n");
        }
        out.push_str("\t<cvList>\n");
        out.push_str("\t<cv uri=\"http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" ID=\"psi_cv_ref\" fullName=\"PSI-MS\" version=\"3.41.0\"/>\n");
        out.push_str("\t<cv uri=\"https://github.com/qcML/qcML-development/blob/master/cv/qc-cv.obo\" ID=\"qc_cv_ref\" fullName=\"QC-CV\" version=\"0.1.1\"/>\n");
        out.push_str("\t<cv uri=\"http://obo.cvs.sourceforge.net/viewvc/obo/obo/ontology/phenotype/unit.obo\" ID=\"uo_cv_ref\" fullName=\"unit\" version=\"1.0.0\"/>\n");
        out.push_str("\t</cvList>\n");
        if let Some(sheet) = &options.stylesheet {
            out.push_str(&sheet.xslt);
            out.push('\n');
        }
        out.push_str("</qcML>\n");
        if options.source_encoding_declaration && !out.is_ascii() {
            return Err(Error::Unsupported(
                "qcML source encoding declaration is ISO-8859-1 and cannot label non-ASCII text"
                    .into(),
            ));
        }
        if out.len() > Self::MAX_OUTPUT_BYTES {
            return Err(limit("qcML output byte limit exceeded"));
        }
        Ok(out)
    }

    /// Write the document to `path`, under the native default options.
    ///
    /// The bytes are serialised in full and published atomically, so a failure
    /// never leaves a partial report behind. The source streams straight into an
    /// `ofstream` and throws `Exception::UnableToCreateFile` if the stream
    /// cannot be opened, leaving whatever it had already written.
    ///
    /// # Errors
    ///
    /// As [`to_xml_string_with_options`](Self::to_xml_string_with_options), plus
    /// [`Error::Io`].
    pub fn store(&self, path: impl AsRef<Path>) -> Result<()> {
        self.store_with_options(path, &WriteOptions::default())
    }

    /// Write the document to `path` under explicit options.
    ///
    /// # Errors
    ///
    /// As [`store`](Self::store).
    pub fn store_with_options(&self, path: impl AsRef<Path>, options: &WriteOptions) -> Result<()> {
        let text = self.to_xml_string_with_options(options)?;
        super::path_io::write(path.as_ref(), |writer| {
            writer.write_all(text.as_bytes())?;
            Ok(())
        })
    }

    // A run's or set's registered name is written nowhere of its own: the reader
    // recovers it from an MS:1000577 parameter inside a runQuality and from a
    // QC:0000058 parameter inside a setQuality, and falls back to the entry's ID
    // when neither is present. The source therefore loses every name that is not
    // also carried by such a parameter, silently.
    fn check_names_are_persisted(&self) -> Result<()> {
        for (names, data, accession, kind) in [
            (&self.run_name_ids, &self.run_qps, RUN_NAME_ACCESSION, "run"),
            (&self.set_name_ids, &self.set_qps, SET_NAME_ACCESSION, "set"),
        ] {
            for (name, id) in names {
                if name == id {
                    continue;
                }
                let carried = data.get(id).is_some_and(|list| {
                    list.iter()
                        .any(|qp| qp.cv_acc == accession && &qp.value == name)
                });
                if !carried {
                    return Err(Error::MissingInformation(format!(
                        "qcML {kind} {id:?} is named {name:?} but carries no {accession} \
                         parameter with that value, so the name would not survive the write"
                    )));
                }
            }
        }
        Ok(())
    }

    fn preflight_output(&self, options: &WriteOptions) -> Result<()> {
        let mut bytes = 512
            + options
                .stylesheet
                .as_ref()
                .map_or(0, |s| s.xslt.len() + s.id.len() + 256);
        for (qps, ats) in [
            (&self.run_qps, &self.run_ats),
            (&self.set_qps, &self.set_ats),
        ] {
            for (id, list) in qps {
                bytes = bytes
                    .checked_add(id.len().saturating_add(64))
                    .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
                for qp in list {
                    bytes = bytes
                        .checked_add(qp.weight())
                        .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
                }
                if bytes > Self::MAX_OUTPUT_BYTES {
                    return Err(limit("qcML output byte limit exceeded"));
                }
            }
            for (id, list) in ats {
                bytes = bytes
                    .checked_add(id.len().saturating_add(64))
                    .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
                for at in list {
                    bytes = bytes
                        .checked_add(at.weight())
                        .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
                }
                if bytes > Self::MAX_OUTPUT_BYTES {
                    return Err(limit("qcML output byte limit exceeded"));
                }
            }
        }
        for (id, members) in &self.set_member_names {
            bytes = bytes
                .checked_add(id.len().saturating_add(32))
                .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
            for member in members {
                bytes = bytes
                    .checked_add(member.len().saturating_add(160))
                    .ok_or_else(|| limit("qcML output byte limit exceeded"))?;
            }
            if bytes > Self::MAX_OUTPUT_BYTES {
                return Err(limit("qcML output byte limit exceeded"));
            }
        }
        Ok(())
    }
}

/// Render a two-level string table as delimiter-separated text.
///
/// The outer key becomes the first cell of each line, under the literal header
/// `qp`. Column names are taken from the **first** row's keys, in lexical order,
/// and every line — the header included — ends with a trailing separator. An
/// empty table renders as the empty string.
///
/// The source declares this as a `const` member function of `QcMLFile` that
/// never touches `this`, so it is a free function here.
///
/// # Warning
///
/// A row missing one of the header's columns emits neither the cell nor its
/// separator, so the remaining cells shift left and the line no longer lines up
/// with the header; a row with extra keys loses them entirely. Both are the
/// source's behaviour — its own comment on the missing case reads "TODO else
/// throw error" — and both are preserved because
/// [`QcMLFile::export_id_stats`] and the qcML report stylesheet consume exactly
/// this output. Pass rows with identical key sets to get a well-formed table.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `separator` is empty, and
/// [`Error::InvalidRange`] when the text would exceed
/// [`QcMLFile::MAX_OUTPUT_BYTES`].
pub fn map_to_csv(
    table: &BTreeMap<String, BTreeMap<String, String>>,
    separator: &str,
) -> Result<String> {
    if separator.is_empty() {
        return Err(Error::InvalidValue(
            "qcML CSV separator must not be empty".into(),
        ));
    }
    let mut out = String::new();
    let Some((_, first)) = table.iter().next() else {
        return Ok(out);
    };
    let columns: Vec<&String> = first.keys().collect();
    let mut budget = 0usize;
    for (row, cells) in table {
        budget = budget
            .checked_add(row.len().saturating_add(separator.len()))
            .ok_or_else(|| limit("qcML CSV byte limit exceeded"))?;
        for (name, value) in cells {
            budget = budget
                .checked_add(
                    name.len()
                        .saturating_add(value.len())
                        .saturating_add(2 * separator.len()),
                )
                .ok_or_else(|| limit("qcML CSV byte limit exceeded"))?;
        }
        if budget > QcMLFile::MAX_OUTPUT_BYTES {
            return Err(limit("qcML CSV byte limit exceeded"));
        }
    }
    out.push_str("qp");
    out.push_str(separator);
    for column in &columns {
        out.push_str(column);
        out.push_str(separator);
    }
    out.push('\n');
    for (row, cells) in table {
        out.push_str(row);
        out.push_str(separator);
        for column in &columns {
            if let Some(value) = cells.get(*column) {
                out.push_str(value);
                out.push_str(separator);
            }
        }
        out.push('\n');
    }
    Ok(out)
}

/// Read a qcML document from `path` under the default [`Limits`].
///
/// The source's member `load` clears the object, parses into it and leaves it
/// empty when the parse throws. This returns a new document instead, so a
/// failed read cannot damage one a caller already holds; assign the result to
/// replace a document.
///
/// Compression is detected from the leading bytes, not the file name, so a
/// gzip- or bzip2-compressed qcML reads as well as a plain one.
///
/// # Errors
///
/// As [`read_with_limits`], plus [`Error::Io`].
pub fn load(path: impl AsRef<Path>) -> Result<QcMLFile> {
    load_with_limits(path, &Limits::default())
}

/// Read a qcML document from `path` under explicit [`Limits`].
///
/// # Errors
///
/// As [`read_with_limits`], plus [`Error::Io`].
pub fn load_with_limits(path: impl AsRef<Path>, limits: &Limits) -> Result<QcMLFile> {
    read_with_limits(super::path_io::open(path.as_ref())?, limits)
}

/// Read a qcML document from `input` under the default [`Limits`].
///
/// # Errors
///
/// As [`read_with_limits`].
pub fn read(input: impl BufRead) -> Result<QcMLFile> {
    read_with_limits(input, &Limits::default())
}

/// Read a qcML document from `input` under explicit [`Limits`].
///
/// Bytes are decoded before parsing: UTF-8 with or without a byte-order mark,
/// UTF-16 of either endianness, and — because the source's writer declares
/// `ISO-8859-1` and then emits raw `std::string` bytes — Latin-1 when the
/// declaration says so and the bytes are not valid UTF-8.
///
/// A DOCTYPE is accepted only when it declares no entity and names no external
/// subset, which is exactly the shape the source's own writer emits alongside a
/// report stylesheet. Anything else is refused rather than expanded.
///
/// # Errors
///
/// [`Error::Parse`] for malformed XML, a missing required attribute, a
/// `qualityParameter`, `attachment`, `binary`, `tableColumnTypes` or
/// `tableRowValues` element outside the entry it belongs to, a nested or
/// duplicated `runQuality`/`setQuality`, or any exceeded [`Limits`] ceiling.
/// [`Error::Unsupported`] for a document that is not one of the decodable
/// encodings, and for a DTD this port will not process.
pub fn read_with_limits(mut input: impl BufRead, limits: &Limits) -> Result<QcMLFile> {
    let mut bytes: Vec<u8> = Vec::new();
    loop {
        let buffer = input.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        let take = buffer.len().min(8192);
        if bytes.len().saturating_add(take) > limits.max_input_bytes {
            return Err(parse_error(0, "qcML input byte limit exceeded"));
        }
        bytes.extend_from_slice(&buffer[..take]);
        input.consume(take);
    }
    let text = decode(&bytes, limits)?;
    Parser::new(&text, limits).run()
}

fn decode(bytes: &[u8], limits: &Limits) -> Result<String> {
    let utf16 = if bytes.starts_with(&[0xff, 0xfe]) {
        Some((true, 2))
    } else if bytes.starts_with(&[0xfe, 0xff]) {
        Some((false, 2))
    } else if bytes.starts_with(&[b'<', 0, b'?', 0]) {
        Some((true, 0))
    } else if bytes.starts_with(&[0, b'<', 0, b'?']) {
        Some((false, 0))
    } else {
        None
    };
    let text = if let Some((little, offset)) = utf16 {
        let body = bytes.get(offset..).unwrap_or(&[]);
        let expected: &str = if little { "utf-16le" } else { "utf-16be" };
        if body.len() % 2 != 0 {
            return Err(parse_error(0, "odd UTF-16 byte count in qcML"));
        }
        let units = body.chunks_exact(2).map(|pair| {
            if little {
                u16::from_le_bytes([pair[0], pair[1]])
            } else {
                u16::from_be_bytes([pair[0], pair[1]])
            }
        });
        let mut out = String::new();
        for c in char::decode_utf16(units) {
            let c = c.map_err(|_| parse_error(0, "invalid UTF-16 in qcML"))?;
            if out.len().saturating_add(c.len_utf8()) > limits.max_input_bytes {
                return Err(parse_error(0, "qcML input byte limit exceeded"));
            }
            out.push(c);
        }
        let declared = declaration(out.as_bytes())?.and_then(|(_, encoding)| encoding);
        if declared
            .as_deref()
            .is_some_and(|v| !matches!(v, "utf-16" | "utf16") && v != expected)
        {
            return Err(Error::Unsupported(
                "qcML XML declaration conflicts with its UTF-16 bytes".into(),
            ));
        }
        out
    } else {
        let body = bytes.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(bytes);
        let declared = declaration(body)?.and_then(|(_, encoding)| encoding);
        let latin1 = match declared.as_deref() {
            None | Some("utf-8" | "utf8") => false,
            Some("us-ascii" | "ascii") => {
                if !body.is_ascii() {
                    return Err(Error::Unsupported(
                        "qcML declares US-ASCII but holds non-ASCII bytes".into(),
                    ));
                }
                false
            }
            Some("iso-8859-1" | "iso8859-1" | "latin1" | "latin-1") => true,
            Some(other) => {
                return Err(Error::Unsupported(format!("qcML XML encoding {other}")));
            }
        };
        match std::str::from_utf8(body) {
            // A document the source wrote from UTF-8 std::strings is valid UTF-8
            // under its ISO-8859-1 declaration, so a decodable UTF-8 body always
            // wins over the byte-wise Latin-1 reading.
            Ok(text) => text.to_owned(),
            Err(_) if latin1 => {
                // The source's writer declares ISO-8859-1 and writes bytes
                // unchanged, so its own non-UTF-8 output is Latin-1.
                if body.len().saturating_mul(2) > limits.max_input_bytes {
                    return Err(parse_error(0, "qcML input byte limit exceeded"));
                }
                body.iter().map(|b| char::from(*b)).collect()
            }
            Err(_) => {
                return Err(Error::Unsupported(
                    "qcML input requires UTF-8, UTF-16 or a declared ISO-8859-1 encoding".into(),
                ));
            }
        }
    };
    if text.len() > limits.max_input_bytes {
        return Err(parse_error(0, "qcML input byte limit exceeded"));
    }
    if !text.chars().all(xml_char) {
        return Err(parse_error(0, "invalid XML 1.0 character in qcML"));
    }
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// The XML declaration's version and lower-cased encoding, when the bytes open
/// with one. `<?xml-stylesheet ...?>` is not a declaration and yields `None`.
type DeclaredText = (Option<String>, Option<String>);
fn declaration(bytes: &[u8]) -> Result<Option<DeclaredText>> {
    if !bytes.starts_with(b"<?xml") {
        return Ok(None);
    }
    if !bytes.get(5).is_some_and(u8::is_ascii_whitespace) {
        return Ok(None);
    }
    let head = bytes.get(..bytes.len().min(512)).unwrap_or(bytes);
    let end = head
        .windows(2)
        .position(|w| w == b"?>")
        .ok_or_else(|| parse_error(1, "unterminated qcML XML declaration"))?;
    let body = head.get(5..end).unwrap_or(&[]);
    if !body.is_ascii() {
        return Err(parse_error(1, "non-ASCII qcML XML declaration"));
    }
    let body = std::str::from_utf8(body)
        .map_err(|_| parse_error(1, "invalid qcML XML declaration"))?
        .replace('\'', "\"");
    let mut version = None;
    let mut encoding = None;
    let mut fields = body.split('"');
    while let Some(key) = fields.next() {
        let Some(value) = fields.next() else { break };
        let key = key.trim().trim_end_matches('=').trim();
        match key {
            "version" => version = Some(value.to_owned()),
            "encoding" => encoding = Some(value.to_ascii_lowercase()),
            "standalone" | "" => {}
            other => {
                return Err(parse_error(
                    1,
                    format!("unexpected qcML XML declaration attribute {other:?}"),
                ));
            }
        }
    }
    if version.as_deref().is_some_and(|v| v != "1.0") {
        return Err(Error::Unsupported(
            "qcML XML declaration needs version 1.0".into(),
        ));
    }
    Ok(Some((version, encoding)))
}

fn union_keys<'a, A, B>(
    first: &'a BTreeMap<String, A>,
    second: &'a BTreeMap<String, B>,
) -> Vec<&'a str> {
    let mut keys: BTreeSet<&str> = first.keys().map(String::as_str).collect();
    keys.extend(second.keys().map(String::as_str));
    keys.into_iter().collect()
}

fn collapse<T: Ord + Clone>(list: &mut Vec<T>, options: &MergeOptions, same: fn(&T, &T) -> bool) {
    list.sort();
    if options.collapse_by_name {
        list.dedup_by(|a, b| same(a, b));
    } else {
        list.dedup();
    }
}

fn preflight_merge(target: &QcMLFile, addendum: &QcMLFile) -> Result<()> {
    let mut bytes = 0usize;
    for file in [target, addendum] {
        for (qps, ats) in [
            (&file.run_qps, &file.run_ats),
            (&file.set_qps, &file.set_ats),
        ] {
            for (id, list) in qps {
                bytes = bytes
                    .checked_add(id.len().saturating_add(64))
                    .ok_or_else(|| limit("qcML merge byte limit exceeded"))?;
                for qp in list {
                    bytes = bytes
                        .checked_add(qp.weight())
                        .ok_or_else(|| limit("qcML merge byte limit exceeded"))?;
                }
                if bytes > QcMLFile::MAX_OUTPUT_BYTES {
                    return Err(limit("qcML merge byte limit exceeded"));
                }
            }
            for (id, list) in ats {
                bytes = bytes
                    .checked_add(id.len().saturating_add(64))
                    .ok_or_else(|| limit("qcML merge byte limit exceeded"))?;
                for at in list {
                    bytes = bytes
                        .checked_add(at.weight())
                        .ok_or_else(|| limit("qcML merge byte limit exceeded"))?;
                }
                if bytes > QcMLFile::MAX_OUTPUT_BYTES {
                    return Err(limit("qcML merge byte limit exceeded"));
                }
            }
        }
    }
    Ok(())
}

fn indent(level: u32) -> Result<String> {
    if level > QualityParameter::MAX_INDENTATION {
        return Err(limit("qcML indentation level limit exceeded"));
    }
    Ok("\t".repeat(level as usize))
}

fn check_parameter(qp: &QualityParameter) -> Result<()> {
    for field in qp.fields() {
        check_text(field, QualityParameter::MAX_TEXT_BYTES)?;
    }
    Ok(())
}

fn check_attachment(at: &Attachment) -> Result<()> {
    for field in at.scalar_fields() {
        check_text(field, Attachment::MAX_TEXT_BYTES)?;
    }
    check_text(&at.binary, Attachment::MAX_TEXT_BYTES)?;
    at.preflight_table()
}

fn check_text(value: &str, maximum: usize) -> Result<()> {
    if value.len() > maximum {
        return Err(limit("qcML field byte limit exceeded"));
    }
    if !value.chars().all(xml_char) {
        return Err(Error::InvalidValue(
            "qcML text holds a character XML 1.0 cannot represent".into(),
        ));
    }
    Ok(())
}

fn trimmed(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r')
}

fn join_substituted(out: &mut String, cells: &[String], separator: &str, replacement: &str) {
    for (index, cell) in cells.iter().enumerate() {
        if index > 0 {
            out.push_str(separator);
        }
        out.push_str(&cell.replace(separator, replacement));
    }
}

fn attribute(out: &mut String, name: &str, value: &str) -> Result<()> {
    out.push(' ');
    out.push_str(name);
    out.push_str("=\"");
    push_attribute_value(out, value)?;
    out.push('"');
    Ok(())
}

fn push_attribute_value(out: &mut String, value: &str) -> Result<()> {
    for c in value.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            // Literal whitespace in an attribute is normalized to a space by
            // any conforming reader, so it is written as a character reference.
            '\t' => out.push_str("&#9;"),
            '\n' => out.push_str("&#10;"),
            '\r' => out.push_str("&#13;"),
            _ if xml_char(c) => out.push(c),
            _ => {
                return Err(Error::InvalidValue(
                    "qcML attribute holds a character XML 1.0 cannot represent".into(),
                ));
            }
        }
    }
    Ok(())
}

fn push_text(out: &mut String, value: &str) -> Result<()> {
    for c in value.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            _ if xml_char(c) => out.push(c),
            _ => {
                return Err(Error::InvalidValue(
                    "qcML text holds a character XML 1.0 cannot represent".into(),
                ));
            }
        }
    }
    Ok(())
}

fn push_cells(
    out: &mut String,
    cells: &[String],
    options: &WriteOptions,
    label: &str,
    substitute_spaces: bool,
) -> Result<()> {
    let mut line = String::new();
    for (index, cell) in cells.iter().enumerate() {
        if index > 0 {
            line.push(' ');
        }
        if options.source_table_text {
            // Column types are substituted, row values are not: the source
            // builds a substituted copy of each row and then concatenates the
            // original.
            if substitute_spaces {
                line.push_str(&cell.replace(' ', "_"));
            } else {
                line.push_str(cell);
            }
            continue;
        }
        if cell.is_empty() {
            return Err(Error::InvalidValue(format!(
                "qcML {label} must not be empty; the space-delimited table cannot represent it"
            )));
        }
        if cell.contains(trimmed) {
            return Err(Error::InvalidValue(format!(
                "qcML {label} must not contain whitespace; the table is space-delimited"
            )));
        }
        line.push_str(cell);
    }
    push_text(out, line.trim_matches(trimmed))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum EntryKind {
    Run,
    Set,
}

struct Entry {
    kind: EntryKind,
    id: String,
    name: String,
    members: BTreeSet<String>,
    qps: Vec<QualityParameter>,
    ats: Vec<Attachment>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Capture {
    Binary,
    ColumnTypes,
    RowValues,
}

struct Parser<'a> {
    text: &'a str,
    limits: Limits,
    out: QcMLFile,
    depth: usize,
    elements: usize,
    /// Bytes of `text` the line counter has already walked.
    scanned: usize,
    /// One-based line number at byte `scanned`.
    scanned_line: usize,
    /// Cells of the table being read, counted as its rows arrive.
    table_cells: usize,
    root: bool,
    entry: Option<Entry>,
    open_qp: Option<QualityParameter>,
    open_at: Option<Attachment>,
    capture: Option<Capture>,
    buffer: String,
}

impl<'a> Parser<'a> {
    fn new(text: &'a str, limits: &Limits) -> Self {
        Self {
            text,
            limits: *limits,
            out: QcMLFile::default(),
            depth: 0,
            elements: 0,
            scanned: 0,
            scanned_line: 1,
            table_cells: 0,
            root: false,
            entry: None,
            open_qp: None,
            open_at: None,
            capture: None,
            buffer: String::new(),
        }
    }

    /// The one-based line number at byte `offset`, carried forward.
    ///
    /// The parser's byte position only ever advances, so the newlines between
    /// the previous event and this one are counted once and added to the
    /// running total: every byte of the document is examined exactly once over
    /// the whole parse. Counting over the prefix on each event instead cost
    /// O(events x document bytes), which made a document of many small
    /// elements - or many comments - quadratic in its own size. A position
    /// that did not advance (`get` yields `None`) reuses the last line.
    fn line_at(&mut self, offset: usize) -> usize {
        let offset = offset.min(self.text.len());
        if let Some(span) = self.text.as_bytes().get(self.scanned..offset) {
            self.scanned_line += span.iter().filter(|b| **b == b'\n').count();
            self.scanned = offset;
        }
        self.scanned_line
    }

    /// Charge one markup node against [`Limits::max_elements`].
    fn count_node(&mut self, line: usize) -> Result<()> {
        self.elements += 1;
        if self.elements > self.limits.max_elements {
            return Err(parse_error(line, "qcML element limit exceeded"));
        }
        Ok(())
    }

    fn run(mut self) -> Result<QcMLFile> {
        let mut reader = Reader::from_str(self.text);
        reader.config_mut().check_end_names = true;
        reader.config_mut().check_comments = true;
        let mut stack: Vec<String> = Vec::new();
        let mut doctype = false;
        loop {
            let offset = reader.buffer_position() as usize;
            let line = self.line_at(offset);
            let event = reader
                .read_event()
                .map_err(|e| parse_error(line, e.to_string()))?;
            let self_closing = matches!(&event, Event::Empty(_));
            match event {
                Event::Start(e) | Event::Empty(e) => {
                    let tag = tag_name(&e, line)?.to_owned();
                    self.count_node(line)?;
                    if self.depth + 1 > self.limits.max_depth {
                        return Err(parse_error(line, "qcML nesting limit exceeded"));
                    }
                    if stack.is_empty() && self.root {
                        return Err(parse_error(line, "more than one qcML root"));
                    }
                    let attributes = self.attributes(&e, line)?;
                    let parent = stack.last().map(String::as_str).unwrap_or("");
                    self.start(&tag, parent, &attributes, line)?;
                    if self_closing {
                        // A self-closing element opens and closes at once; the
                        // source's SAX driver also calls both callbacks.
                        let parent = stack.last().map(String::as_str).unwrap_or("");
                        self.end(&tag, parent, line)?;
                    } else {
                        stack.push(tag);
                        self.depth += 1;
                    }
                }
                Event::End(e) => {
                    let tag = std::str::from_utf8(e.name().0)
                        .map_err(|_| parse_error(line, "invalid qcML tag name"))?
                        .to_owned();
                    let popped = stack.pop();
                    if popped.as_deref() != Some(tag.as_str()) {
                        return Err(parse_error(line, "mismatched qcML closing tag"));
                    }
                    self.depth = self.depth.saturating_sub(1);
                    let parent = stack.last().map(String::as_str).unwrap_or("");
                    self.end(&tag, parent, line)?;
                }
                Event::Text(e) => {
                    let raw = std::str::from_utf8(e.as_ref())
                        .map_err(|_| parse_error(line, "invalid UTF-8 in qcML text"))?;
                    if stack.is_empty() && !raw.trim_matches(trimmed).is_empty() {
                        return Err(parse_error(line, "text outside the qcML root"));
                    }
                    if self.capture.is_some() {
                        let decoded = quick_xml::escape::unescape(raw)
                            .map_err(|e| parse_error(line, e.to_string()))?;
                        self.push_captured(&decoded, line)?;
                    }
                }
                Event::CData(e) => {
                    if stack.is_empty() {
                        return Err(parse_error(line, "CDATA outside the qcML root"));
                    }
                    if self.capture.is_some() {
                        let raw = std::str::from_utf8(e.as_ref())
                            .map_err(|_| parse_error(line, "invalid UTF-8 in qcML CDATA"))?;
                        let owned = raw.to_owned();
                        self.push_captured(&owned, line)?;
                    }
                }
                Event::GeneralRef(e) => {
                    if stack.is_empty() {
                        return Err(parse_error(line, "entity outside the qcML root"));
                    }
                    let raw = std::str::from_utf8(e.as_ref())
                        .map_err(|_| parse_error(line, "invalid qcML entity"))?;
                    let encoded = format!("&{raw};");
                    let decoded = quick_xml::escape::unescape(&encoded)
                        .map_err(|e| parse_error(line, e.to_string()))?;
                    if self.capture.is_some() {
                        let owned = decoded.into_owned();
                        self.push_captured(&owned, line)?;
                    }
                }
                Event::Decl(_) => {}
                Event::PI(e) => {
                    self.count_node(line)?;
                    let target = std::str::from_utf8(e.target())
                        .map_err(|_| parse_error(line, "invalid qcML processing instruction"))?;
                    if target.eq_ignore_ascii_case("xml") {
                        return Err(parse_error(line, "reserved qcML processing instruction"));
                    }
                }
                Event::DocType(e) => {
                    if doctype {
                        return Err(parse_error(line, "more than one qcML DOCTYPE"));
                    }
                    doctype = true;
                    if e.len() > self.limits.max_doctype_bytes {
                        return Err(parse_error(line, "qcML DOCTYPE byte limit exceeded"));
                    }
                    let raw = std::str::from_utf8(e.as_ref())
                        .map_err(|_| parse_error(line, "invalid qcML DOCTYPE"))?;
                    let upper = raw.to_ascii_uppercase();
                    if upper.contains("ENTITY")
                        || upper.contains("SYSTEM")
                        || upper.contains("PUBLIC")
                    {
                        return Err(Error::Unsupported(
                            "DTD entities or an external subset in qcML".into(),
                        ));
                    }
                }
                // A comment holds no data but costs the parser the same walk as
                // an element, so it draws on the same ceiling; charged to no
                // ceiling at all, a document of nothing but comments was
                // bounded only by `max_input_bytes`.
                Event::Comment(_) => self.count_node(line)?,
                Event::Eof => {
                    if !stack.is_empty() {
                        return Err(parse_error(line, "incomplete qcML document"));
                    }
                    if !self.root {
                        return Err(parse_error(line, "no qcML root element"));
                    }
                    break;
                }
            }
        }
        Ok(self.out)
    }

    fn attributes(&self, element: &BytesStart<'_>, line: usize) -> Result<Vec<(String, String)>> {
        let mut out: Vec<(String, String)> = Vec::new();
        for attribute in element.attributes().with_checks(false) {
            let attribute = attribute.map_err(|e| parse_error(line, e.to_string()))?;
            // The duplicate scan below is quadratic in the attribute count, and
            // no qcML element declares more than nine attributes.
            if out.len() >= MAX_ATTRIBUTES {
                return Err(parse_error(line, "qcML attribute count limit exceeded"));
            }
            let key = std::str::from_utf8(attribute.key.as_ref())
                .map_err(|_| parse_error(line, "invalid qcML attribute name"))?;
            if out.iter().any(|(previous, _)| previous == key) {
                return Err(parse_error(line, "duplicate qcML attribute"));
            }
            if attribute.value.contains(&b'<') {
                return Err(parse_error(line, "raw '<' in a qcML attribute"));
            }
            let raw = std::str::from_utf8(&attribute.value)
                .map_err(|_| parse_error(line, "invalid UTF-8 in a qcML attribute"))?;
            if raw.len() > self.limits.max_text_bytes {
                return Err(parse_error(line, "qcML attribute byte limit exceeded"));
            }
            let mut normalized = String::with_capacity(raw.len());
            for c in raw.chars() {
                normalized.push(if matches!(c, '\t' | '\n' | '\r') {
                    ' '
                } else {
                    c
                });
            }
            let value = quick_xml::escape::unescape(&normalized)
                .map_err(|e| parse_error(line, e.to_string()))?;
            out.push((key.to_owned(), value.into_owned()));
        }
        Ok(out)
    }

    fn push_captured(&mut self, text: &str, line: usize) -> Result<()> {
        if self.buffer.len().saturating_add(text.len()) > self.limits.max_text_bytes {
            return Err(parse_error(line, "qcML character data limit exceeded"));
        }
        self.buffer.push_str(text);
        Ok(())
    }

    fn start(
        &mut self,
        tag: &str,
        parent: &str,
        attributes: &[(String, String)],
        line: usize,
    ) -> Result<()> {
        match tag {
            "qcML" => {
                if self.root {
                    return Err(parse_error(line, "more than one qcML root"));
                }
                if !parent.is_empty() {
                    return Err(parse_error(line, "nested qcML element"));
                }
                self.root = true;
            }
            "runQuality" | "setQuality" => {
                if self.entry.is_some() {
                    return Err(parse_error(line, "nested qcML run or set"));
                }
                let kind = if tag == "runQuality" {
                    EntryKind::Run
                } else {
                    EntryKind::Set
                };
                self.entry = Some(Entry {
                    kind,
                    id: required(attributes, "ID", tag, line)?,
                    name: String::new(),
                    members: BTreeSet::new(),
                    qps: Vec::new(),
                    ats: Vec::new(),
                });
            }
            "qualityParameter" => {
                if self.entry.is_none() {
                    return Err(parse_error(
                        line,
                        "qualityParameter outside a runQuality or setQuality",
                    ));
                }
                if self.open_qp.is_some() || self.open_at.is_some() {
                    return Err(parse_error(line, "nested qcML qualityParameter"));
                }
                self.open_qp = Some(QualityParameter {
                    name: required(attributes, "name", tag, line)?,
                    id: required(attributes, "ID", tag, line)?,
                    value: optional(attributes, "value"),
                    cv_ref: required(attributes, "cvRef", tag, line)?,
                    cv_acc: required(attributes, "accession", tag, line)?,
                    unit_ref: unit(attributes, "unitCvRef", "unitRef"),
                    unit_acc: unit(attributes, "unitAccession", "unitAcc"),
                    flag: optional(attributes, "flag"),
                });
            }
            "attachment" => {
                if self.entry.is_none() {
                    return Err(parse_error(
                        line,
                        "attachment outside a runQuality or setQuality",
                    ));
                }
                if self.open_qp.is_some() || self.open_at.is_some() {
                    return Err(parse_error(line, "nested qcML attachment"));
                }
                self.open_at = Some(Attachment {
                    name: required(attributes, "name", tag, line)?,
                    id: required(attributes, "ID", tag, line)?,
                    value: optional(attributes, "value"),
                    cv_ref: required(attributes, "cvRef", tag, line)?,
                    cv_acc: required(attributes, "accession", tag, line)?,
                    unit_ref: unit(attributes, "unitCvRef", "unitRef"),
                    unit_acc: unit(attributes, "unitAccession", "unitAcc"),
                    binary: String::new(),
                    // Required by the source's reader but omitted by its own
                    // writer when empty, so it is optional here.
                    quality_ref: optional(attributes, "qualityParameterRef"),
                    col_types: Vec::new(),
                    table_rows: Vec::new(),
                });
                self.table_cells = 0;
            }
            "binary" | "tableColumnTypes" | "tableRowValues" => {
                if self.open_at.is_none() {
                    return Err(parse_error(
                        line,
                        "table or binary content outside an attachment",
                    ));
                }
                if self.capture.is_some() {
                    return Err(parse_error(line, "nested qcML character content"));
                }
                self.capture = Some(match tag {
                    "binary" => Capture::Binary,
                    "tableColumnTypes" => Capture::ColumnTypes,
                    _ => Capture::RowValues,
                });
                self.buffer.clear();
            }
            _ => {}
        }
        Ok(())
    }

    fn end(&mut self, tag: &str, _parent: &str, line: usize) -> Result<()> {
        match tag {
            "binary" | "tableColumnTypes" | "tableRowValues" => {
                let capture = self
                    .capture
                    .take()
                    .ok_or_else(|| parse_error(line, "unbalanced qcML character content"))?;
                let text = std::mem::take(&mut self.buffer);
                let at = self
                    .open_at
                    .as_mut()
                    .ok_or_else(|| parse_error(line, "character content outside an attachment"))?;
                // Every ceiling below is decided on the captured text, which is
                // already bounded by `max_text_bytes`, and before the `Vec` and
                // `String`s of the cells exist. Reading them off the assembled
                // table instead - as `Attachment::preflight_table` does when
                // the attachment is finally committed - lets a document build
                // the whole oversized table in memory first.
                match capture {
                    // The source concatenates every notification, so repeated
                    // <binary> elements in one attachment accumulate; the
                    // per-payload ceiling therefore has to be charged here and
                    // not once at the end.
                    Capture::Binary => {
                        if at.binary.len().saturating_add(text.len()) > Attachment::MAX_TEXT_BYTES {
                            return Err(parse_error(line, "qcML attachment byte limit exceeded"));
                        }
                        at.binary.push_str(&text);
                    }
                    Capture::ColumnTypes => {
                        let columns = count_cells(&text);
                        if columns > Attachment::MAX_COLUMNS {
                            return Err(parse_error(line, "qcML attachment column limit exceeded"));
                        }
                        // The column types are part of the cell count, as in
                        // `preflight_table`. They are added rather than
                        // assigned so that rows read before them - and a
                        // repeated <tableColumnTypes>, which overwrites the
                        // list the source also overwrites - still count.
                        let total = self
                            .table_cells
                            .checked_add(columns)
                            .ok_or_else(|| limit("qcML attachment cell limit exceeded"))?;
                        if total > Attachment::MAX_TABLE_CELLS {
                            return Err(parse_error(line, "qcML attachment cell limit exceeded"));
                        }
                        at.col_types = split_cells(&text);
                        self.table_cells = total;
                    }
                    Capture::RowValues => {
                        let cells = count_cells(&text);
                        if cells > 0 {
                            if at.table_rows.len() >= Attachment::MAX_ROWS {
                                return Err(parse_error(
                                    line,
                                    "qcML attachment row limit exceeded",
                                ));
                            }
                            let total = self
                                .table_cells
                                .checked_add(cells)
                                .ok_or_else(|| limit("qcML attachment cell limit exceeded"))?;
                            if total > Attachment::MAX_TABLE_CELLS {
                                return Err(parse_error(
                                    line,
                                    "qcML attachment cell limit exceeded",
                                ));
                            }
                            at.table_rows.push(split_cells(&text));
                            self.table_cells = total;
                        }
                    }
                }
            }
            "qualityParameter" => {
                let qp = self
                    .open_qp
                    .take()
                    .ok_or_else(|| parse_error(line, "unbalanced qcML qualityParameter"))?;
                let entry = self
                    .entry
                    .as_mut()
                    .ok_or_else(|| parse_error(line, "qualityParameter outside an entry"))?;
                match entry.kind {
                    EntryKind::Run => {
                        if qp.cv_acc == RUN_NAME_ACCESSION {
                            entry.name = qp.value.clone();
                        }
                        if entry.qps.len() >= QcMLFile::MAX_PARAMETERS_PER_ENTRY {
                            return Err(parse_error(line, "qcML quality parameter limit exceeded"));
                        }
                        entry.qps.push(qp);
                    }
                    EntryKind::Set => {
                        if qp.cv_acc == RUN_NAME_ACCESSION {
                            if entry.members.len() >= QcMLFile::MAX_SET_MEMBERS {
                                return Err(parse_error(line, "qcML set member limit exceeded"));
                            }
                            entry.members.insert(qp.value);
                            return Ok(());
                        }
                        if qp.cv_acc == SET_MEMBER_ACCESSION {
                            // Native recovery: `store` documents each member as
                            // a QC:0000005 parameter whose ID is the member's
                            // run identifier, but the source's reader does not
                            // recognise it, so C++ loses set membership on every
                            // store/load round trip. The parameter is still kept
                            // as an ordinary set parameter, as in the source.
                            if entry.members.len() >= QcMLFile::MAX_SET_MEMBERS {
                                return Err(parse_error(line, "qcML set member limit exceeded"));
                            }
                            entry.members.insert(qp.id.clone());
                        }
                        if qp.cv_acc == SET_NAME_ACCESSION {
                            entry.name = qp.value.clone();
                        }
                        if entry.qps.len() >= QcMLFile::MAX_PARAMETERS_PER_ENTRY {
                            return Err(parse_error(line, "qcML quality parameter limit exceeded"));
                        }
                        entry.qps.push(qp);
                    }
                }
            }
            "attachment" => {
                let at = self
                    .open_at
                    .take()
                    .ok_or_else(|| parse_error(line, "unbalanced qcML attachment"))?;
                let entry = self
                    .entry
                    .as_mut()
                    .ok_or_else(|| parse_error(line, "attachment outside an entry"))?;
                if entry.ats.len() >= QcMLFile::MAX_ATTACHMENTS_PER_ENTRY {
                    return Err(parse_error(line, "qcML attachment limit exceeded"));
                }
                entry.ats.push(at);
            }
            "runQuality" | "setQuality" => {
                let entry = self
                    .entry
                    .take()
                    .ok_or_else(|| parse_error(line, "unbalanced qcML run or set"))?;
                if self.open_qp.is_some() || self.open_at.is_some() {
                    return Err(parse_error(line, "unclosed qcML child element"));
                }
                self.commit(entry, line)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn commit(&mut self, entry: Entry, line: usize) -> Result<()> {
        // The source defaults a nameless run or set to its own identifier; its
        // comment reads "TODO give warning that a run should have a name cv!!!".
        let name = if entry.name.is_empty() {
            entry.id.clone()
        } else {
            entry.name.clone()
        };
        let known = match entry.kind {
            EntryKind::Run => self.out.run_qps.contains_key(&entry.id),
            EntryKind::Set => self.out.set_qps.contains_key(&entry.id),
        };
        if known {
            // registerRun/registerSet reset the entry, so a repeated ID would
            // silently discard everything read for the first occurrence.
            return Err(parse_error(
                line,
                "duplicate qcML runQuality or setQuality ID",
            ));
        }
        let map_err = |e: Error| match e {
            Error::InvalidValue(message) => parse_error(line, message),
            Error::MissingInformation(message) => parse_error(line, message),
            other => other,
        };
        // Validate children once, then move their complete lists. Calling the
        // public single-child setters here rescanned and copied the enclosing
        // ID for every child: a long ID multiplied by many tiny children.
        for qp in &entry.qps {
            check_parameter(qp).map_err(map_err)?;
        }
        for at in &entry.ats {
            check_attachment(at).map_err(map_err)?;
        }
        match entry.kind {
            EntryKind::Run => {
                self.out.register_run(&entry.id, &name).map_err(map_err)?;
                self.out.run_qps.insert(entry.id.clone(), entry.qps);
                self.out.run_ats.insert(entry.id, entry.ats);
            }
            EntryKind::Set => {
                self.out
                    .register_set(&entry.id, &name, &entry.members)
                    .map_err(map_err)?;
                self.out.set_qps.insert(entry.id.clone(), entry.qps);
                self.out.set_ats.insert(entry.id, entry.ats);
            }
        }
        Ok(())
    }
}

fn tag_name<'a>(element: &'a BytesStart<'_>, line: usize) -> Result<&'a str> {
    std::str::from_utf8(element.name().0).map_err(|_| parse_error(line, "invalid qcML tag name"))
}

fn required(attributes: &[(String, String)], name: &str, tag: &str, line: usize) -> Result<String> {
    attributes
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value.clone())
        .ok_or_else(|| parse_error(line, format!("{tag} requires the attribute {name}")))
}

fn optional(attributes: &[(String, String)], name: &str) -> String {
    attributes
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value.clone())
        .unwrap_or_default()
}

fn unit(attributes: &[(String, String)], schema: &str, written: &str) -> String {
    let value = optional(attributes, schema);
    if value.is_empty() {
        return optional(attributes, written);
    }
    value
}

fn split_cells(text: &str) -> Vec<String> {
    let trimmed_text = text.trim_matches(trimmed);
    if trimmed_text.is_empty() {
        return Vec::new();
    }
    trimmed_text.split(' ').map(str::to_owned).collect()
}

/// The number of cells [`split_cells`] would return, without allocating one.
///
/// This is what lets the table ceilings be decided before the cells exist.
fn count_cells(text: &str) -> usize {
    let trimmed_text = text.trim_matches(trimmed);
    if trimmed_text.is_empty() {
        return 0;
    }
    trimmed_text.split(' ').count()
}

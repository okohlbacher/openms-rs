// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PEFF reading: `readStart`'s header parse, `readNext`/`readEntry_` and
//! `parseAnnotations_` with the tuple parsers (`PEFFFile.cpp:1360-1685`,
//! `1751-1808`, `1887-2958`).

use super::text::{
    extract_annotation_id, is_unsigned_integer, parse_custom_key_def, peff_unescape, split,
    split_by_pipe_escape_aware, to_int32, to_size, trim, validate_custom_key_field_type,
};
use super::{
    PEFFDatabaseMetadata, PEFFDisulfideBond, PEFFEntry, PEFFModification, PEFFModificationType,
    PEFFProcessedRegion, PEFFSequenceType, PEFFVariantComplex, PEFFVariantSimple, ReadLimits,
};
use crate::concept::boost_regex::BoostRegex;
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::io::{BufRead, Read};

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("PEFF {what} limit exceeded"))
}

fn utf8(bytes: Vec<u8>, what: &str) -> Result<String> {
    String::from_utf8(bytes).map_err(|_| parse_error(format!("PEFF {what} is not valid UTF-8")))
}

/// A streaming PEFF reader over any buffered input.
///
/// Construction reads the header (source `readStart`); each
/// [`next_entry`](Self::next_entry) reads one entry (source `readNext`).
/// Warnings the source writes to `OPENMS_LOG_WARN` accumulate in
/// [`warnings`](Self::warnings), up to [`ReadLimits::max_warnings`].
pub struct PEFFReader<R> {
    input: R,
    limits: ReadLimits,
    headers: Vec<PEFFDatabaseMetadata>,
    entries_read: usize,
    warnings: Vec<String>,
    suppressed: usize,
    /// The source stream's `eofbit`, set when an entry read reached the end.
    eof: bool,
}

impl<R: BufRead> PEFFReader<R> {
    /// Read the header of `input` with the default [`ReadLimits`].
    ///
    /// # Errors
    ///
    /// As [`with_limits`](Self::with_limits).
    pub fn new(input: R) -> Result<Self> {
        Self::with_limits(input, ReadLimits::default())
    }

    /// Read the header of `input`.
    ///
    /// Source `readStart` after it has opened the file: leading white space
    /// is skipped, every `#` line is a header line (carriage returns
    /// dropped), and the first other character ends the header. A block
    /// without `DbName` that closes before any database block is the file
    /// description block: its version and general comments are carried into
    /// every later database block (the comments only when the block has none
    /// of its own). A block still open when the header ends counts, unless
    /// the input ends inside it. Without any block there is one default
    /// header. Every block missing a key the specification requires is a
    /// warning, as is `ProteoformDb` together with
    /// `HasAnnotationIdentifiers`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when `NumberOfEntries` does not convert (the source's
    /// `ConversionError`) or a line is not UTF-8; [`Error::InvalidValue`]
    /// when a line or the header exceeds its bound; [`Error::Io`] when the
    /// input fails.
    pub fn with_limits(input: R, limits: ReadLimits) -> Result<Self> {
        let mut reader = Self {
            input,
            limits,
            headers: Vec::new(),
            entries_read: 0,
            warnings: Vec::new(),
            suppressed: 0,
            eof: false,
        };
        reader.read_header()?;
        Ok(reader)
    }

    /// The database headers, one per block.
    pub fn headers(&self) -> &[PEFFDatabaseMetadata] {
        &self.headers
    }

    /// The warnings so far, in order.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Remove and return the warnings so far.
    pub fn take_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.warnings)
    }

    /// Warnings beyond [`ReadLimits::max_warnings`] that were counted, not
    /// kept.
    pub fn suppressed_warnings(&self) -> usize {
        self.suppressed
    }

    /// Entries read so far.
    pub fn entries_read(&self) -> usize {
        self.entries_read
    }

    /// The headers and warnings, consuming the reader.
    pub fn into_parts(self) -> (Vec<PEFFDatabaseMetadata>, Vec<String>) {
        (self.headers, self.warnings)
    }

    fn warn(&mut self, message: String) {
        if self.warnings.len() < self.limits.max_warnings {
            self.warnings.push(message);
        } else {
            self.suppressed += 1;
        }
    }

    fn peek(&mut self) -> Result<Option<u8>> {
        Ok(self.input.fill_buf()?.first().copied())
    }

    fn bump(&mut self) -> Result<Option<u8>> {
        let byte = self.peek()?;
        if byte.is_some() {
            self.input.consume(1);
        }
        Ok(byte)
    }

    /// Whether the input has no bytes left. Source `atEnd`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the input fails.
    pub fn at_end(&mut self) -> Result<bool> {
        Ok(self.peek()?.is_none())
    }

    fn read_header(&mut self) -> Result<()> {
        let mut current = PEFFDatabaseMetadata::default();
        let mut has_header = false;
        let mut file_version = String::new();
        let mut file_comments: Vec<String> = Vec::new();
        let mut lines = 0usize;
        let merge = |header: &mut PEFFDatabaseMetadata, version: &str, comments: &[String]| {
            if !version.is_empty() {
                header.version = version.to_owned();
            }
            if !comments.is_empty() && header.general_comments.is_empty() {
                header.general_comments = comments.to_vec();
            }
        };
        while let Some(c) = self.peek()? {
            if matches!(c, b' ' | b'\t' | b'\n' | b'\r') {
                self.bump()?;
                continue;
            }
            if c != b'#' {
                if has_header {
                    merge(&mut current, &file_version, &file_comments);
                    self.headers.push(std::mem::take(&mut current));
                }
                break;
            }
            lines += 1;
            if lines > self.limits.max_header_lines {
                return Err(limit("header line"));
            }
            let mut line = Vec::new();
            loop {
                match self.bump()? {
                    None | Some(b'\n') => break,
                    Some(b'\r') => {}
                    Some(byte) => {
                        if line.len() == self.limits.max_line_bytes {
                            return Err(limit("line byte"));
                        }
                        line.push(byte);
                    }
                }
            }
            let line = utf8(line, "header line")?;
            if self.parse_header_line(&line, &mut current)? {
                if has_header {
                    if current.db_name.is_empty() && self.headers.is_empty() {
                        // The file description block.
                        file_version = current.version.clone();
                        file_comments = current.general_comments.clone();
                        current = PEFFDatabaseMetadata::default();
                    } else {
                        merge(&mut current, &file_version, &file_comments);
                        self.headers.push(std::mem::take(&mut current));
                    }
                    has_header = false;
                }
            } else {
                has_header = true;
            }
        }
        if self.headers.is_empty() {
            self.headers.push(PEFFDatabaseMetadata::default());
        }
        let mut warnings = Vec::new();
        for (index, header) in self.headers.iter().enumerate() {
            let label = if header.db_name.is_empty() {
                format!("database block {}", index + 1)
            } else {
                header.db_name.clone()
            };
            if header.db_name.is_empty() {
                warnings.push(format!(
                    "PEFF: Database block should start with 'DbName' key in {label}."
                ));
            }
            if header.prefix.is_empty() {
                warnings.push(format!("PEFF: Required key 'Prefix' missing in {label}."));
            }
            if header.db_version.is_empty() {
                warnings.push(format!(
                    "PEFF: Required key 'DbVersion' missing in {label}."
                ));
            }
            if header.db_sources.is_empty() {
                warnings.push(format!("PEFF: Required key 'DbSource' missing in {label}."));
            }
            if header.number_of_entries == 0 {
                warnings.push(format!(
                    "PEFF: Required key 'NumberOfEntries' missing or zero in {label}."
                ));
            }
            if header.is_proteoform_db && header.has_annotation_identifiers {
                warnings.push(format!(
                    "PEFF: ProteoformDb and HasAnnotationIdentifiers MUST NOT both be true in {label}."
                ));
            }
        }
        for warning in warnings {
            self.warn(warning);
        }
        Ok(())
    }

    /// Source `parseHeaderLine_`; `true` for a `# //` block separator.
    fn parse_header_line(&mut self, line: &str, header: &mut PEFFDatabaseMetadata) -> Result<bool> {
        let trimmed = trim(line);
        if trimmed == "# //" || trimmed == "#//" {
            return Ok(true);
        }
        let Some(content) = trimmed.strip_prefix('#') else {
            return Ok(false);
        };
        let content = trim(content);
        if let Some(version) = content.strip_prefix("PEFF") {
            let version = trim(version);
            if !version.is_empty() {
                header.version = version.to_owned();
            }
            return Ok(false);
        }
        let Some((key, value)) = content.split_once('=') else {
            return Ok(false);
        };
        let (key, value) = (trim(key), trim(value));
        let flag = |value: &str| value.to_ascii_lowercase() == "true" || value == "1";
        match key {
            "DbName" => header.db_name = value.to_owned(),
            "Prefix" => header.prefix = value.to_owned(),
            "DbDescription" => header.db_description = value.to_owned(),
            "Decoy" => header.is_decoy = flag(value),
            "DbSource" => header.db_sources.push(value.to_owned()),
            "DbVersion" => header.db_version = value.to_owned(),
            "NumberOfEntries" => header.number_of_entries = to_size(value)?,
            "SequenceType" => {
                header.sequence_type = if value == "NA" {
                    PEFFSequenceType::Na
                } else {
                    PEFFSequenceType::Aa
                };
            }
            "GeneralComment" => {
                if value.is_empty() {
                    self.warn(
                        "PEFF: GeneralComment MUST NOT be empty if present (spec section 3.3.1)."
                            .to_owned(),
                    );
                } else {
                    header.general_comments.push(value.to_owned());
                }
            }
            "Conversion" => header.conversion = value.to_owned(),
            "SpecificKey" => {
                if let Some((name, description)) = value.split_once(':') {
                    header
                        .specific_keys
                        .insert(name.to_owned(), description.to_owned());
                }
            }
            "SpecificValue" => {
                if let Some((name, definition)) = value.split_once(':') {
                    header
                        .specific_values
                        .insert(name.to_owned(), definition.to_owned());
                }
            }
            "OptionalTagDef" => header.optional_tag_defs.push(value.to_owned()),
            "CustomKeyDef" => {
                let mut warnings = Vec::new();
                header
                    .custom_key_defs
                    .push(parse_custom_key_def(value, &mut warnings));
                for warning in warnings {
                    self.warn(warning);
                }
            }
            "HasAnnotationIdentifiers" => header.has_annotation_identifiers = flag(value),
            "ProteoformDb" | "ProteoformDB" | "IsProteoformDB" => {
                header.is_proteoform_db = flag(value);
            }
            _ => {
                header
                    .unrecognized_keys
                    .insert(key.to_owned(), value.to_owned());
            }
        }
        Ok(false)
    }

    /// Read the next entry; `None` at the end of the input.
    ///
    /// Source `readNext` and `readEntry_`. After leading white space, an
    /// entry is `>`, the identifier (up to the first space or tab after a
    /// non-blank character, or the line end), the description (the rest of
    /// the line) and the sequence: every following byte except space, tab
    /// and carriage return, up to a line end directly followed by `>` or the
    /// end of the input. Carriage returns are dropped everywhere. The prefix
    /// is the identifier up to its first `:`, and `sequence_length` the
    /// sequence's length unless `\Length` sets it.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the source's messages when the entry cannot be
    /// read: a line beginning with `;` (which also warns), no `>`, an empty
    /// identifier or sequence, or an input that ends inside the identifier
    /// or description line. [`Error::Parse`] also for a number that does not
    /// convert or text that is not UTF-8, and [`Error::InvalidValue`] when a
    /// bound is exceeded.
    pub fn next_entry(&mut self) -> Result<Option<PEFFEntry>> {
        if self.eof || self.at_end()? {
            return Ok(None);
        }
        let Some((identifier, description, sequence)) = self.read_entry()? else {
            return Err(parse_error(if self.entries_read == 0 {
                "The first PEFF entry could not be read!".to_owned()
            } else {
                format!(
                    "Only {} entries could be read. Parsing next record failed.",
                    self.entries_read
                )
            }));
        };
        self.entries_read += 1;
        let mut entry = PEFFEntry {
            prefix: identifier
                .split_once(':')
                .map(|(prefix, _)| prefix.to_owned())
                .unwrap_or_default(),
            sequence_length: sequence.len(),
            identifier,
            sequence,
            ..PEFFEntry::default()
        };
        self.parse_annotations(&description, &mut entry)?;
        Ok(Some(entry))
    }

    fn push_bounded(&self, buffer: &mut Vec<u8>, byte: u8, bound: usize, what: &str) -> Result<()> {
        if buffer.len() == bound {
            return Err(limit(what));
        }
        buffer.push(byte);
        Ok(())
    }

    fn read_entry(&mut self) -> Result<Option<(String, String, String)>> {
        while let Some(c) = self.peek()? {
            if matches!(c, b' ' | b'\t' | b'\n' | b'\r') {
                self.bump()?;
            } else {
                break;
            }
        }
        if self.peek()? == Some(b';') {
            self.warn("PEFF: Lines beginning with ';' are not permitted.".to_owned());
            while let Some(c) = self.bump()? {
                if c == b'\n' {
                    break;
                }
            }
            return Ok(None);
        }
        if self.bump()? != Some(b'>') {
            return Ok(None);
        }
        let line_bound = self.limits.max_line_bytes;
        let mut identifier = Vec::new();
        let mut has_description = true;
        loop {
            match self.bump()? {
                Some(b' ' | b'\t') => {
                    if !identifier.is_empty() {
                        break;
                    }
                }
                Some(b'\n') => {
                    has_description = false;
                    break;
                }
                Some(b'\r') => {}
                None => {
                    self.eof = true;
                    return Ok(None);
                }
                Some(byte) => self.push_bounded(&mut identifier, byte, line_bound, "line byte")?,
            }
        }
        if identifier.is_empty() {
            return Ok(None);
        }
        let mut description = Vec::new();
        if has_description {
            loop {
                match self.bump()? {
                    Some(b'\n') => break,
                    Some(b'\r') => {}
                    None => {
                        self.eof = true;
                        return Ok(None);
                    }
                    Some(byte) => {
                        let bound = line_bound.saturating_sub(identifier.len());
                        self.push_bounded(&mut description, byte, bound, "line byte")?;
                    }
                }
            }
        }
        let mut sequence = Vec::new();
        loop {
            match self.bump()? {
                Some(b'\n') => {
                    if self.peek()? == Some(b'>') {
                        break;
                    }
                }
                Some(b'\r' | b' ' | b'\t') => {}
                None => {
                    self.eof = true;
                    break;
                }
                Some(byte) => {
                    let bound = self.limits.max_sequence_bytes;
                    self.push_bounded(&mut sequence, byte, bound, "sequence byte")?;
                }
            }
        }
        if sequence.is_empty() {
            return Ok(None);
        }
        Ok(Some((
            utf8(identifier, "identifier")?,
            utf8(description, "description line")?,
            utf8(sequence, "sequence")?,
        )))
    }

    /// Source `parseAnnotations_`.
    fn parse_annotations(&mut self, description: &str, entry: &mut PEFFEntry) -> Result<()> {
        if description.is_empty() {
            return Ok(());
        }
        let bytes = description.as_bytes();
        let Some(first) = next_key_separator(bytes, 0) else {
            // No annotations: the whole description is the protein name.
            entry.protein_names.push(trim(description).to_owned());
            return Ok(());
        };
        let before = trim(&description[..first]);
        if !before.is_empty() {
            entry.protein_names.push(before.to_owned());
        }

        // The header of this entry's database, for CustomKeyDef validation.
        let by_prefix = if entry.prefix.is_empty() {
            None
        } else {
            self.headers.iter().position(|h| h.prefix == entry.prefix)
        };
        let header = by_prefix.or_else(|| (self.headers.len() == 1).then_some(0));
        let nucleotide =
            header.is_some_and(|index| self.headers[index].sequence_type == PEFFSequenceType::Na);
        let sequence_length = entry.sequence.len();
        let identifier = entry.identifier.clone();
        let mut items = 0usize;
        let max_items = self.limits.max_items_per_entry;

        let mut seen = BTreeSet::new();
        let mut previous = first;
        while previous < bytes.len() {
            previous += 1;
            let position = next_key_separator(bytes, previous).unwrap_or(bytes.len());
            let annotation = trim(&description[previous..position]);
            previous = position;
            if annotation.is_empty() {
                continue;
            }
            let Some((key, value)) = annotation.split_once('=') else {
                continue;
            };
            if !seen.insert(key.to_owned()) {
                self.warn(format!(
                    "PEFF: Duplicate key '{key}' in entry '{identifier}'."
                ));
            }
            let mut count = |n: usize| -> Result<()> {
                items = items.saturating_add(n);
                if items > max_items {
                    Err(limit("annotation item"))
                } else {
                    Ok(())
                }
            };
            match key {
                "PName" => {
                    let value = trim(value);
                    if is_list_form(value.as_bytes()) {
                        let names = parse_paren_list(value);
                        count(names.len())?;
                        entry
                            .protein_names
                            .extend(names.into_iter().map(peff_unescape));
                    } else if !value.is_empty() {
                        count(1)?;
                        entry.protein_names.push(peff_unescape(value));
                    }
                }
                "GName" => entry.gene_name = value.to_owned(),
                "NcbiTaxId" | "OX" => entry.ncbi_tax_id = to_int32(value)?,
                "TaxName" => entry.taxonomy_name = value.to_owned(),
                "Length" => entry.sequence_length = to_size(value)?,
                "SV" => entry.sequence_version = value.to_owned(),
                "EV" => entry.entry_version = value.to_owned(),
                "PE" => entry.protein_existence = to_int32(value)?,
                "DbUniqueId" => entry.db_unique_id = value.to_owned(),
                "ID" => entry.entry_id = value.to_owned(),
                "AltAC" => {
                    let accessions = parse_paren_list(value);
                    count(accessions.len())?;
                    entry
                        .alt_accessions
                        .extend(accessions.into_iter().map(peff_unescape));
                }
                "ModRes" | "ModResPsi" | "ModResUnimod" => {
                    let tuples = parse_paren_list(value);
                    count(tuples.len())?;
                    self.modifications(key, &tuples, nucleotide, sequence_length, entry)?;
                }
                "VariantSimple" => {
                    let tuples = parse_paren_list(value);
                    count(tuples.len())?;
                    for tuple in &tuples {
                        let variant =
                            self.variant_simple(tuple, nucleotide, sequence_length, &identifier)?;
                        entry.simple_variants.push(variant);
                    }
                }
                "VariantComplex" => {
                    let tuples = parse_paren_list(value);
                    count(tuples.len())?;
                    for tuple in &tuples {
                        let variant =
                            self.variant_complex(tuple, nucleotide, sequence_length, &identifier)?;
                        entry.complex_variants.push(variant);
                    }
                }
                "Variant" => {
                    self.warn(format!(
                        "PEFF: Deprecated key 'Variant' encountered in entry '{identifier}'. \
                         Readers should use VariantSimple/VariantComplex."
                    ));
                    entry
                        .custom_annotations
                        .insert(key.to_owned(), value.to_owned());
                }
                "Processed" => {
                    let tuples = parse_paren_list(value);
                    count(tuples.len())?;
                    for tuple in &tuples {
                        let region = self.processed_region(tuple, sequence_length, &identifier)?;
                        entry.processed_regions.push(region);
                    }
                }
                "DisulfideBond" => {
                    let tuples = parse_paren_list(value);
                    count(tuples.len())?;
                    for tuple in &tuples {
                        entry.disulfide_bonds.push(parse_disulfide_bond(tuple)?);
                    }
                }
                "Proteoform" => {
                    let forms = parse_paren_list(value);
                    count(forms.len())?;
                    entry
                        .proteoforms
                        .extend(forms.into_iter().map(peff_unescape));
                }
                _ => {
                    if let Some(index) = header {
                        self.validate_custom(index, key, value, &identifier)?;
                    }
                    count(1)?;
                    entry
                        .custom_annotations
                        .insert(key.to_owned(), value.to_owned());
                }
            }
        }
        Ok(())
    }

    /// The `ModRes`, `ModResPsi` and `ModResUnimod` branch of
    /// `parseAnnotations_`.
    fn modifications(
        &mut self,
        key: &str,
        tuples: &[&str],
        nucleotide: bool,
        sequence_length: usize,
        entry: &mut PEFFEntry,
    ) -> Result<()> {
        let identifier = entry.identifier.clone();
        if nucleotide {
            self.warn(format!(
                "PEFF: {key} annotations are not applicable to nucleotide databases (entry '{identifier}')."
            ));
        }
        let typed = key == "ModResPsi" || key == "ModResUnimod";
        for tuple in tuples {
            let parts = split_by_pipe_escape_aware(tuple);
            if parts.len() < 3 {
                self.warn(format!(
                    "PEFF: {key} tuple must have at least 3 fields (position|accession|name) in entry '{identifier}'."
                ));
            } else if parts[2].is_empty() {
                self.warn(format!(
                    "PEFF: {key} tuple is missing required name field in entry '{identifier}'."
                ));
            }
            if typed && parts.len() >= 2 && parts[1].is_empty() {
                self.warn(format!(
                    "PEFF: {key} tuple is missing required accession field in entry '{identifier}'."
                ));
            }
            if parts.len() >= 2 && !parts[1].is_empty() {
                if key == "ModResPsi" && !parts[1].starts_with("MOD:") {
                    self.warn(format!(
                        "PEFF: {key} accession should start with 'MOD:' in entry '{identifier}'."
                    ));
                } else if key == "ModResUnimod" && !parts[1].starts_with("UNIMOD:") {
                    self.warn(format!(
                        "PEFF: {key} accession should start with 'UNIMOD:' in entry '{identifier}'."
                    ));
                }
            }
            // Validate the known positions.
            let (_, positions) = extract_annotation_id(&parts[0])?;
            let tokens = if positions.contains(',') {
                split(positions, ',')
            } else {
                vec![positions]
            };
            for token in tokens {
                if token == "?" || token.is_empty() {
                    continue;
                }
                if !is_unsigned_integer(token) || token == "0" {
                    self.warn(format!(
                        "PEFF: {key} position must be a positive integer or '?' in entry '{identifier}'."
                    ));
                    continue;
                }
                let position = to_size(token)?;
                if sequence_length > 0 && position > sequence_length {
                    self.warn(format!(
                        "PEFF: {key} position {position} exceeds sequence length {sequence_length} in entry '{identifier}'."
                    ));
                }
            }

            let mut modification = self.parse_modification(tuple)?;
            if key == "ModResPsi" {
                modification.modification_type = PEFFModificationType::PsiMod;
            } else if key == "ModResUnimod" {
                modification.modification_type = PEFFModificationType::Unimod;
            }
            if modification.modification_type != PEFFModificationType::Generic
                && modification.name.is_empty()
            {
                self.warn(format!(
                    "PEFF: Modification with accession '{}' is missing required name field.",
                    modification.accession
                ));
            }

            // Comma-separated positions expand into one modification each.
            let (annotation_id, positions) = extract_annotation_id(&parts[0])?;
            if positions.contains(',') {
                for token in split(positions, ',') {
                    let unknown = token == "?" || token.is_empty();
                    if !unknown && (!is_unsigned_integer(token) || token == "0") {
                        self.warn(format!(
                            "PEFF: {key} position must be a positive integer or '?' in entry '{identifier}'."
                        ));
                    }
                    let mut expanded = modification.clone();
                    expanded.annotation_id = annotation_id;
                    expanded.position = if unknown { 0 } else { to_size(token)? };
                    if expanded.position > 0
                        && sequence_length > 0
                        && expanded.position > sequence_length
                    {
                        self.warn(format!(
                            "PEFF: {key} position {} exceeds sequence length {sequence_length} in entry '{identifier}'.",
                            expanded.position
                        ));
                    }
                    entry.modifications.push(expanded);
                }
                continue;
            }
            if modification.position > 0
                && sequence_length > 0
                && modification.position > sequence_length
            {
                self.warn(format!(
                    "PEFF: {key} position {} exceeds sequence length {sequence_length} in entry '{identifier}'.",
                    modification.position
                ));
            }
            entry.modifications.push(modification);
        }
        Ok(())
    }

    /// Source `parseModification_`: `[id:]pos|accession|name|OptionalTag`,
    /// the position `?` or empty for unknown; with comma-separated positions
    /// the first is stored.
    fn parse_modification(&mut self, tuple: &str) -> Result<PEFFModification> {
        let mut modification = PEFFModification::default();
        let parts = split_by_pipe_escape_aware(tuple);
        if parts.len() < 2 {
            return Ok(modification);
        }
        let (annotation_id, positions) = extract_annotation_id(&parts[0])?;
        modification.annotation_id = annotation_id;
        if positions == "?" || positions.is_empty() {
            modification.position = 0;
        } else if positions.contains(',') {
            // The source converts the first token as it is, so "?,5" is a
            // conversion error there as well.
            if let Some(first) = split(positions, ',').first().filter(|t| !t.is_empty()) {
                modification.position = to_size(first)?;
            }
        } else {
            modification.position = to_size(positions)?;
        }
        modification.modification_type = PEFFModificationType::of_accession(&parts[1]);
        modification.accession = parts[1].clone();
        if parts.len() >= 3 {
            modification.name = parts[2].clone();
        }
        if parts.len() >= 4 {
            modification.optional_tag = parts[3].clone();
        }
        if parts.len() > 4 {
            self.warn(
                "PEFF: Modification tuple has more than 4 fields; treating last field as OptionalTag."
                    .to_owned(),
            );
            modification.optional_tag = parts[parts.len() - 1].clone();
        }
        Ok(modification)
    }

    /// The `VariantSimple` branch with `parseVariantSimple_`:
    /// `[id:]pos|aa|OptionalTag`.
    fn variant_simple(
        &mut self,
        tuple: &str,
        nucleotide: bool,
        sequence_length: usize,
        identifier: &str,
    ) -> Result<PEFFVariantSimple> {
        let parts = split_by_pipe_escape_aware(tuple);
        if parts.len() >= 2 && parts[1].len() != 1 {
            self.warn(format!(
                "PEFF: VariantSimple newAminoAcid must be exactly 1 character in entry '{identifier}'."
            ));
        }
        let (annotation_id, position) = extract_annotation_id(&parts[0])?;
        if !is_unsigned_integer(position) || position == "0" {
            self.warn(format!(
                "PEFF: VariantSimple position must be a positive integer in entry '{identifier}'."
            ));
        }
        let mut variant = PEFFVariantSimple::default();
        if parts.len() >= 2 {
            variant.annotation_id = annotation_id;
            variant.position = to_size(position)?;
            variant.variant_aa = parts[1].chars().next().unwrap_or('\0');
            if parts.len() >= 3 {
                variant.optional_tag = parts[2].clone();
            }
        }
        if variant.position == 0 {
            self.warn(format!(
                "PEFF: VariantSimple position must be > 0 in entry '{identifier}'."
            ));
        }
        if sequence_length > 0 && variant.position > sequence_length {
            self.warn(format!(
                "PEFF: VariantSimple position {} exceeds sequence length {sequence_length} in entry '{identifier}'.",
                variant.position
            ));
        }
        if variant.variant_aa == '\0' {
            self.warn(format!(
                "PEFF: VariantSimple variant_aa must not be empty in entry '{identifier}'."
            ));
        } else if !valid_variant_code(variant.variant_aa, nucleotide) {
            self.warn(format!(
                "PEFF: VariantSimple newAminoAcid '{}' is not valid in entry '{identifier}'.",
                variant.variant_aa
            ));
        }
        Ok(variant)
    }

    /// The `VariantComplex` branch with `parseVariantComplex_`:
    /// `[id:]start|end|replacement|OptionalTag`.
    fn variant_complex(
        &mut self,
        tuple: &str,
        nucleotide: bool,
        sequence_length: usize,
        identifier: &str,
    ) -> Result<PEFFVariantComplex> {
        let parts = split_by_pipe_escape_aware(tuple);
        if parts.len() < 3 {
            self.warn(format!(
                "PEFF: VariantComplex tuple must have at least 3 fields (start|end|newSequence) in entry '{identifier}'."
            ));
        }
        // The source extracts the annotation identifier only from tuples of
        // two or more fields, so only those can fail its conversion.
        let (annotation_id, start) = if parts.len() >= 2 {
            extract_annotation_id(&parts[0])?
        } else {
            (super::ANNOTATION_ID_NOT_SET, parts[0].as_str())
        };
        if parts.len() >= 2
            && (!is_unsigned_integer(start)
                || start == "0"
                || !is_unsigned_integer(&parts[1])
                || parts[1] == "0")
        {
            self.warn(format!(
                "PEFF: VariantComplex positions must be positive integers in entry '{identifier}'."
            ));
        }
        let mut variant = PEFFVariantComplex::default();
        if parts.len() >= 3 {
            variant.annotation_id = annotation_id;
            variant.start_position = to_size(start)?;
            variant.end_position = to_size(&parts[1])?;
            variant.replacement = parts[2].clone();
            if parts.len() >= 4 {
                variant.optional_tag = parts[3].clone();
            }
        }
        if variant.start_position == variant.end_position && variant.replacement.len() == 1 {
            self.warn(format!(
                "PEFF: VariantComplex at position {} with single AA replacement must be encoded as VariantSimple in entry '{identifier}'.",
                variant.start_position
            ));
        }
        if variant.start_position == 0 || variant.end_position == 0 {
            self.warn(format!(
                "PEFF: VariantComplex positions must be > 0 in entry '{identifier}'."
            ));
        }
        if variant.start_position > variant.end_position {
            self.warn(format!(
                "PEFF: VariantComplex startPosition must be <= endPosition in entry '{identifier}'."
            ));
        }
        if sequence_length > 0
            && (variant.start_position > sequence_length || variant.end_position > sequence_length)
        {
            self.warn(format!(
                "PEFF: VariantComplex positions exceed sequence length {sequence_length} in entry '{identifier}'."
            ));
        }
        if let Some(bad) = variant
            .replacement
            .chars()
            .find(|&c| !valid_variant_code(c, nucleotide))
        {
            self.warn(format!(
                "PEFF: VariantComplex newSequence contains invalid character '{bad}' in entry '{identifier}'."
            ));
        }
        Ok(variant)
    }

    /// The `Processed` branch with `parseProcessedRegion_`:
    /// `[id:]start|end|accession|name|OptionalTag`.
    fn processed_region(
        &mut self,
        tuple: &str,
        sequence_length: usize,
        identifier: &str,
    ) -> Result<PEFFProcessedRegion> {
        let parts = split_by_pipe_escape_aware(tuple);
        let (annotation_id, start) = if parts.len() >= 2 {
            extract_annotation_id(&parts[0])?
        } else {
            (super::ANNOTATION_ID_NOT_SET, parts[0].as_str())
        };
        if parts.len() >= 2
            && (!is_unsigned_integer(start)
                || start == "0"
                || !is_unsigned_integer(&parts[1])
                || parts[1] == "0")
        {
            self.warn(format!(
                "PEFF: Processed positions must be positive integers in entry '{identifier}'."
            ));
        }
        if parts.len() < 4 {
            self.warn(format!(
                "PEFF: Processed tuple must have at least 4 fields (start|end|accession|name) in entry '{identifier}'."
            ));
        } else {
            if parts[2].is_empty() {
                self.warn(format!(
                    "PEFF: Processed tuple is missing required accession field in entry '{identifier}'."
                ));
            }
            if parts[3].is_empty() {
                self.warn(format!(
                    "PEFF: Processed tuple is missing required name field in entry '{identifier}'."
                ));
            }
        }
        let mut region = PEFFProcessedRegion::default();
        if parts.len() >= 3 {
            region.annotation_id = annotation_id;
            region.start_position = to_size(start)?;
            region.end_position = to_size(&parts[1])?;
            region.accession = parts[2].clone();
            if parts.len() >= 4 {
                region.name = parts[3].clone();
            }
            if parts.len() >= 5 {
                region.optional_tag = parts[4].clone();
            }
        }
        if region.start_position == 0 || region.end_position == 0 {
            self.warn(format!(
                "PEFF: Processed positions must be > 0 in entry '{identifier}'."
            ));
        }
        if region.start_position > region.end_position {
            self.warn(format!(
                "PEFF: Processed startPosition must be <= endPosition in entry '{identifier}'."
            ));
        }
        if sequence_length > 0
            && (region.start_position > sequence_length || region.end_position > sequence_length)
        {
            self.warn(format!(
                "PEFF: Processed positions exceed sequence length {sequence_length} in entry '{identifier}'."
            ));
        }
        if !region.accession.is_empty() && !region.accession.starts_with("PEFF:") {
            self.warn(format!(
                "PEFF: Processed accession should start with 'PEFF:' in entry '{identifier}'."
            ));
        }
        Ok(region)
    }

    /// The custom-key branch: validation against the `CustomKeyDef` of the
    /// entry's database, as warnings.
    fn validate_custom(
        &mut self,
        header: usize,
        key: &str,
        value: &str,
        identifier: &str,
    ) -> Result<()> {
        let Some(definition) = self.headers[header]
            .custom_key_defs
            .iter()
            .find(|d| d.key_name == key)
            .cloned()
        else {
            return Ok(());
        };
        let expected = if definition.field_names.is_empty() {
            definition.field_types.len()
        } else {
            definition.field_names.len()
        };
        if !definition.field_names.is_empty()
            && !definition.field_types.is_empty()
            && definition.field_names.len() != definition.field_types.len()
        {
            self.warn(format!(
                "PEFF: CustomKeyDef for '{key}' has mismatched FieldNames/FieldTypes lengths."
            ));
        }
        let optional_tag_last = definition
            .field_names
            .last()
            .is_some_and(|name| name == "OptionalTag");
        let expression = if definition.regexp.is_empty() {
            None
        } else {
            match BoostRegex::new(&definition.regexp) {
                Ok(expression) => Some(expression),
                Err(_) => {
                    self.warn(format!(
                        "PEFF: CustomKeyDef RegExp for '{key}' is not a valid C++ regex."
                    ));
                    None
                }
            }
        };
        for item in parse_paren_list(value) {
            if let Some(expression) = &expression {
                if !expression.is_full_match(item.as_bytes())? {
                    self.warn(format!(
                        "PEFF: Value of custom key '{key}' does not match RegExp in entry '{identifier}'."
                    ));
                }
            }
            let components = split_by_pipe_escape_aware(item);
            if expected > 0
                && components.len() != expected
                && !(optional_tag_last && components.len() + 1 == expected)
            {
                self.warn(format!(
                    "PEFF: Custom key '{key}' has {} component(s), expected {expected} in entry '{identifier}'.",
                    components.len()
                ));
            }
            for (index, (field_type, component)) in
                definition.field_types.iter().zip(&components).enumerate()
            {
                if !validate_custom_key_field_type(field_type, component) {
                    self.warn(format!(
                        "PEFF: Custom key '{key}' field {} does not match type '{field_type}' in entry '{identifier}'.",
                        index + 1
                    ));
                }
            }
        }
        Ok(())
    }
}

/// `find_next_key_sep` of `parseAnnotations_`: the next unescaped backslash
/// followed by an ASCII letter, digit or underscore.
fn next_key_separator(bytes: &[u8], from: usize) -> Option<usize> {
    let mut i = from;
    while i < bytes.len() {
        if bytes[i] != b'\\' || i + 1 >= bytes.len() {
            i += 1;
            continue;
        }
        let next = bytes[i + 1];
        if matches!(next, b'\\' | b'|' | b'(' | b')') {
            i += 2;
            continue;
        }
        if next.is_ascii_alphanumeric() || next == b'_' {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// `isListForm` of the `PName` branch: the value is one or more adjacent,
/// balanced `(…)` groups with nothing between them.
fn is_list_form(bytes: &[u8]) -> bool {
    if bytes.first() != Some(&b'(') || bytes.last() != Some(&b')') {
        return false;
    }
    let n = bytes.len();
    let mut depth = 0i64;
    let mut zero_in_interior = false;
    let mut i = 0;
    while i < n {
        if bytes[i] == b'\\' && i + 1 < n {
            i += 2;
            continue;
        }
        if bytes[i] == b'(' {
            depth += 1;
        } else if bytes[i] == b')' {
            depth -= 1;
        }
        if depth == 0 && i + 1 < n {
            zero_in_interior = true;
        }
        if depth < 0 {
            return false;
        }
        i += 1;
    }
    if depth != 0 {
        return false;
    }
    if !zero_in_interior {
        return true;
    }
    i = 0;
    while i < n {
        if bytes[i] == b'\\' && i + 1 < n {
            i += 2;
            continue;
        }
        if bytes[i] == b'(' {
            depth += 1;
        } else if bytes[i] == b')' {
            depth -= 1;
            if depth == 0 && i + 1 < n && bytes[i + 1] != b'(' {
                return false;
            }
        }
        i += 1;
    }
    true
}

/// Source `parseParenList_`: the contents of each top-level `(…)` group,
/// escapes kept; with no group, the whole value. An unclosed group runs to
/// the end; an empty group yields nothing.
pub(super) fn parse_paren_list(value: &str) -> Vec<&str> {
    let bytes = value.as_bytes();
    let mut result = Vec::new();
    let mut position = 0;
    while position < bytes.len() {
        if bytes[position] != b'(' {
            position += 1;
            continue;
        }
        let mut depth = 1usize;
        let start = position + 1;
        let mut end = start;
        while end < bytes.len() && depth > 0 {
            if bytes[end] == b'\\' && end + 1 < bytes.len() {
                end += 2;
                continue;
            }
            if bytes[end] == b'(' {
                depth += 1;
            } else if bytes[end] == b')' {
                depth -= 1;
            }
            if depth > 0 {
                end += 1;
            }
        }
        // `end` is at the closing ')' or at the end of the value, both
        // character boundaries; an escape step can overshoot only to the end.
        let end = end.min(bytes.len());
        if end > start {
            result.push(&value[start..end]);
        }
        position = end + 1;
    }
    if result.is_empty() && !value.is_empty() {
        result.push(value);
    }
    result
}

/// Source `parseDisulfideBond_`: `[id:]id1,id2|description`, split at the
/// first `|` without escape handling; the description is un-escaped.
fn parse_disulfide_bond(tuple: &str) -> Result<PEFFDisulfideBond> {
    let mut bond = PEFFDisulfideBond::default();
    let identifiers = match tuple.split_once('|') {
        Some((identifiers, description)) => {
            bond.optional_tag = peff_unescape(description);
            identifiers
        }
        None => tuple,
    };
    let (annotation_id, rest) = extract_annotation_id(identifiers)?;
    bond.annotation_id = annotation_id;
    if let Some((first, second)) = rest.split_once(',') {
        bond.id1 = first.to_owned();
        bond.id2 = second.to_owned();
    }
    Ok(bond)
}

/// `is_valid_variant_code`: `*`, or an ASCII letter; in a nucleotide
/// database only the IUPAC nucleotide codes.
fn valid_variant_code(c: char, nucleotide: bool) -> bool {
    if c == '*' {
        return true;
    }
    if !c.is_ascii_alphabetic() {
        return false;
    }
    !nucleotide || "ACGTURYSWKMBDHVN".contains(c.to_ascii_uppercase())
}

/// Source `isPEFFFile` over an open input.
///
/// Reads lines until one decides: a `#` line whose content (after `#` and
/// spaces or tabs) starts with `PEFF ` and a digit is PEFF, other `#` lines
/// are skipped, a `>` line is PEFF exactly when it holds a backslash, and any
/// other non-blank line is not. At most 64 MiB are inspected; a longer
/// undecided prefix is not PEFF.
pub(super) fn is_peff(mut input: impl BufRead) -> bool {
    const MAX_SCAN: usize = 64 * 1024 * 1024;
    let mut scanned = 0usize;
    let mut line = Vec::new();
    loop {
        line.clear();
        let remaining = MAX_SCAN.saturating_sub(scanned);
        let read = match input
            .by_ref()
            .take(remaining as u64 + 1)
            .read_until(b'\n', &mut line)
        {
            Ok(read) => read,
            Err(_) => return false,
        };
        if read == 0 {
            return false;
        }
        scanned += read;
        if scanned > MAX_SCAN {
            return false;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        let Some(start) = line
            .iter()
            .position(|b| !matches!(b, b' ' | b'\t' | b'\r' | b'\n'))
        else {
            continue;
        };
        match line[start] {
            b'#' => {
                let mut content = start + 1;
                while matches!(line.get(content), Some(b' ' | b'\t')) {
                    content += 1;
                }
                if content + 5 < line.len()
                    && &line[content..content + 5] == b"PEFF "
                    && line[content + 5].is_ascii_digit()
                {
                    return true;
                }
            }
            b'>' => return line.contains(&b'\\'),
            _ => return false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paren_lists_follow_the_source_scanner() {
        assert_eq!(parse_paren_list("(a)(b|c)"), ["a", "b|c"]);
        assert_eq!(parse_paren_list("(a\\)b)"), ["a\\)b"]);
        assert_eq!(parse_paren_list("(a(b)c)"), ["a(b)c"]);
        assert_eq!(parse_paren_list("single"), ["single"]);
        assert_eq!(parse_paren_list("()"), ["()"]);
        assert_eq!(parse_paren_list("(open"), ["open"]);
        assert_eq!(parse_paren_list("(é\\"), ["é\\"]);
        assert!(parse_paren_list("").is_empty());
    }

    #[test]
    fn list_form_needs_adjacent_balanced_groups() {
        for (text, expected) in [
            ("(a)", true),
            ("(a)(b)", true),
            ("(R)-amino acid hydroxylase", false),
            ("Insulin (Fragment)", false),
            ("(a) (b)", false),
            ("(a\\))", true),
            ("((a)", false),
        ] {
            assert_eq!(is_list_form(text.as_bytes()), expected, "{text}");
        }
    }

    #[test]
    fn key_separators_skip_escapes_and_need_a_key_character() {
        let text = br"name \| x \(y\) \PName=1 \ \GName=2";
        let first = next_key_separator(text, 0).unwrap();
        assert_eq!(&text[first..first + 6], br"\PName");
        let second = next_key_separator(text, first + 1).unwrap();
        assert_eq!(&text[second..second + 6], br"\GName");
    }

    #[test]
    fn disulfide_bonds_split_at_the_first_pipe() {
        let bond = parse_disulfide_bond(r"4:1,2|between \| chains").unwrap();
        assert_eq!(bond.annotation_id, 4);
        assert_eq!((bond.id1.as_str(), bond.id2.as_str()), ("1", "2"));
        assert_eq!(bond.optional_tag, "between | chains");
        let lone = parse_disulfide_bond("7").unwrap();
        assert!(lone.id1.is_empty() && lone.id2.is_empty());
    }
}

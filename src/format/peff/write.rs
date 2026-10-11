// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PEFF writing: `formatHeader_`, `appendDatabaseDescriptionBlock_`,
//! `formatEntry_` and `toProForma` (`PEFFFile.cpp:1810-1885`, `2960-3462`).

use super::text::{escape_custom_key_def_quoted_value, escape_peff};
use super::{
    ANNOTATION_ID_NOT_SET, PEFFDatabaseMetadata, PEFFEntry, PEFFModification, PEFFModificationType,
    PEFFSequenceType,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Source `appendDatabaseDescriptionBlock_`: the required-key warnings, then
/// the block. Optional keys are written only when set; `Decoy` and
/// `SequenceType` always are.
fn database_block(out: &mut String, header: &PEFFDatabaseMetadata, warnings: &mut Vec<String>) {
    if header.db_name.is_empty() {
        warnings.push("PEFF: Required key 'DbName' is empty on write.".to_owned());
    }
    if header.prefix.is_empty() {
        warnings.push("PEFF: Required key 'Prefix' is empty on write.".to_owned());
    }
    if header.db_version.is_empty() {
        warnings.push("PEFF: Required key 'DbVersion' is empty on write.".to_owned());
    }
    if header.db_sources.is_empty() {
        warnings.push("PEFF: Required key 'DbSource' is missing on write.".to_owned());
    }
    if header.number_of_entries == 0 {
        warnings.push("PEFF: Required key 'NumberOfEntries' is zero on write.".to_owned());
    }
    // `write!` into a String cannot fail.
    if !header.db_name.is_empty() {
        let _ = writeln!(out, "# DbName={}", header.db_name);
    }
    if !header.prefix.is_empty() {
        let _ = writeln!(out, "# Prefix={}", header.prefix);
    }
    if !header.db_description.is_empty() {
        let _ = writeln!(out, "# DbDescription={}", header.db_description);
    }
    let _ = writeln!(
        out,
        "# Decoy={}",
        if header.is_decoy { "true" } else { "false" }
    );
    for source in &header.db_sources {
        let _ = writeln!(out, "# DbSource={source}");
    }
    if !header.db_version.is_empty() {
        let _ = writeln!(out, "# DbVersion={}", header.db_version);
    }
    if header.number_of_entries > 0 {
        let _ = writeln!(out, "# NumberOfEntries={}", header.number_of_entries);
    }
    let _ = writeln!(
        out,
        "# SequenceType={}",
        match header.sequence_type {
            PEFFSequenceType::Aa => "AA",
            PEFFSequenceType::Na => "NA",
        }
    );
    if !header.conversion.is_empty() {
        let _ = writeln!(out, "# Conversion={}", header.conversion);
    }
    for (key, description) in &header.specific_keys {
        let _ = writeln!(out, "# SpecificKey={key}:{description}");
    }
    for (key, definition) in &header.specific_values {
        let _ = writeln!(out, "# SpecificValue={key}:{definition}");
    }
    for definition in &header.optional_tag_defs {
        let _ = writeln!(out, "# OptionalTagDef={definition}");
    }
    for definition in &header.custom_key_defs {
        let _ = write!(out, "# CustomKeyDef=(KeyName={}", definition.key_name);
        if !definition.description.is_empty() {
            let _ = write!(
                out,
                "|Description=\"{}\"",
                escape_custom_key_def_quoted_value(&definition.description)
            );
        }
        if !definition.concept_curie.is_empty() {
            let _ = write!(out, "|ConceptCURIE={}", definition.concept_curie);
        }
        if !definition.regexp.is_empty() {
            let _ = write!(
                out,
                "|RegExp=\"{}\"",
                escape_custom_key_def_quoted_value(&definition.regexp)
            );
        }
        if !definition.field_names.is_empty() {
            let _ = write!(out, "|FieldNames={}", definition.field_names.join(","));
        }
        if !definition.field_types.is_empty() {
            let _ = write!(out, "|FieldTypes={}", definition.field_types.join(","));
        }
        out.push_str(")\n");
    }
    if header.is_proteoform_db {
        out.push_str("# ProteoformDb=true\n");
    }
    if header.has_annotation_identifiers {
        out.push_str("# HasAnnotationIdentifiers=true\n");
    }
    for (key, value) in &header.unrecognized_keys {
        let _ = writeln!(out, "# {key}={value}");
    }
    out.push_str("# //\n");
}

/// Source `formatHeader_(header)`: the file description block (version,
/// general comments) and one database block.
pub(super) fn format_header(header: &PEFFDatabaseMetadata, warnings: &mut Vec<String>) -> String {
    let mut out = String::new();
    let version = if header.version.is_empty() {
        warnings.push("PEFF: Missing PEFF version on write; defaulting to 1.0.".to_owned());
        "1.0"
    } else {
        header.version.as_str()
    };
    let _ = writeln!(out, "# PEFF {version}");
    for comment in header.general_comments.iter().filter(|c| !c.is_empty()) {
        let _ = writeln!(out, "# GeneralComment={comment}");
    }
    out.push_str("# //\n");
    database_block(&mut out, header, warnings);
    out
}

/// Source `formatHeader_(headers)`: one file description block with the
/// first non-empty version and the merged general comments, then one
/// database block per header.
pub(super) fn format_headers(
    headers: &[PEFFDatabaseMetadata],
    warnings: &mut Vec<String>,
) -> Result<String> {
    if headers.is_empty() {
        return Err(Error::InvalidValue(
            "PEFF: formatHeader() requires at least one database header".into(),
        ));
    }
    let version = match headers.iter().find(|h| !h.version.is_empty()) {
        Some(header) => header.version.as_str(),
        None => {
            warnings.push("PEFF: Missing PEFF version on write; defaulting to 1.0.".to_owned());
            "1.0"
        }
    };
    if headers
        .iter()
        .any(|h| !h.version.is_empty() && h.version != version)
    {
        warnings.push(format!(
            "PEFF: Multiple PEFF version values provided; writing '{version}'."
        ));
    }
    let mut seen = BTreeSet::new();
    let mut out = String::new();
    let _ = writeln!(out, "# PEFF {version}");
    for comment in headers.iter().flat_map(|h| &h.general_comments) {
        if !comment.is_empty() && seen.insert(comment.as_str()) {
            let _ = writeln!(out, "# GeneralComment={comment}");
        }
    }
    out.push_str("# //\n");
    for header in headers {
        database_block(&mut out, header, warnings);
    }
    Ok(out)
}

fn annotation_id(out: &mut String, id: u32) {
    if id != ANNOTATION_ID_NOT_SET {
        let _ = write!(out, "{id}:");
    }
}

fn modifications(out: &mut String, modifications: &[&PEFFModification], key: &str) {
    if modifications.is_empty() {
        return;
    }
    let _ = write!(out, " \\{key}=");
    for modification in modifications {
        out.push('(');
        annotation_id(out, modification.annotation_id);
        if modification.position == 0 {
            out.push('?');
        } else {
            let _ = write!(out, "{}", modification.position);
        }
        let _ = write!(out, "|{}", escape_peff(&modification.accession));
        if !modification.name.is_empty() || !modification.optional_tag.is_empty() {
            let _ = write!(out, "|{}", escape_peff(&modification.name));
        }
        if !modification.optional_tag.is_empty() {
            let _ = write!(out, "|{}", escape_peff(&modification.optional_tag));
        }
        out.push(')');
    }
}

/// Source `formatEntry_`: the description line with every set annotation in
/// a fixed key order, then the sequence in lines of 80.
///
/// Modifications are grouped by type into `\ModResPsi`, `\ModResUnimod` and
/// `\ModRes`, keeping their order within each; custom annotations are
/// written verbatim, in key order. Text fields are PEFF-escaped except the
/// identifier, the disulfide identifiers' separator and custom values.
pub(super) fn format_entry(entry: &PEFFEntry) -> String {
    let mut out = String::with_capacity(entry.sequence.len() + entry.sequence.len() / 80 + 128);
    out.push('>');
    out.push_str(&entry.identifier);
    if !entry.protein_names.is_empty() {
        out.push_str(" \\PName=");
        for name in &entry.protein_names {
            let _ = write!(out, "({})", escape_peff(name));
        }
    }
    if !entry.gene_name.is_empty() {
        let _ = write!(out, " \\GName={}", escape_peff(&entry.gene_name));
    }
    if entry.ncbi_tax_id != 0 {
        let _ = write!(out, " \\NcbiTaxId={}", entry.ncbi_tax_id);
    }
    if !entry.taxonomy_name.is_empty() {
        let _ = write!(out, " \\TaxName={}", escape_peff(&entry.taxonomy_name));
    }
    if entry.sequence_length > 0 {
        let _ = write!(out, " \\Length={}", entry.sequence_length);
    }
    if !entry.sequence_version.is_empty() {
        let _ = write!(out, " \\SV={}", escape_peff(&entry.sequence_version));
    }
    if !entry.entry_version.is_empty() {
        let _ = write!(out, " \\EV={}", escape_peff(&entry.entry_version));
    }
    if entry.protein_existence > 0 {
        let _ = write!(out, " \\PE={}", entry.protein_existence);
    }
    if !entry.db_unique_id.is_empty() {
        let _ = write!(out, " \\DbUniqueId={}", escape_peff(&entry.db_unique_id));
    }
    if !entry.entry_id.is_empty() {
        let _ = write!(out, " \\ID={}", escape_peff(&entry.entry_id));
    }
    if !entry.alt_accessions.is_empty() {
        out.push_str(" \\AltAC=");
        for accession in &entry.alt_accessions {
            let _ = write!(out, "({})", escape_peff(accession));
        }
    }
    if !entry.modifications.is_empty() {
        let of = |kind| {
            entry
                .modifications
                .iter()
                .filter(|m| m.modification_type == kind)
                .collect::<Vec<_>>()
        };
        modifications(&mut out, &of(PEFFModificationType::PsiMod), "ModResPsi");
        modifications(&mut out, &of(PEFFModificationType::Unimod), "ModResUnimod");
        modifications(&mut out, &of(PEFFModificationType::Generic), "ModRes");
    }
    if !entry.simple_variants.is_empty() {
        out.push_str(" \\VariantSimple=");
        for variant in &entry.simple_variants {
            out.push('(');
            annotation_id(&mut out, variant.annotation_id);
            let _ = write!(out, "{}|{}", variant.position, variant.variant_aa);
            if !variant.optional_tag.is_empty() {
                let _ = write!(out, "|{}", escape_peff(&variant.optional_tag));
            }
            out.push(')');
        }
    }
    if !entry.complex_variants.is_empty() {
        out.push_str(" \\VariantComplex=");
        for variant in &entry.complex_variants {
            out.push('(');
            annotation_id(&mut out, variant.annotation_id);
            let _ = write!(
                out,
                "{}|{}|{}",
                variant.start_position,
                variant.end_position,
                escape_peff(&variant.replacement)
            );
            if !variant.optional_tag.is_empty() {
                let _ = write!(out, "|{}", escape_peff(&variant.optional_tag));
            }
            out.push(')');
        }
    }
    if !entry.processed_regions.is_empty() {
        out.push_str(" \\Processed=");
        for region in &entry.processed_regions {
            out.push('(');
            annotation_id(&mut out, region.annotation_id);
            let _ = write!(
                out,
                "{}|{}|{}",
                region.start_position,
                region.end_position,
                escape_peff(&region.accession)
            );
            if !region.name.is_empty() || !region.optional_tag.is_empty() {
                let _ = write!(out, "|{}", escape_peff(&region.name));
            }
            if !region.optional_tag.is_empty() {
                let _ = write!(out, "|{}", escape_peff(&region.optional_tag));
            }
            out.push(')');
        }
    }
    if !entry.disulfide_bonds.is_empty() {
        out.push_str(" \\DisulfideBond=");
        for bond in &entry.disulfide_bonds {
            out.push('(');
            annotation_id(&mut out, bond.annotation_id);
            let _ = write!(out, "{},{}", escape_peff(&bond.id1), escape_peff(&bond.id2));
            if !bond.optional_tag.is_empty() {
                let _ = write!(out, "|{}", escape_peff(&bond.optional_tag));
            }
            out.push(')');
        }
    }
    if !entry.proteoforms.is_empty() {
        out.push_str(" \\Proteoform=");
        for proteoform in &entry.proteoforms {
            let _ = write!(out, "({})", escape_peff(proteoform));
        }
    }
    for (key, value) in &entry.custom_annotations {
        let _ = write!(out, " \\{key}={value}");
    }
    out.push('\n');
    // The source writes 80-byte chunks; a sequence is ASCII in practice, and
    // chunks here never split a character, so a non-ASCII sequence gets
    // lines of 80 characters where the source would cut bytes.
    let mut count = 0;
    for c in entry.sequence.chars() {
        out.push(c);
        count += 1;
        if count == 80 {
            out.push('\n');
            count = 0;
        }
    }
    if count > 0 {
        out.push('\n');
    }
    out
}

/// Source `toProForma`.
pub(super) fn to_proforma(entry: &PEFFEntry) -> String {
    if entry.modifications.is_empty() {
        return entry.sequence.clone();
    }
    let format = |modification: &PEFFModification| -> String {
        if !modification.accession.is_empty() {
            format!("[{}]", modification.accession)
        } else if !modification.name.is_empty() {
            format!("[{}]", modification.name)
        } else {
            String::new()
        }
    };
    let mut unlocalised = Vec::new();
    let mut by_position: BTreeMap<usize, Vec<&PEFFModification>> = BTreeMap::new();
    for modification in &entry.modifications {
        if modification.position == 0 {
            unlocalised.push(modification);
        } else {
            by_position
                .entry(modification.position)
                .or_default()
                .push(modification);
        }
    }
    let mut out = String::new();
    if !unlocalised.is_empty() {
        out.push_str("<?");
        for modification in unlocalised {
            out.push_str(&format(modification));
        }
        out.push('>');
    }
    for (index, residue) in entry.sequence.chars().enumerate() {
        out.push(residue);
        if let Some(modifications) = by_position.get(&(index + 1)) {
            for modification in modifications {
                out.push_str(&format(modification));
            }
        }
    }
    out
}

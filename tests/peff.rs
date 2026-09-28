// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/PEFFFile.h`: reading, writing and the record types.
//!
//! Every expectation is transcribed from `PEFFFile_test.cpp` at the pinned
//! core revision, over its eleven unchanged `PEFFFile_*.peff` fixtures (and
//! `FASTAFile_test.fasta`, here `tests/data/fasta_source.fasta`). Section
//! names are quoted in each test's comment. The sequence-producing
//! `PEFFEntry` members are tested in `tests/peff_sequences.rs`. See
//! `tests/data/peff_provenance.json`.

use openms::Error;
use openms::chemistry::proforma::{Peptidoform, SequenceSection};
use openms::format::fasta::FASTAEntry;
use openms::format::peff::{
    self, ANNOTATION_ID_NOT_SET, PEFFDatabaseMetadata, PEFFDisulfideBond, PEFFEntry, PEFFFile,
    PEFFModification, PEFFModificationType, PEFFProcessedRegion, PEFFReader, PEFFSequenceType,
    PEFFVariantComplex, PEFFVariantSimple, ReadLimits,
};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// A fresh directory per test, so parallel tests never share a file.
fn scratch(test: &str) -> PathBuf {
    let directory =
        std::env::temp_dir().join(format!("openms-peff-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir_all(&directory).unwrap();
    directory
}

fn load(path: impl AsRef<Path>) -> (Vec<PEFFEntry>, Vec<PEFFDatabaseMetadata>) {
    let document = PEFFFile::new().load(path).unwrap();
    (document.entries, document.headers)
}

/// A PEFF text written to a file, as the `[EXTRA]` sections do.
fn load_text(test: &str, text: &str) -> Vec<PEFFEntry> {
    let path = scratch(test).join("input.peff");
    std::fs::write(&path, text).unwrap();
    load(&path).0
}

const HEADER: &str = "# PEFF 1.0\n# //\n# DbName=Test\n# Prefix=sp\n# Decoy=false\n\
                      # DbSource=local\n# DbVersion=test\n# NumberOfEntries=1\n\
                      # SequenceType=AA\n# //\n";

/// Section `load`.
#[test]
fn load_reads_headers_and_entries() {
    let mut file = PEFFFile::new();
    assert!(matches!(
        file.load("PEFFFile_test_this_file_does_not_exist.peff"),
        Err(Error::Io(_))
    ));
    let (entries, headers) = load(data("PEFFFile_test.peff"));

    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0].db_name, "TestDB");
    assert_eq!(headers[0].prefix, "test");
    assert_eq!(headers[0].db_description, "Test database for PEFF file parsing");
    assert!(!headers[0].is_decoy);
    assert_eq!(headers[0].db_sources.len(), 1);
    assert_eq!(headers[0].db_version, "1.0");
    assert_eq!(headers[0].number_of_entries, 3);
    assert_eq!(headers[0].sequence_type, PEFFSequenceType::Aa);
    assert_eq!(headers[0].general_comments.len(), 1);
    assert_eq!(headers[0].general_comments[0], "This is a test file");

    assert_eq!(entries.len(), 3);
    let e = &entries[0];
    assert_eq!(e.identifier, "test:P12345");
    assert_eq!(e.sequence, "MSEQVENCEOFTHEFIRSTPROTEINWITHSOMEMODIFICATIONSOK");
    assert_eq!(e.protein_names.len(), 1);
    assert_eq!(e.protein_names[0], "Test Protein 1");
    assert_eq!(e.gene_name, "TEST1");
    assert_eq!(e.ncbi_tax_id, 9606);
    assert_eq!(e.taxonomy_name, "Homo sapiens");
    assert_eq!(e.sequence_length, 49);
    assert_eq!(e.sequence_version, "1");
    assert_eq!(e.entry_version, "1");
    assert_eq!(e.protein_existence, 1);

    assert_eq!(e.modifications.len(), 2);
    assert_eq!(e.modifications[0].position, 10);
    assert_eq!(e.modifications[0].accession, "UNIMOD:35");
    assert_eq!(e.modifications[0].name, "Oxidation");
    assert_eq!(e.modifications[0].modification_type, PEFFModificationType::Unimod);
    assert_eq!(e.modifications[1].position, 25);
    assert_eq!(e.modifications[1].accession, "MOD:00046");
    assert_eq!(e.modifications[1].name, "Phosphorylation");
    assert_eq!(e.modifications[1].modification_type, PEFFModificationType::PsiMod);

    assert_eq!(e.simple_variants.len(), 1);
    assert_eq!(e.simple_variants[0].position, 5);
    assert_eq!(e.simple_variants[0].variant_aa, 'R');
    assert_eq!(e.simple_variants[0].optional_tag, "dbSNP:rs123456");

    assert_eq!(e.processed_regions.len(), 1);
    assert_eq!(e.processed_regions[0].start_position, 1);
    assert_eq!(e.processed_regions[0].end_position, 20);
    assert_eq!(e.processed_regions[0].accession, "PEFF:0001021");
    assert_eq!(e.processed_regions[0].name, "signal peptide");

    let e = &entries[1];
    assert_eq!(e.identifier, "test:P12346");
    assert_eq!(e.protein_names, ["Test Protein 2", "Alternative Name"]);
    assert_eq!(e.complex_variants.len(), 1);
    assert_eq!(e.complex_variants[0].start_position, 5);
    assert_eq!(e.complex_variants[0].end_position, 10);
    assert_eq!(e.complex_variants[0].replacement, "LONGER");

    let e = &entries[2];
    assert_eq!(e.identifier, "test:P12347");
    assert_eq!(e.protein_names, ["Simple Protein"]);
}

/// Section `store`: store, reload and compare; an uncreatable path fails.
#[test]
fn store_round_trips_the_test_file() {
    let directory = scratch("store");
    let (entries, headers) = load(data("PEFFFile_test.peff"));
    let mut file = PEFFFile::new();
    assert!(matches!(
        file.store("/bla/bluff/blblb/sdfhsdjf/test.peff", &entries, &headers[0]),
        Err(Error::Io(_))
    ));
    let path = directory.join("store.peff");
    file.store(&path, &entries, &headers[0]).unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries.len(), entries2.len());
    for (a, b) in entries.iter().zip(&entries2) {
        assert_eq!(a.identifier, b.identifier);
        assert_eq!(a.sequence, b.sequence);
        assert_eq!(a.protein_names.len(), b.protein_names.len());
        assert_eq!(a.modifications.len(), b.modifications.len());
        assert_eq!(a.simple_variants.len(), b.simple_variants.len());
        assert_eq!(a.complex_variants.len(), b.complex_variants.len());
        assert_eq!(a.processed_regions.len(), b.processed_regions.len());
    }
}

/// Sections `readStart`, `readNext`, `getHeaders`.
#[test]
fn streaming_reads_one_entry_at_a_time() {
    let mut file = PEFFFile::new();
    assert!(matches!(file.read_start("nonexistent.peff"), Err(Error::Io(_))));

    file.read_start(data("PEFFFile_test.peff")).unwrap();
    assert_eq!(file.headers().len(), 1);
    assert_eq!(file.headers()[0].db_name, "TestDB");
    assert_eq!(file.headers()[0].version, "1.0");

    let mut entry = PEFFEntry::default();
    assert!(file.read_next(&mut entry).unwrap());
    assert_eq!(entry.identifier, "test:P12345");
    assert!(file.read_next(&mut entry).unwrap());
    assert_eq!(entry.identifier, "test:P12346");
    assert!(file.read_next(&mut entry).unwrap());
    assert_eq!(entry.identifier, "test:P12347");
    assert!(!file.read_next(&mut entry).unwrap());
}

/// Section `atEnd`.
#[test]
fn at_end_after_the_last_entry() {
    let mut file = PEFFFile::new();
    let mut entry = PEFFEntry::default();
    file.read_start(data("PEFFFile_minimal.peff")).unwrap();
    assert!(!file.at_end().unwrap());
    file.read_next(&mut entry).unwrap();
    assert!(file.at_end().unwrap());
}

/// Sections `writeStart`, `writeNext`, `writeEnd`.
#[test]
fn streaming_writes() {
    let directory = scratch("write");
    let mut file = PEFFFile::new();
    let path = directory.join("empty.peff");
    file.write_start(&path, &PEFFDatabaseMetadata::default())
        .unwrap();
    file.write_end().unwrap();
    assert!(path.exists());

    let header = PEFFDatabaseMetadata {
        db_name: "WriteTest".into(),
        ..PEFFDatabaseMetadata::default()
    };
    let entry = PEFFEntry {
        identifier: "test:WRITE001".into(),
        sequence: "TESTSEQUENCE".into(),
        protein_names: vec!["Write Test Protein".into()],
        ..PEFFEntry::default()
    };
    let path = directory.join("one.peff");
    file.write_start(&path, &header).unwrap();
    file.write_next(&entry).unwrap();
    file.write_end().unwrap();
    let (entries, _) = load(&path);
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].identifier, "test:WRITE001");
    assert_eq!(entries[0].sequence, "TESTSEQUENCE");
}

/// Section `isPEFFFile`.
#[test]
fn is_peff_file() {
    assert!(PEFFFile::is_peff_file(data("PEFFFile_test.peff")));
    assert!(!PEFFFile::is_peff_file(data("fasta_source.fasta")));
    assert!(!PEFFFile::is_peff_file("nonexistent.peff"));
}

/// Section `toProForma`, tests 1-6. Test 7 parses the result with the
/// ProForma parser; see `to_proforma_output_parses`.
#[test]
fn to_proforma() {
    let entry = |modifications: Vec<PEFFModification>| PEFFEntry {
        sequence: "PEPTIDE".into(),
        modifications,
        ..PEFFEntry::default()
    };
    assert_eq!(PEFFFile::to_proforma(&entry(vec![])), "PEPTIDE");
    assert_eq!(
        PEFFFile::to_proforma(&entry(vec![PEFFModification::new(3, "UNIMOD:35", "Oxidation")])),
        "PEP[UNIMOD:35]TIDE"
    );
    assert_eq!(
        PEFFFile::to_proforma(&entry(vec![
            PEFFModification::new(1, "UNIMOD:1", "Acetyl"),
            PEFFModification::new(5, "MOD:00046", "Phosphorylation"),
        ])),
        "P[UNIMOD:1]EPTI[MOD:00046]DE"
    );
    assert_eq!(
        PEFFFile::to_proforma(&entry(vec![PEFFModification::new(0, "UNIMOD:35", "Oxidation")])),
        "<?[UNIMOD:35]>PEPTIDE"
    );
    assert_eq!(
        PEFFFile::to_proforma(&entry(vec![
            PEFFModification::new(0, "UNIMOD:35", "Oxidation"),
            PEFFModification::new(3, "MOD:00046", "Phosphorylation"),
        ])),
        "<?[UNIMOD:35]>PEP[MOD:00046]TIDE"
    );
    let named = PEFFModification {
        position: 4,
        name: "Carbamidomethyl".into(),
        ..PEFFModification::default()
    };
    assert_eq!(
        PEFFFile::to_proforma(&entry(vec![named])),
        "PEPT[Carbamidomethyl]IDE"
    );
}

/// Section `toProForma`, test 7: the output parses as ProForma with seven
/// residues and one modification on residues 3 and 5.
#[test]
fn to_proforma_output_parses() {
    let entry = PEFFEntry {
        sequence: "PEPTIDE".into(),
        modifications: vec![
            PEFFModification::new(3, "UNIMOD:35", "Oxidation"),
            PEFFModification::new(5, "MOD:00046", "Phosphorylation"),
        ],
        ..PEFFEntry::default()
    };
    let proforma = PEFFFile::to_proforma(&entry);
    let parsed = Peptidoform::parse(&proforma).unwrap();
    assert_eq!(parsed.sequence.len(), 7);
    let modifications = |index: usize| match &parsed.sequence[index] {
        SequenceSection::Element(element) => element.modifications.len(),
        other => panic!("residue {index} is {other:?}"),
    };
    assert_eq!(modifications(2), 1);
    assert_eq!(modifications(4), 1);
}

/// Sections `[PEFFEntry] toFASTAEntry` and `fromFASTAEntry`.
#[test]
fn fasta_conversion() {
    let entry = PEFFEntry {
        identifier: "test:P12345".into(),
        sequence: "SEQUENCE".into(),
        protein_names: vec!["Test Protein".into()],
        ..PEFFEntry::default()
    };
    let fasta = entry.to_fasta_entry();
    assert_eq!(fasta.identifier, "test:P12345");
    assert_eq!(fasta.sequence, "SEQUENCE");
    assert_eq!(fasta.description, "Test Protein");

    let fasta = FASTAEntry::new("test:P12345", "Test Protein", "SEQUENCE");
    let entry = PEFFEntry::from_fasta_entry(&fasta);
    assert_eq!(entry.identifier, "test:P12345");
    assert_eq!(entry.sequence, "SEQUENCE");
    assert_eq!(entry.protein_names, ["Test Protein"]);
    assert_eq!(entry.sequence_length, 8);
}

/// Sections `[PEFFModification]`, `[PEFFVariantSimple]`,
/// `[PEFFVariantComplex]`, `[PEFFProcessedRegion]`, `[PEFFDisulfideBond]`
/// and `[PEFFDatabaseMetadata]`: default and argument construction.
#[test]
fn record_construction() {
    let modification = PEFFModification::default();
    assert_eq!(modification.position, 0);
    assert!(modification.accession.is_empty());
    assert!(modification.name.is_empty());
    assert!(modification.optional_tag.is_empty());
    assert_eq!(modification.modification_type, PEFFModificationType::Generic);

    let modification =
        PEFFModification::with_details(10, "UNIMOD:35", "Oxidation", "experimental", ANNOTATION_ID_NOT_SET);
    assert_eq!(modification.position, 10);
    assert_eq!(modification.accession, "UNIMOD:35");
    assert_eq!(modification.name, "Oxidation");
    assert_eq!(modification.optional_tag, "experimental");
    assert_eq!(modification.modification_type, PEFFModificationType::Unimod);
    assert_eq!(
        PEFFModification::new(5, "MOD:00046", "Phospho").modification_type,
        PEFFModificationType::PsiMod
    );

    let simple = PEFFVariantSimple::default();
    assert_eq!(simple.position, 0);
    assert_eq!(simple.variant_aa, '\0');
    assert!(simple.optional_tag.is_empty());

    let complex = PEFFVariantComplex::default();
    assert_eq!((complex.start_position, complex.end_position), (0, 0));
    assert!(complex.replacement.is_empty() && complex.optional_tag.is_empty());

    let region = PEFFProcessedRegion::default();
    assert_eq!((region.start_position, region.end_position), (0, 0));
    assert!(region.accession.is_empty() && region.name.is_empty() && region.optional_tag.is_empty());

    let bond = PEFFDisulfideBond::default();
    assert!(bond.id1.is_empty() && bond.id2.is_empty() && bond.optional_tag.is_empty());

    let metadata = PEFFDatabaseMetadata::default();
    assert_eq!(metadata.version, "1.0");
    assert!(metadata.db_name.is_empty());
    assert!(!metadata.is_decoy);
    assert_eq!(metadata.sequence_type, PEFFSequenceType::Aa);
    assert!(!metadata.has_annotation_identifiers);
}

/// Section `test_annotation_identifiers`.
#[test]
fn annotation_identifiers() {
    let (entries, headers) = load(data("PEFFFile_annotid.peff"));
    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0].db_name, "AnnotIDTest");
    assert!(headers[0].has_annotation_identifiers);
    assert_eq!(entries.len(), 2);

    let e = &entries[0];
    assert_eq!(e.identifier, "test:P12345");
    assert_eq!(e.modifications.len(), 2);
    assert_eq!(e.modifications[0].annotation_id, 0);
    assert_eq!(e.modifications[0].position, 10);
    assert_eq!(e.modifications[0].accession, "MOD:00046");
    assert_eq!(e.modifications[1].annotation_id, 1);
    assert_eq!(e.modifications[1].position, 25);
    assert_eq!(e.modifications[1].accession, "UNIMOD:35");
    assert_eq!(e.simple_variants.len(), 1);
    assert_eq!(e.simple_variants[0].annotation_id, 2);
    assert_eq!(e.simple_variants[0].position, 5);
    assert_eq!(e.simple_variants[0].variant_aa, 'R');
    assert_eq!(e.processed_regions.len(), 1);
    assert_eq!(e.processed_regions[0].annotation_id, 3);
    assert_eq!(e.processed_regions[0].start_position, 1);
    assert_eq!(e.processed_regions[0].end_position, 20);

    let e = &entries[1];
    assert_eq!(e.identifier, "test:P12346");
    assert_eq!(e.complex_variants.len(), 1);
    assert_eq!(e.complex_variants[0].annotation_id, 4);
    assert_eq!(e.complex_variants[0].start_position, 5);
    assert_eq!(e.complex_variants[0].end_position, 10);
    assert_eq!(e.complex_variants[0].replacement, "LONGER");
}

/// Section `test_annotation_identifiers_roundtrip`.
#[test]
fn annotation_identifiers_round_trip() {
    let (entries, headers) = load(data("PEFFFile_annotid.peff"));
    let path = scratch("annotid").join("annotid.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (entries2, headers2) = load(&path);
    assert!(headers2[0].has_annotation_identifiers);
    assert_eq!(entries.len(), entries2.len());
    for (a, b) in entries.iter().zip(&entries2) {
        assert_eq!(a.identifier, b.identifier);
        assert_eq!(a.modifications.len(), b.modifications.len());
        for (x, y) in a.modifications.iter().zip(&b.modifications) {
            assert_eq!(x.annotation_id, y.annotation_id);
            assert_eq!(x.position, y.position);
            assert_eq!(x.accession, y.accession);
        }
        assert_eq!(a.simple_variants.len(), b.simple_variants.len());
        for (x, y) in a.simple_variants.iter().zip(&b.simple_variants) {
            assert_eq!(x.annotation_id, y.annotation_id);
            assert_eq!(x.position, y.position);
        }
        assert_eq!(a.complex_variants.len(), b.complex_variants.len());
        for (x, y) in a.complex_variants.iter().zip(&b.complex_variants) {
            assert_eq!(x.annotation_id, y.annotation_id);
            assert_eq!(x.start_position, y.start_position);
        }
        assert_eq!(a.processed_regions.len(), b.processed_regions.len());
        for (x, y) in a.processed_regions.iter().zip(&b.processed_regions) {
            assert_eq!(x.annotation_id, y.annotation_id);
            assert_eq!(x.start_position, y.start_position);
        }
    }
}

/// Section `test_mixed_annotation_formats`.
#[test]
fn annotation_identifiers_unset_without_prefix() {
    let (entries, headers) = load(data("PEFFFile_test.peff"));
    assert!(!headers[0].has_annotation_identifiers);
    assert_eq!(entries[0].modifications[0].annotation_id, u32::MAX);
    assert_eq!(entries[0].modifications[1].annotation_id, u32::MAX);
    assert_eq!(entries[0].simple_variants[0].annotation_id, u32::MAX);
    assert_eq!(entries[0].processed_regions[0].annotation_id, u32::MAX);
}

/// Section `test_minimal_file`.
#[test]
fn minimal_file() {
    let (entries, _) = load(data("PEFFFile_minimal.peff"));
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].identifier, "sp:P00001");
    assert_eq!(entries[0].protein_names, ["Minimal Test"]);
    assert_eq!(entries[0].sequence, "SEQUENCE");
}

/// Section `test_uniprot_format`.
#[test]
fn uniprot_format() {
    let (entries, headers) = load(data("PEFFFile_uniprot.peff"));
    assert_eq!(headers.len(), 1);
    assert_eq!(
        headers[0].general_comments,
        ["this is a hand crafted peff file", "This is a second comment line."]
    );
    assert_eq!(headers[0].conversion, "manual extraction of entries");
    assert_eq!(headers[0].specific_keys.len(), 1);
    assert_eq!(
        headers[0].specific_keys["isoform"],
        "description of a specific isoform"
    );
    assert_eq!(headers[0].specific_values.len(), 1);
    assert_eq!(headers[0].specific_values["isoform"], "(xsd:type=string)");

    assert_eq!(entries.len(), 4);
    let e = &entries[0];
    assert_eq!(e.identifier, "sp:P06748");
    assert_eq!(e.entry_id, "NPM_HUMAN");
    assert_eq!(e.db_unique_id, "P06748");
    assert_eq!(e.alt_accessions, ["Q5SRR7"]);
    assert_eq!(e.gene_name, "NPM1");
    assert_eq!(e.modifications.len(), 2);

    let e = &entries[1];
    assert_eq!(e.identifier, "sp:P02144_CHAIN0");
    assert_eq!(e.entry_id, "MYG_HUMAN");
    assert_eq!(e.ncbi_tax_id, 9606);
    assert_eq!(e.simple_variants.len(), 4);

    let e = &entries[2];
    assert_eq!(e.identifier, "sp:P00761");
    assert_eq!(e.processed_regions.len(), 2);
    assert_eq!(e.processed_regions[0].accession, "PEFF:0001021");
    assert_eq!(e.processed_regions[1].accession, "PEFF:0001020");
    assert_eq!(e.simple_variants.len(), 2);
    assert_eq!(e.simple_variants[1].variant_aa, '*');
    assert_eq!(e.simple_variants[1].optional_tag, "[ESP][ExAC]");

    let e = &entries[3];
    assert_eq!(e.identifier, "sp:NX_P01308-1");
    assert_eq!(e.disulfide_bonds.len(), 3);
    assert_eq!(e.disulfide_bonds[0].id1, "1");
    assert_eq!(e.disulfide_bonds[0].id2, "2");
    assert_eq!(e.disulfide_bonds[0].optional_tag, "between chains");
    assert_eq!(e.disulfide_bonds[2].id1, "5");
    assert_eq!(e.disulfide_bonds[2].id2, "6");
    assert_eq!(e.disulfide_bonds[2].optional_tag, "A chain only");
}

/// Section `test_uniprot_roundtrip`.
#[test]
fn uniprot_round_trip() {
    let (entries, headers) = load(data("PEFFFile_uniprot.peff"));
    let path = scratch("uniprot").join("uniprot.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (entries2, headers2) = load(&path);
    assert_eq!(headers2[0].general_comments.len(), 2);
    assert_eq!(headers2[0].conversion, headers[0].conversion);
    assert_eq!(headers2[0].specific_keys.len(), headers[0].specific_keys.len());
    assert_eq!(headers2[0].specific_values.len(), headers[0].specific_values.len());
    assert_eq!(entries.len(), entries2.len());
    for (a, b) in entries.iter().zip(&entries2) {
        assert_eq!(a.identifier, b.identifier);
        assert_eq!(a.entry_id, b.entry_id);
        assert_eq!(a.db_unique_id, b.db_unique_id);
        assert_eq!(a.alt_accessions.len(), b.alt_accessions.len());
        assert_eq!(a.disulfide_bonds.len(), b.disulfide_bonds.len());
    }
}

/// Section `test_custom_key_def_parsing_and_roundtrip`.
#[test]
fn custom_key_definitions() {
    let (entries, headers) = load(data("PEFFFile_customkeydef.peff"));
    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0].custom_key_defs.len(), 2);
    let definition = &headers[0].custom_key_defs[0];
    assert_eq!(definition.key_name, "SecondaryStructure");
    assert_eq!(definition.description, "Secondary structure information");
    assert_eq!(definition.concept_curie, "BAO:0000014");
    assert_eq!(
        definition.regexp,
        "([0-9]+)\\|([0-9]+)\\|([A-Za-z]+:[0-9]+)?\\|(.*) \\|?(.+)?"
    );
    assert_eq!(definition.field_names.len(), 5);
    assert_eq!(definition.field_names[0], "StartPosition");
    assert_eq!(definition.field_names[4], "OptionalTag");
    assert_eq!(definition.field_types.len(), 5);
    assert_eq!(definition.field_types[0], "integer");
    assert_eq!(definition.field_types[4], "string");
    let definition = &headers[0].custom_key_defs[1];
    assert_eq!(definition.key_name, "EnumTest");
    assert_eq!(definition.description, "Enumeration test");
    assert_eq!(definition.field_names, ["Kind"]);
    assert_eq!(definition.field_types, ["enumeration(alpha|beta|gamma)"]);

    let path = scratch("customkeydef").join("custom.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (_, headers2) = load(&path);
    assert_eq!(headers2.len(), 1);
    assert_eq!(headers2[0].custom_key_defs, headers[0].custom_key_defs);
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.lines().any(|line| line.starts_with("# CustomKeyDef=(")));
}

/// Section `test_custom_key_def_usage_roundtrip`.
#[test]
fn custom_key_usage_round_trip() {
    let (entries, headers) = load(data("PEFFFile_customkeydef_usage.peff"));
    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0].custom_key_defs.len(), 2);
    assert_eq!(entries.len(), 1);
    assert_eq!(
        entries[0].custom_annotations["SecondaryStructure"],
        "(1|2|ncithesaurus:C47937|Helix|Abcg2\\|meta\\\\x10)"
    );
    assert_eq!(entries[0].custom_annotations["EnumTest"], "alpha");

    let path = scratch("customusage").join("usage.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries2.len(), 1);
    assert_eq!(
        entries2[0].custom_annotations["SecondaryStructure"],
        entries[0].custom_annotations["SecondaryStructure"]
    );
    assert_eq!(
        entries2[0].custom_annotations["EnumTest"],
        entries[0].custom_annotations["EnumTest"]
    );
    let text = std::fs::read_to_string(&path).unwrap();
    let line = text
        .lines()
        .find(|line| line.contains("\\SecondaryStructure="))
        .unwrap();
    assert!(line.contains("\\SecondaryStructure=(1|2|ncithesaurus:C47937|Helix|Abcg2\\|meta\\\\x10)"));
}

/// Section `test_p53_uniprot`.
#[test]
fn p53_uniprot() {
    let (entries, headers) = load(data("PEFFFile_P53_uniprot.peff"));
    assert_eq!(headers.len(), 1);
    assert_eq!(headers[0].db_name, "UniProt");
    assert_eq!(headers[0].general_comments.len(), 2);
    assert_eq!(entries.len(), 1);
    let e = &entries[0];
    assert_eq!(e.identifier, "up:P04637");
    assert_eq!(e.gene_name, "TP53");
    assert_eq!(e.ncbi_tax_id, 9606);
    assert_eq!(e.taxonomy_name, "Homo sapiens");
    assert_eq!(e.sequence_length, 393);
    assert_eq!(e.protein_existence, 1);
    assert_eq!(e.sequence_version, "4");
    assert!(e.simple_variants.len() > 1000);
    assert_eq!(e.simple_variants[0].position, 2);
    assert_eq!(e.simple_variants[0].variant_aa, '*');
    assert!(e.simple_variants[0].optional_tag.contains("ExAC"));
    assert_eq!(e.sequence.len(), 393);
}

/// Section `test_modification_type_roundtrip`.
#[test]
fn modification_types_and_comma_positions() {
    let (entries, headers) = load(data("PEFFFile_modtypes.peff"));
    assert_eq!(entries.len(), 2);
    assert_eq!(entries[0].identifier, "mt:P00001");
    assert_eq!(entries[0].modifications.len(), 5);
    let count = |kind| {
        entries[0]
            .modifications
            .iter()
            .filter(|m| m.modification_type == kind)
            .count()
    };
    assert_eq!(count(PEFFModificationType::PsiMod), 2);
    assert_eq!(count(PEFFModificationType::Unimod), 2);
    assert_eq!(count(PEFFModificationType::Generic), 1);

    assert_eq!(entries[1].identifier, "mt:P00002");
    let positions: Vec<_> = entries[1].modifications.iter().map(|m| m.position).collect();
    assert_eq!(positions, [5, 10, 15]);
    for modification in &entries[1].modifications {
        assert_eq!(modification.modification_type, PEFFModificationType::Unimod);
        assert_eq!(modification.accession, "UNIMOD:21");
        assert_eq!(modification.name, "Phospho");
    }

    let path = scratch("modtypes").join("modtypes.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries2[0].modifications.len(), entries[0].modifications.len());
    for (a, b) in entries[0].modifications.iter().zip(&entries2[0].modifications) {
        assert_eq!(a.modification_type, b.modification_type);
        assert_eq!(a.position, b.position);
        assert_eq!(a.accession, b.accession);
    }
}

/// Section `test_modification_optional_tag_roundtrip`.
#[test]
fn modification_optional_tag_round_trip() {
    let (entries, headers) = load(data("PEFFFile_mod_optionaltag.peff"));
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].modifications.len(), 1);
    assert_eq!(entries[0].modifications[0].optional_tag, "invitro");
    let path = scratch("modtag").join("modtag.peff");
    PEFFFile::new().store(&path, &entries, &headers[0]).unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries2.len(), 1);
    assert_eq!(entries2[0].modifications.len(), 1);
    assert_eq!(entries2[0].modifications[0].optional_tag, "invitro");
    let text = std::fs::read_to_string(&path).unwrap();
    let line = text.lines().find(|l| l.contains("\\ModResUnimod=")).unwrap();
    assert!(line.contains("|Phospho|invitro)"));
}

fn complete_header(name: &str, prefix: &str) -> PEFFDatabaseMetadata {
    PEFFDatabaseMetadata {
        db_name: name.into(),
        prefix: prefix.into(),
        db_version: "1.0".into(),
        db_sources: vec!["test".into()],
        number_of_entries: 1,
        ..PEFFDatabaseMetadata::default()
    }
}

/// Section `test_proteoform_db_header`.
#[test]
fn proteoform_db_header() {
    let header = PEFFDatabaseMetadata {
        is_proteoform_db: true,
        ..complete_header("ProteoformTest", "pf")
    };
    let entry = PEFFEntry {
        identifier: "pf:P00001".into(),
        sequence: "PEPTIDE".into(),
        sequence_length: 7,
        ..PEFFEntry::default()
    };
    let path = scratch("proteoform").join("pf.peff");
    PEFFFile::new().store(&path, &[entry], &header).unwrap();
    let (_, headers2) = load(&path);
    assert_eq!(headers2.len(), 1);
    assert!(headers2[0].is_proteoform_db);
}

/// Section `test_altac_parenthesized_list`.
#[test]
fn alt_accession_list() {
    let entry = PEFFEntry {
        identifier: "sp:P12345".into(),
        sequence: "PEPTIDE".into(),
        sequence_length: 7,
        alt_accessions: vec!["Q11111".into(), "Q22222".into(), "Q33333".into()],
        ..PEFFEntry::default()
    };
    let path = scratch("altac").join("altac.peff");
    PEFFFile::new()
        .store(&path, &[entry], &complete_header("AltACTest", "sp"))
        .unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries2[0].alt_accessions, ["Q11111", "Q22222", "Q33333"]);
}

/// Section `test_disulfide_bond_annotation_id`.
#[test]
fn disulfide_bond_annotation_id() {
    let entry = PEFFEntry {
        identifier: "sp:P12345".into(),
        sequence: "PEPTIDE".into(),
        sequence_length: 7,
        disulfide_bonds: vec![PEFFDisulfideBond::new("1", "2", "between chains", 5)],
        ..PEFFEntry::default()
    };
    let header = PEFFDatabaseMetadata {
        has_annotation_identifiers: true,
        ..complete_header("DSBondTest", "sp")
    };
    let path = scratch("dsbond").join("dsbond.peff");
    PEFFFile::new().store(&path, &[entry], &header).unwrap();
    let (entries2, _) = load(&path);
    assert_eq!(entries2[0].disulfide_bonds.len(), 1);
    assert_eq!(entries2[0].disulfide_bonds[0].annotation_id, 5);
    assert_eq!(entries2[0].disulfide_bonds[0].id1, "1");
    assert_eq!(entries2[0].disulfide_bonds[0].id2, "2");
    assert_eq!(entries2[0].disulfide_bonds[0].optional_tag, "between chains");
}

/// Section `test_multi_database_store_roundtrip`.
#[test]
fn multi_database_store_round_trip() {
    let db1 = complete_header("Database1", "db1");
    let db2 = PEFFDatabaseMetadata {
        is_decoy: true,
        db_version: "2.0".into(),
        ..complete_header("Database2", "db2")
    };
    let entry = |identifier: &str, sequence: &str, name: &str| PEFFEntry {
        identifier: identifier.into(),
        sequence: sequence.into(),
        sequence_length: 10,
        protein_names: vec![name.into()],
        ..PEFFEntry::default()
    };
    let entries = [
        entry("db1:P00001", "SEQUENCEONE", "Protein from DB1"),
        entry("db2:P00002", "SEQUENCETWO", "Protein from DB2"),
    ];
    let path = scratch("multidb").join("multi.peff");
    PEFFFile::new()
        .store_multi(&path, &entries, &[db1, db2])
        .unwrap();
    let (entries2, headers2) = load(&path);
    assert_eq!(headers2.len(), 2);
    assert_eq!(headers2[0].prefix, "db1");
    assert_eq!(headers2[0].db_name, "Database1");
    assert!(!headers2[0].is_decoy);
    assert_eq!(headers2[1].prefix, "db2");
    assert_eq!(headers2[1].db_name, "Database2");
    assert!(headers2[1].is_decoy);
    assert_eq!(entries2.len(), 2);
    assert_eq!(entries2[0].identifier, "db1:P00001");
    assert_eq!(entries2[1].identifier, "db2:P00002");
    let text = std::fs::read_to_string(&path).unwrap();
    assert_eq!(text.lines().filter(|l| l.starts_with("# PEFF ")).count(), 1);
    assert!(matches!(
        PEFFFile::new().store_multi(&path, &entries, &[]),
        Err(Error::InvalidValue(_))
    ));
}

/// Section `[EXTRA] round-trip: PEFFFile must consume the byte-exact output
/// of the UniPEFF TOPP tool`, apart from its two `AASequence` checks, which
/// are in `tests/peff_sequences.rs`.
#[test]
fn unipeff_output_loads() {
    let (entries, headers) = load(data("PEFFFile_UniPEFF_test.peff"));
    assert_eq!(entries.len(), 2);
    assert_eq!(headers.len(), 2);
    assert_eq!(headers[0].prefix, "sp");
    assert_eq!(headers[1].prefix, "tr");
    assert_eq!(headers[0].number_of_entries, 1);
    assert_eq!(headers[1].number_of_entries, 1);
    assert!(!headers[0].has_annotation_identifiers);

    let e1 = &entries[0];
    assert_eq!(e1.identifier, "sp:P12345");
    assert_eq!(e1.prefix, "sp");
    assert_eq!(e1.sequence_length, 100);
    assert_eq!(e1.sequence.len(), 100);
    assert_eq!(e1.ncbi_tax_id, 9606);
    assert_eq!(e1.taxonomy_name, "Homo sapiens");
    assert_eq!(e1.gene_name, "GENE1");
    assert_eq!(e1.entry_id, "KSINK_HUMAN");
    assert_eq!(e1.alt_accessions, ["P99999"]);
    assert_eq!(
        e1.protein_names,
        ["Kitchen-sink protein (test) | alpha\\beta :-)"]
    );
    assert_eq!(e1.modifications.len(), 18);
    assert_eq!(
        e1.modifications
            .iter()
            .filter(|m| m.accession == "MOD:00798")
            .count(),
        5
    );
    assert_eq!(e1.simple_variants.len(), 2);
    assert_eq!(e1.complex_variants.len(), 5);
    assert_eq!(e1.processed_regions.len(), 5);
    assert_eq!(e1.disulfide_bonds.len(), 0);

    let e2 = &entries[1];
    assert_eq!(e2.identifier, "tr:Q67890");
    assert_eq!(e2.prefix, "tr");
    assert_eq!(e2.sequence_length, 31);
}

/// Section `[EXTRA] escaped '\|' inside a structured tuple must not be split
/// as a field separator`.
#[test]
fn escaped_pipe_stays_in_its_field() {
    let entries = load_text(
        "pipe",
        &format!(
            "{HEADER}>sp:P00001 \\PName=(Test) \\Length=5 \\ModResPsi=(1|MOD:00046|name with \\| pipe)\nINSUL\n"
        ),
    );
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].modifications.len(), 1);
    assert_eq!(entries[0].modifications[0].accession, "MOD:00046");
    assert_eq!(entries[0].modifications[0].name, "name with | pipe");
    assert_eq!(entries[0].modifications[0].optional_tag, "");
}

/// Section `[EXTRA] literal backslash in a structured tuple text field must
/// not double-decode`: four backslashes on the wire are two in the name.
#[test]
fn literal_backslashes_decode_once() {
    let entries = load_text(
        "backslash",
        &format!(
            "{HEADER}>sp:P00001 \\PName=(Test) \\Length=5 \\ModResPsi=(1|MOD:00046|name with \\\\\\\\ backslashes)\nINSUL\n"
        ),
    );
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].modifications.len(), 1);
    assert_eq!(entries[0].modifications[0].name, "name with \\\\ backslashes");
}

/// Section `[EXTRA] scalar \PName= with embedded balanced parens must
/// round-trip in full`.
#[test]
fn scalar_pname_with_parentheses() {
    let entries = load_text(
        "pname",
        &format!("{HEADER}>sp:P00001 \\PName=Insulin (Fragment) \\GName=INS \\Length=5\nINSUL\n"),
    );
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].protein_names, ["Insulin (Fragment)"]);
    assert_eq!(entries[0].gene_name, "INS");
}

/// The source's warnings are returned, verbatim, in order: the required-key
/// check of `readStart` (`PEFFFile.cpp:1497-1526`) on the minimal fixture,
/// which has only a version and one entry.
#[test]
fn header_warnings_are_the_source_lines() {
    let document = PEFFFile::new().load(data("PEFFFile_minimal.peff")).unwrap();
    assert_eq!(
        document.warnings,
        [
            "PEFF: Database block should start with 'DbName' key in database block 1.",
            "PEFF: Required key 'Prefix' missing in database block 1.",
            "PEFF: Required key 'DbVersion' missing in database block 1.",
            "PEFF: Required key 'DbSource' missing in database block 1.",
            "PEFF: Required key 'NumberOfEntries' missing or zero in database block 1.",
        ]
    );
}

/// `readNext`'s two `ParseError` messages (`PEFFFile.cpp:1542-1552`), and a
/// `;` line, which the source warns about and then fails on.
#[test]
fn unreadable_entries_are_parse_errors() {
    let error = |text: &str| peff::read(text.as_bytes()).unwrap_err();
    match error("# PEFF 1.0\n>\n") {
        Error::Parse { message, .. } => {
            assert_eq!(message, "The first PEFF entry could not be read!");
        }
        other => panic!("{other}"),
    }
    match error(">a\nSEQ\n>b desc\n") {
        Error::Parse { message, .. } => assert_eq!(
            message,
            "Only 1 entries could be read. Parsing next record failed."
        ),
        other => panic!("{other}"),
    }
    let mut reader = PEFFReader::new(b"; comment\n>a\nSEQ\n".as_slice()).unwrap();
    assert!(matches!(reader.next_entry(), Err(Error::Parse { .. })));
    assert_eq!(
        reader.warnings().last().unwrap(),
        "PEFF: Lines beginning with ';' are not permitted."
    );
}

#[test]
fn limits_refuse_instead_of_truncating() {
    let text = std::fs::read(data("PEFFFile_P53_uniprot.peff")).unwrap();
    let limits = ReadLimits {
        max_items_per_entry: 100,
        ..ReadLimits::default()
    };
    assert!(matches!(
        peff::read_with_limits(text.as_slice(), &limits),
        Err(Error::InvalidValue(_))
    ));
    let limits = ReadLimits {
        max_sequence_bytes: 100,
        ..ReadLimits::default()
    };
    assert!(matches!(
        peff::read_with_limits(text.as_slice(), &limits),
        Err(Error::InvalidValue(_))
    ));
    let limits = ReadLimits {
        max_warnings: 0,
        ..ReadLimits::default()
    };
    let minimal = std::fs::read(data("PEFFFile_minimal.peff")).unwrap();
    let mut reader = PEFFReader::with_limits(minimal.as_slice(), limits).unwrap();
    assert!(reader.warnings().is_empty());
    assert_eq!(reader.suppressed_warnings(), 5);
    assert!(reader.next_entry().unwrap().is_some());
}

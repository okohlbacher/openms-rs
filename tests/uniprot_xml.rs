// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/UniProtXMLFile.h` and `FORMAT/HANDLERS/UniProtXMLHandler.h`.
//!
//! Every expectation in the two class-test sections comes from
//! `UniProtXMLFile_test.cpp`, over its unchanged fixture
//! `UniProtXMLFile_test.xml`. See `tests/data/uniprot_xml_provenance.json`.

#![cfg(feature = "kb-xml")]

use openms::Error;
use openms::format::uniprot_xml::{self, ReadLimits, UniProtEntry};
use std::io::Write;
use std::path::PathBuf;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/UniProtXMLFile_test.xml")
}

/// `UniProtXMLFile_test.cpp`, section `load`.
#[test]
fn load_reads_the_class_test_fixture() {
    let entries = uniprot_xml::load(fixture()).unwrap();
    // Two entries: one Swiss-Prot kitchen-sink, one minimal TrEMBL.
    assert_eq!(entries.len(), 2);

    let e1 = &entries[0];
    assert_eq!(e1.dataset, "Swiss-Prot");
    assert_eq!(e1.entry_version, "42");
    assert_eq!(e1.accession, "P00001");
    assert_eq!(e1.alt_accessions.len(), 1);
    assert_eq!(e1.alt_accessions[0], "P99999");
    assert_eq!(e1.name, "TEST1_HUMAN");
    assert_eq!(e1.full_name, "Test protein one");
    // The primary gene name is captured even though a synonym appears first.
    assert_eq!(e1.primary_gene, "PRIMGENE");
    assert_eq!(e1.ncbi_tax_id, "9606");
    // The scientific organism name is captured even though a common name
    // appears first.
    assert_eq!(e1.tax_name, "Homo sapiens");
    assert_eq!(e1.protein_existence, "evidence at protein level");
    assert_eq!(e1.sequence_version, "2");
    assert_eq!(e1.sequence, "MSEQNCEMTFQIQRIYNCKDISFEAPNAPC");
    // Isoform sequences inside <alternative products> are skipped, so the
    // canonical 30-residue sequence is the one read.
    assert_eq!(e1.sequence.len(), 30);

    // Four features: signal peptide, modified residue, disulfide bond,
    // sequence variant.
    assert_eq!(e1.features.len(), 4);

    let f_sig = &e1.features[0];
    assert_eq!(f_sig.feature_type, "signal peptide");
    assert!(f_sig.has_range);
    assert_eq!(f_sig.begin, 1);
    assert_eq!(f_sig.end, 10);

    let f_mod = &e1.features[1];
    assert_eq!(f_mod.feature_type, "modified residue");
    assert_eq!(f_mod.description, "Phosphoserine");
    assert!(f_mod.has_position);
    assert_eq!(f_mod.position, 20);

    let f_ss = &e1.features[2];
    assert_eq!(f_ss.feature_type, "disulfide bond");
    assert!(f_ss.has_range);
    assert_eq!(f_ss.begin, 15);
    assert_eq!(f_ss.end, 25);

    let f_var = &e1.features[3];
    assert_eq!(f_var.feature_type, "sequence variant");
    assert_eq!(f_var.original, "A");
    assert_eq!(f_var.variation, "V");
    assert!(f_var.has_position);
    assert_eq!(f_var.position, 30);

    let e2 = &entries[1];
    assert_eq!(e2.dataset, "TrEMBL");
    assert_eq!(e2.accession, "Q00002");
    assert_eq!(e2.name, "TEST2_HUMAN");
    assert_eq!(e2.sequence, "MSEQE");
    assert_eq!(e2.features.len(), 0);
}

/// `UniProtXMLFile_test.cpp`, section `loadStreaming`.
#[test]
fn load_streaming_delivers_each_entry() {
    let mut entry_count = 0;
    let mut accessions = Vec::new();
    uniprot_xml::load_streaming(fixture(), |entry: UniProtEntry| {
        entry_count += 1;
        accessions.push(entry.accession);
        Ok(())
    })
    .unwrap();
    assert_eq!(entry_count, 2);
    assert_eq!(accessions.len(), 2);
    assert_eq!(accessions[0], "P00001");
    assert_eq!(accessions[1], "Q00002");
}

/// The header documents `.xml.gz` input, detected by its magic bytes rather
/// than by the file name.
#[test]
fn gzip_input_is_detected_by_content() {
    let bytes = std::fs::read(fixture()).unwrap();
    let path = std::env::temp_dir().join(format!("openms-uniprot-{}.dat", std::process::id()));
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&path).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(&bytes).unwrap();
    encoder.finish().unwrap();
    assert_eq!(
        uniprot_xml::load(&path).unwrap(),
        uniprot_xml::load(fixture()).unwrap()
    );
    std::fs::remove_file(&path).unwrap();
}

#[test]
fn failures_and_bounds() {
    assert!(matches!(
        uniprot_xml::load(fixture().with_extension("missing")),
        Err(Error::Io(_))
    ));
    assert!(matches!(
        uniprot_xml::read(b"<uniprot><entry></uniprot>".as_slice()),
        Err(Error::Parse { .. })
    ));
    let limits = ReadLimits {
        max_features_per_entry: 3,
        ..ReadLimits::default()
    };
    assert!(matches!(
        uniprot_xml::load_with_limits(fixture(), &limits),
        Err(Error::InvalidValue(_))
    ));
    let limits = ReadLimits {
        max_entries: 1,
        ..ReadLimits::default()
    };
    assert!(matches!(
        uniprot_xml::load_with_limits(fixture(), &limits),
        Err(Error::InvalidValue(_))
    ));
}

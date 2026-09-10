// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::file_types::*;

#[test]
fn complete_registry_and_source_property_literal() {
    assert_eq!(FileType::ALL.len(), 73);
    let unique: std::collections::BTreeSet<_> = FileType::ALL.iter().collect();
    assert_eq!(unique.len(), 73);
    for &t in FileType::ALL {
        assert_eq!(FileType::from_name(t.name()), t);
        assert_eq!(FileType::from_name(&t.name().to_ascii_uppercase()), t);
        assert!(!t.description().contains("  "));
    }
    assert_eq!(FileType::from_name("pqt"), FileType::Parquet);
    assert_eq!(FileType::from_name(".mzML"), FileType::Unknown);
    let readable = FileTypeList::with_source_properties(&[FileProperty::Readable]);
    assert_eq!(readable.0.len(), 50); // Pinned FileTypes_test.cpp expectation.
    assert!(readable.contains(FileType::PeakMapParquet));
    assert!(!readable.contains(FileType::Yaml));
    assert_eq!(FileTypeList::with_source_properties(&[]).0, FileType::ALL);
    let both = FileTypeList::with_source_properties(&[
        FileProperty::Readable,
        FileProperty::ProvidesSpectrum,
    ]);
    assert_eq!(both.0, [FileType::Dta, FileType::Xmass]);
    for &t in FileType::ALL {
        assert_eq!(
            t.is_directory(),
            [
                FileType::BrukerTdf,
                FileType::IdParquet,
                FileType::FeatureParquet,
                FileType::ConsensusParquet
            ]
            .contains(&t)
        );
    }
    assert_eq!(FileType::ImzMl.mzml_name(), "mzML file");
    assert_eq!(FileType::Dta2d.mzml_name(), "DTA file");
    assert_eq!(FileType::IdXml.mzml_name(), "");
}

#[test]
fn source_dialog_filter_literals_and_exact_inverse() {
    let list = FileTypeList(vec![FileType::MzMl, FileType::Bz2]);
    assert_eq!(
        list.to_file_dialog_filter(FilterLayout::Both, true),
        "all readable files (*.mzML *.bz2);;mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2);;all files (*)"
    );
    for layout in [
        FilterLayout::Both,
        FilterLayout::Compact,
        FilterLayout::OneByOne,
    ] {
        for all in [true, false] {
            for item in list.to_file_dialog_filter(layout, all).split(";;") {
                let expected = if item.starts_with("all ") {
                    FileType::IdXml
                } else if item.starts_with("mzML") {
                    FileType::MzMl
                } else {
                    FileType::Bz2
                };
                assert_eq!(
                    list.from_file_dialog_filter(item, FileType::IdXml).unwrap(),
                    expected
                );
            }
        }
    }
    assert!(
        list.from_file_dialog_filter("mzML files (*.mzML)", FileType::Unknown)
            .is_err()
    );
    assert_eq!(
        FileTypeList::default().to_file_dialog_filter(FilterLayout::Compact, false),
        "all readable files ()"
    );
}

#[test]
fn source_filename_vectors_and_compression_without_recursion() {
    for (name, kind, stem) in [
        ("", FileType::Unknown, ""),
        ("fid", FileType::Xmass, "fid"),
        ("sample.mzML.gz", FileType::MzMl, "sample"),
        ("sample.mzML.bz2", FileType::MzMl, "sample"),
        ("sample.d.zip", FileType::BrukerTdf, "sample"),
        ("sample.pep.xml", FileType::PepXml, "sample.pep"),
        ("sample.prot.xml", FileType::ProtXml, "sample.prot"),
        ("sample.pep.xml.gz", FileType::PepXml, "sample.pep.xml"),
        ("sample.newEnding", FileType::Unknown, "sample"),
        (
            "/dotted.directory/sample",
            FileType::Unknown,
            "/dotted.directory/sample",
        ),
        (
            r"C:\dotted.directory\sample.mzML.gz",
            FileType::MzMl,
            r"C:\dotted.directory\sample",
        ),
        (
            r"C:\dotted.directory\sample",
            FileType::Unknown,
            r"C:\dotted.directory\sample",
        ),
    ] {
        assert_eq!(type_by_file_name(name), kind, "{name}");
        assert_eq!(strip_extension(name), stem, "{name}");
        assert_eq!(swap_extension(name, FileType::MzMl), format!("{stem}.mzML"));
    }
    assert_eq!(
        type_by_file_name(&format!("x.mzML{}", ".gz".repeat(10_000))),
        FileType::MzMl
    );
    assert_eq!(type_by_file_name("x.PEP.XML"), FileType::Xml); // Source special aliases are case-sensitive.
    assert_eq!(type_by_file_name("x.XQUEST.XML"), FileType::Xml);
    assert_eq!(type_by_file_name("x.pqt"), FileType::Parquet);
    assert!(has_valid_extension("out.unknownSuffix", FileType::MzMl));
    assert!(!has_valid_extension("out.idXML", FileType::MzMl));
    assert_eq!(
        consistent_output_type("out.idXML", "mzML"),
        FileType::Unknown
    );
    assert_eq!(consistent_output_type("out", "MzML"), FileType::MzMl);
    assert_eq!(consistent_output_type("out.mzML", ""), FileType::MzMl);
    assert_eq!(strip_extension("unknown.foo"), "unknown.foo"); // Checked source underflow.
}

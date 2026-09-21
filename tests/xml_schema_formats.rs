// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Each ported `XMLFile`'s `isValid`, through the format's own entry point, and
//! the `isValid` sections of the class tests that validate what the writer
//! stored. Each test needs its format's feature as well as `xml-schema`.
//!
//! A source section such as `FeatureXMLFile_test.cpp:394-403` stores with the
//! C++ writer and asserts that the file validates; here the same content is
//! stored with this port's writer, and the same assertion is made of it. That
//! is the contract the class test states for a writer, not a comparison with
//! C++ output.
#![cfg(feature = "xml-schema")]
#![allow(unused_imports)]

use openms::Error;
use openms::format::xml_schema::{SchemaKind, SchemaValidationReport};
use std::path::{Path, PathBuf};

fn data(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative)
}

fn temp() -> openms::system::file::TempDir {
    openms::system::file::TempDir::new_in(std::env::temp_dir(), false).unwrap()
}

#[track_caller]
fn valid(r: SchemaValidationReport, kind: SchemaKind) {
    assert!(r.is_valid(), "{r:?}");
    assert!(r.diagnostics.is_empty(), "{r:?}");
    assert_eq!(r.schema, kind);
}

// ===========================================================================
// FeatureXMLFile_test.cpp:385-405
// ===========================================================================

#[cfg(feature = "featurexml")]
#[test]
fn featurexml_fixtures_and_stored_maps_are_valid() {
    use openms::format::featurexml;
    use openms::kernel::FeatureMap;
    for fixture in [
        "featurexml_source_1.featureXML",
        "featurexml_source_options.featureXML",
    ] {
        valid(
            featurexml::is_valid(data(fixture)).unwrap(),
            SchemaKind::FeatureXML,
        );
    }
    let directory = temp();
    // "test if empty file is valid"
    let empty = directory.path().join("empty.featureXML");
    featurexml::store(&empty, &FeatureMap::default()).unwrap();
    valid(
        featurexml::is_valid(&empty).unwrap(),
        SchemaKind::FeatureXML,
    );
    // "test if full file is valid"
    let full = directory.path().join("full.featureXML");
    let map = featurexml::load(data("featurexml_source_1.featureXML")).unwrap();
    featurexml::store(&full, &map).unwrap();
    valid(featurexml::is_valid(&full).unwrap(), SchemaKind::FeatureXML);
    assert!(matches!(
        featurexml::is_valid(directory.path().join("missing.featureXML")),
        Err(Error::Io(_))
    ));
}

// ===========================================================================
// ConsensusXMLFile_test.cpp:274-289 and 331-422
// ===========================================================================

#[cfg(feature = "consensusxml")]
#[test]
fn consensusxml_fixtures_and_stored_maps_are_valid() {
    use openms::format::consensusxml;
    for fixture in [
        "consensusxml/ConsensusXMLFile_1.consensusXML",
        "consensusxml/ConsensusXMLFile_2_options.consensusXML",
    ] {
        valid(
            consensusxml::is_valid(data(fixture)).unwrap(),
            SchemaKind::ConsensusXML,
        );
    }
    // "test if written empty file - this is invalid, so it is not tested :)"
    // The source asserts nothing about an empty map, so neither does this.
    let directory = temp();
    let full = directory.path().join("full.consensusXML");
    let map = consensusxml::load(data("consensusxml/ConsensusXMLFile_1.consensusXML")).unwrap();
    consensusxml::store(&full, &map).unwrap();
    valid(
        consensusxml::is_valid(&full).unwrap(),
        SchemaKind::ConsensusXML,
    );
}

/// `[EXTRA] Protein group quantities round-trip`: the group
/// `PeptideAndProteinQuant::annotateQuantificationsToProteins` would attach,
/// stored, and "the file is still schema-valid".
#[cfg(feature = "consensusxml")]
#[test]
fn consensusxml_with_protein_group_quantities_is_valid() {
    use openms::format::consensusxml;
    use openms::identification::ProteinGroup;
    use openms::kernel::DataArray;
    let mut map = consensusxml::load(data("consensusxml/ConsensusXMLFile_1.consensusXML")).unwrap();
    let run = &mut map.protein_identifications[0];
    let accessions: Vec<String> = run.hits.iter().map(|h| h.accession.clone()).collect();
    assert!(!accessions.is_empty());
    run.indistinguishable_groups.push(ProteinGroup {
        probability: 0.75,
        accessions,
        float_data_arrays: vec![
            DataArray::new("psm_count", vec![]),
            DataArray::new("distinct_peptides", vec![]),
            DataArray::new("file_channel_level_abundance", vec![10.0, 20.0, 30.0, 40.0]),
            DataArray::new("fraction_group_level_abundance", vec![1.5, 2.5]),
        ],
        string_data_arrays: vec![
            DataArray::new(
                "file_channel_level_filename",
                vec![
                    "fileA".into(),
                    "fileA".into(),
                    "fileB".into(),
                    "fileB".into(),
                ],
            ),
            DataArray::new("file_level_filename", vec![]),
        ],
        integer_data_arrays: vec![
            DataArray::new("file_channel_level_channel", vec![1, 2, 1, 2]),
            DataArray::new("file_level_psm_count", vec![]),
            DataArray::new("fraction_group_level_fraction_group", vec![1, 2]),
            DataArray::new("fraction_group_level_label", vec![1, 1]),
        ],
    });
    let directory = temp();
    let path = directory.path().join("quantities.consensusXML");
    consensusxml::store(&path, &map).unwrap();
    let loaded = consensusxml::load(&path).unwrap();
    assert_eq!(
        loaded.protein_identifications[0]
            .indistinguishable_groups
            .len(),
        1
    );
    valid(
        consensusxml::is_valid(&path).unwrap(),
        SchemaKind::ConsensusXML,
    );
}

// ===========================================================================
// IdXMLFile_test.cpp:268-291
// ===========================================================================

#[cfg(feature = "idxml")]
#[test]
fn idxml_stored_documents_are_valid() {
    use openms::format::idxml;
    use openms::metadata::MetaValue;
    let directory = temp();
    // "test if empty file is valid" (IdXMLFile_test.cpp:274-277) is not
    // ported: for no protein run the source writes a placeholder
    // SearchParameters and IdentificationRun (IdXMLFile.cpp:186-193 and
    // 416-418), and this port's writer refuses the document instead
    // (src/format/idxml.rs, write_with_registry), so there is no stored file
    // to validate. That writer gap is recorded in docs/XML_SCHEMA_SUPPORT.md;
    // when the writer gains the placeholder, store IdXmlDocument::default()
    // here and assert that it validates.
    // "test if full file is valid", with the three meta values it adds.
    let mut document = idxml::load(data("idxml_upstream_whole.idXML")).unwrap();
    let meta = &mut document.protein_identifications[0].metadata;
    meta.insert("stringvalue".into(), MetaValue::from("bla"));
    meta.insert("intvalue".into(), MetaValue::from(4711));
    meta.insert("floatvalue".into(), MetaValue::try_from(5.3).unwrap());
    let full = directory.path().join("full.idXML");
    idxml::store(&full, &document).unwrap();
    valid(idxml::is_valid(&full).unwrap(), SchemaKind::IdXML);
    // "check if meta information can be loaded"
    let reloaded = idxml::load(&full).unwrap();
    assert_eq!(
        reloaded.protein_identifications[0]
            .metadata
            .get("stringvalue"),
        Some(&MetaValue::from("bla"))
    );
}

// ===========================================================================
// ParamXMLFile_test.cpp:53-61 and 63-243
// ===========================================================================

#[cfg(feature = "paramxml")]
#[test]
fn paramxml_stored_parameters_are_valid() {
    use openms::format::paramxml;
    use openms::param::{Param, ParamValue};
    fn tags(t: &[&str]) -> Vec<String> {
        t.iter().map(|s| (*s).to_owned()).collect()
    }
    let directory = temp();
    let check = |name: &str, param: &Param| {
        let path = directory.path().join(name);
        paramxml::store(&path, param).unwrap();
        valid(paramxml::is_valid(&path).unwrap(), SchemaKind::ParamXML);
    };

    // The test's shared `p`, extended to `p2` in the store section.
    let mut p = Param::new();
    p.set_value("test:float", ParamValue::from(17.4f32), "floatdesc", &[])
        .unwrap();
    p.set_value("test:string", "test,test,test".into(), "stringdesc", &[])
        .unwrap();
    p.set_value("test:int", ParamValue::from(17), "intdesc", &[])
        .unwrap();
    p.set_value("test2:float", ParamValue::from(17.5f32), "", &[])
        .unwrap();
    p.set_value("test2:string", "test2".into(), "", &[])
        .unwrap();
    p.set_value("test2:int", ParamValue::from(18), "", &[])
        .unwrap();
    p.set_section_description("test", "sectiondesc").unwrap();
    p.add_tags("test:float", &tags(&["a", "b", "c"])).unwrap();
    let mut p2 = p.clone();
    p2.set_value(
        "test:a:a1",
        ParamValue::from(47.1),
        "a1desc\"<>\nnewline",
        &[],
    )
    .unwrap();
    p2.set_value("test:b:b1", ParamValue::from(47.1), "", &[])
        .unwrap();
    p2.set_section_description("test:b", "bdesc\"<>\nnewline")
        .unwrap();
    p2.set_value("test2:a:a1", ParamValue::from(47.1), "", &[])
        .unwrap();
    p2.set_value(
        "test2:b:b1",
        ParamValue::from(47.1),
        "",
        &tags(&["advanced"]),
    )
    .unwrap();
    p2.set_section_description("test2:a", "adesc").unwrap();
    check("p2.ini", &p2);

    // "advanced"
    let mut p7 = Param::new();
    p7.set_value("true", ParamValue::from(5), "", &tags(&["advanced"]))
        .unwrap();
    p7.set_value("false", ParamValue::from(5), "", &[]).unwrap();
    check("p7.ini", &p7);

    // "restrictions"
    let mut p5 = Param::new();
    let set =
        |p: &mut Param, key: &str, value: ParamValue| p.set_value(key, value, "", &[]).unwrap();
    set(&mut p5, "int", 5.into());
    set(&mut p5, "int_min", 5.into());
    p5.set_min_int("int_min", 4).unwrap();
    set(&mut p5, "int_max", 5.into());
    p5.set_max_int("int_max", 6).unwrap();
    set(&mut p5, "int_min_max", 5.into());
    p5.set_min_int("int_min_max", 0).unwrap();
    p5.set_max_int("int_min_max", 10).unwrap();
    set(&mut p5, "float", 5.1.into());
    set(&mut p5, "float_min", 5.1.into());
    p5.set_min_float("float_min", 4.1).unwrap();
    set(&mut p5, "float_max", 5.1.into());
    p5.set_max_float("float_max", 6.1).unwrap();
    set(&mut p5, "float_min_max", 5.1.into());
    p5.set_min_float("float_min_max", 0.1).unwrap();
    p5.set_max_float("float_min_max", 10.1).unwrap();
    set(&mut p5, "string", "bli".into());
    set(&mut p5, "string_2", "bla".into());
    p5.set_valid_strings("string_2", &tags(&["bla", "bluff"]))
        .unwrap();
    set(
        &mut p5,
        "stringlist2",
        tags(&["a.txt", "b.xml", "c.pdf"]).into(),
    );
    set(
        &mut p5,
        "stringlist",
        tags(&["aa.C", "bb.h", "c.doxygen"]).into(),
    );
    p5.set_valid_strings("stringlist2", &tags(&["xml", "txt"]))
        .unwrap();
    for key in ["intlist", "intlist2", "intlist3", "intlist4"] {
        set(&mut p5, key, vec![2, 5, 10].into());
    }
    p5.set_min_int("intlist2", 1).unwrap();
    p5.set_max_int("intlist3", 11).unwrap();
    p5.set_min_int("intlist4", 0).unwrap();
    p5.set_max_int("intlist4", 15).unwrap();
    for key in ["doublelist", "doublelist2", "doublelist3", "doublelist4"] {
        set(&mut p5, key, vec![1.2, 3.33, 4.44].into());
    }
    p5.set_min_float("doublelist2", 1.1).unwrap();
    p5.set_max_float("doublelist3", 4.45).unwrap();
    p5.set_min_float("doublelist4", 0.1).unwrap();
    p5.set_max_float("doublelist4", 5.8).unwrap();
    check("p5.ini", &p5);

    // "Test if an empty Param written to a file validates against the schema"
    check("p4.ini", &Param::new());
}

// ===========================================================================
// TransformationXMLFile_test.cpp:38-45
// ===========================================================================

#[cfg(any(feature = "featurexml", feature = "consensusxml"))]
#[test]
fn transformation_xml_fixtures_validate_as_the_class_test_asserts() {
    use openms::format::transformation_xml;
    for fixture in [
        "transformation_xml_1.trafoXML",
        "transformation_xml_2.trafoXML",
        "transformation_xml_4.trafoXML",
    ] {
        valid(
            transformation_xml::is_valid(data(fixture)).unwrap(),
            SchemaKind::TransformationXML,
        );
    }
    // The source's `false` for a document whose end tags do not match.
    assert!(matches!(
        transformation_xml::is_valid(data("transformation_xml_3.trafoXML")),
        Err(Error::Parse { .. })
    ));
}

// ===========================================================================
// MzDataFile_test.cpp:830-847
// ===========================================================================

#[cfg(feature = "mzml")]
#[test]
fn mzdata_stored_experiments_are_valid() {
    use openms::MSExperiment;
    use openms::format::mzdata::{self, MzDataFile, WriteOptions};
    let file = MzDataFile::new();
    let directory = temp();
    // "test if empty file is valid"
    let empty = directory.path().join("empty.mzData");
    file.store(&empty, &MSExperiment::new()).unwrap();
    valid(file.is_valid(&empty).unwrap(), SchemaKind::MzData);
    // "test if filled file is valid": the source store discards what mzData
    // cannot hold, which is WriteOptions::source() here.
    let filled = directory.path().join("filled.mzData");
    let experiment = file.load(data("MzDataFile_1.mzData")).unwrap();
    mzdata::store_with_options(&filled, &experiment, &WriteOptions::source()).unwrap();
    valid(file.is_valid(&filled).unwrap(), SchemaKind::MzData);
}

// ===========================================================================
// MzXMLFile_test.cpp:586-600
// ===========================================================================

#[cfg(feature = "mzml")]
#[test]
fn mzxml_stored_experiment_is_valid() {
    use openms::format::mzxml::MzXMLFile;
    let file = MzXMLFile::new();
    // "Note: empty mzXML files are not valid, thus this test is omitted"
    // "test if full file is valid"
    let directory = temp();
    let path = directory.path().join("full.mzXML");
    let experiment = file.load(data("MzXMLFile_1.mzXML")).unwrap();
    file.store(&path, &experiment).unwrap();
    valid(file.is_valid(&path).unwrap(), SchemaKind::MzXML);
    // The fixture itself is mzXML 2.1, whose namespace the 3.1 schema does
    // not declare, so it is invalid against the schema the source registers.
    assert!(!file.is_valid(data("MzXMLFile_1.mzXML")).unwrap().is_valid());
}

// ===========================================================================
// MzIdentMLFile::isValid(filename, os, used_version), through FileInfo -v
// ===========================================================================

/// TOPP_FileInfo_14 and TOPP_FileInfo_15: "Validating mzid file against XML
/// schema version 1.1.0", then "Success" for the first and a pattern-facet
/// error at line 327 for the second, in the retained C++ outputs. The version
/// comes from `detectVersion`, as FileInfo's does.
#[cfg(feature = "idxml")]
#[test]
fn mzidentml_is_valid_detects_the_version_as_fileinfo_reports_it() {
    use openms::format::mzidentml;
    let r = mzidentml::is_valid(data("file_info/inputs/FileInfo_14_input.mzid")).unwrap();
    assert_eq!(r.schema.version(), Some("1.1.0"));
    valid(r, SchemaKind::MzIdentML1_1_0);
    let r = mzidentml::is_valid(data("file_info/inputs/FileInfo_15_input.mzid")).unwrap();
    assert_eq!(r.schema, SchemaKind::MzIdentML1_1_0);
    assert!(!r.is_valid());
    assert!(r.diagnostics.iter().any(|d| d.line == Some(327)), "{r:?}");
    // detectVersion reads the header within the validation limit.
    let mut o = openms::format::xml_schema::SchemaValidationOptions::default();
    o.limits.max_xml_bytes = 64;
    assert!(
        mzidentml::is_valid_with_options(data("file_info/inputs/FileInfo_14_input.mzid"), &o)
            .is_err()
    );
}

/// `PepXMLFile::isValid` validates against `pepXML_v114.xsd`
/// (`PepXMLFile.cpp:333`). The class test asserts no verdict for any pepXML
/// file, so only the wiring is checked: the entry point uses that schema and
/// answers with a report.
#[cfg(feature = "idxml")]
#[test]
fn pepxml_is_valid_uses_the_registered_schema() {
    use openms::format::pepxml;
    let r = pepxml::is_valid(data("PepXMLFile_test.pepxml")).unwrap();
    assert_eq!(r.schema, SchemaKind::PepXML);
    assert_eq!(r.schema.version(), Some("1.14"));
    // Another format's document has no global declaration in it.
    assert!(
        !pepxml::is_valid(data("transformation_xml_1.trafoXML"))
            .unwrap()
            .is_valid()
    );
}

/// The mzML path shares the engine: `ImzMLFile::isValid` and
/// `MzMLFile::isValid` keep their own tests; this only checks that the
/// re-exported report types are the shared ones.
#[cfg(feature = "mzml-schema")]
#[test]
fn mzml_reports_are_the_shared_report_type() {
    let r: SchemaValidationReport =
        openms::format::mzml::validate_schema(data("mzml_validator/MzMLFile_1.mzML")).unwrap();
    valid(r, SchemaKind::MzML);
}

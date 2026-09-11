#![cfg(feature = "mzml-validation")]
use openms::{
    data_structures::{
        CVMappingRule, CVMappingTerm, CVMappings, CombinationsLogic, RequirementLevel,
    },
    format::{
        controlled_vocabulary::ControlledVocabulary,
        cv_mapping::CVMappingFile,
        mzml,
        mzml_validator::{MzMLValidator, ParsedCVTerm, ValidationReport},
    },
};
use std::io::{Cursor, Write};

fn vocabulary() -> ControlledVocabulary {
    let mut cv = ControlledVocabulary::new();
    cv.load_obo_reader(
        "small",
        br#"[Term]
id: T:A
name: alpha
[Term]
id: T:B
name: beta
is_obsolete: true
[Term]
id: T:I
name: integer
xref: value-type:xsd\:integer "integer"
[Term]
id: T:V
name: value
relationship: has_units U:0
[Term]
id: U:0
name: unit
[Term]
id: U:1
name: child unit
is_a: U:0
[Term]
id: GO:1
name: gene
xref: value-type:xsd\:integer "integer"
[Term]
id: BTO:1
name: tissue
[Term]
id: MS:1000513
name: binary array
[Term]
id: MS:1000518
name: binary type
[Term]
id: A:0
name: array
is_a: MS:1000513
xref: binary-data-type:MS\:1000523 "float"
[Term]
id: MS:1000523
name: float
is_a: MS:1000518
[Term]
id: MS:1000519
name: integer type
is_a: MS:1000518
"#
        .as_slice(),
    )
    .unwrap();
    cv
}
fn mapping(paths: &[(&str, &[&str])], repeat: bool, level: RequirementLevel) -> CVMappings {
    let mut out = CVMappings::default();
    out.mapping_rules = paths
        .iter()
        .enumerate()
        .map(|(i, (path, ids))| CVMappingRule {
            identifier: format!("R{i}"),
            element_path: (*path).into(),
            requirement_level: level,
            combinations_logic: CombinationsLogic::Or,
            terms: ids
                .iter()
                .map(|id| CVMappingTerm {
                    accession: (*id).into(),
                    use_term: true,
                    is_repeatable: repeat,
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        })
        .collect();
    out
}
fn p(id: &str, name: &str) -> String {
    format!("<cvParam accession=\"{id}\" name=\"{name}\"/>")
}
fn validate(v: &MzMLValidator<'_>, xml: &str) -> ValidationReport {
    v.validate_reader(xml.as_bytes()).unwrap()
}

#[test]
fn unchanged_source_class_fixture_counts() {
    // Literal assertions in MzMLFile_test.cpp:1230–1251, not native writer output.
    let mapping = CVMappingFile::default()
        .read(include_bytes!("data/cv_mapping/ms-mapping.xml").as_slice())
        .unwrap();
    let cv = ControlledVocabulary::psi_ms().unwrap();
    let v = MzMLValidator::new(&mapping, cv);
    for bytes in [
        include_bytes!("data/mzml_validator/MzMLFile_1.mzML").as_slice(),
        include_bytes!("data/mzml_validator/MzMLFile_4_indexed.mzML").as_slice(),
    ] {
        let r = v.validate_reader(bytes).unwrap();
        assert_eq!(r, ValidationReport::default());
    }
    let r = v
        .validate_reader(include_bytes!("data/mzml_validator/MzMLFile_3_invalid.mzML").as_slice())
        .unwrap();
    assert!(!r.is_valid());
    assert_eq!((r.errors.len(), r.warnings.len()), (8, 1), "{r:#?}");
}
#[test]
fn units_default_and_all_configured_attribute_controls() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/param/@a", &["T:V", "T:I"])],
        true,
        RequirementLevel::Must,
    );
    let mut v = MzMLValidator::new(&map, &cv);
    assert!(v.options.check_units && v.options.check_term_value_types);
    v.options.tag = "param".into();
    v.options.accession_attribute = "a".into();
    v.options.name_attribute = "n".into();
    v.options.value_attribute = "v".into();
    v.options.unit_accession_attribute = "u".into();
    v.options.unit_name_attribute = "un".into();
    let xml = "<r><param a='T:V' n='value' u='U:1' un='child unit'/><param a='T:I' n='integer' v='42'/></r>";
    assert_eq!(validate(&v, xml), ValidationReport::default());
    assert_eq!(
        validate(&v, "<r><param a='T:V' n='value'/></r>").errors,
        ["CV term must have a unit: T:V - value"]
    );
    v.options.check_units = false;
    v.options.check_term_value_types = false;
    assert_eq!(
        validate(
            &v,
            "<r><param a='T:V' n='value'/><param a='T:I' n='integer' v='garbage'/></r>"
        ),
        ValidationReport::default()
    );
}
#[test]
fn groups_defer_validation_and_warn_at_definition_only() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:B"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let prefix = format!(
        "<r><referenceableParamGroup id='g'>{}{}</referenceableParamGroup>",
        p("T:B", "wrong"),
        p("T:missing", "missing")
    );
    let no_use = validate(&v, &format!("{prefix}</r>"));
    assert!(no_use.errors.is_empty());
    assert_eq!(
        no_use.warnings,
        [
            "Obsolete CV term: 'T:B - wrong' at element '/r/referenceableParamGroup'",
            "Unknown CV term: 'T:missing - missing' at element '/r/referenceableParamGroup'"
        ]
    );
    let r = validate(
        &v,
        &format!(
            "{prefix}<target><referenceableParamGroupRef ref='g'/><referenceableParamGroupRef ref='g'/></target></r>"
        ),
    );
    assert_eq!(r.warnings, no_use.warnings);
    assert_eq!(
        r.errors,
        [
            "Name of CV term not correct: 'T:B - wrong' should be 'beta'",
            "Name of CV term not correct: 'T:B - wrong' should be 'beta'"
        ]
    );
}
#[test]
fn duplicate_definitions_append_and_nonrepeatable_counts_include_every_use() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:A"])],
        false,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let defs = format!(
        "<referenceableParamGroup id='g'>{}</referenceableParamGroup>",
        p("T:A", "alpha")
    );
    let r = validate(
        &v,
        &format!("<r>{defs}{defs}<target><referenceableParamGroupRef ref='g'/></target></r>"),
    );
    assert!(r.warnings.is_empty());
    assert_eq!(
        r.errors,
        ["Violated mapping rule 'R0' number of term repeats at element '/r/target'"]
    );
}
#[test]
fn forward_and_missing_references_are_empty_and_no_state_survives_calls() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let def = format!(
        "<referenceableParamGroup id='g'>{}</referenceableParamGroup>",
        p("T:A", "alpha")
    );
    let missing = "<r><target><referenceableParamGroupRef ref='g'/></target></r>";
    let expected = validate(&v, missing);
    assert_eq!(
        expected.errors,
        ["Violated mapping rule 'R0' at element '/r/target', at least one term must be present!"]
    );
    assert_eq!(
        validate(
            &v,
            &format!("<r>{def}<target><referenceableParamGroupRef ref='g'/></target></r>")
        ),
        ValidationReport::default()
    );
    assert_eq!(validate(&v, missing), expected);
    assert_eq!(
        validate(
            &v,
            &format!("<r><target><referenceableParamGroupRef ref='g'/></target>{def}</r>")
        ),
        expected
    );
    assert!(
        v.validate_reader(
            format!("<r>{def}<target><referenceableParamGroupRef ref='g'/></wrong>").as_bytes()
        )
        .is_err()
    );
    assert_eq!(validate(&v, missing), expected);
}
#[test]
fn nested_group_current_id_is_not_restored_and_user_params_do_not_enter_groups() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:A"])],
        false,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let xml = format!(
        "<r><referenceableParamGroup id='outer'><referenceableParamGroup id='inner'>{}</referenceableParamGroup><userParam/><cvParam accession='T:A' name='alpha'/></referenceableParamGroup><target><referenceableParamGroupRef ref='inner'/></target></r>",
        p("T:A", "alpha")
    );
    assert_eq!(
        validate(&v, &xml).errors,
        ["Violated mapping rule 'R0' number of term repeats at element '/r/target'"]
    );
    assert!(
        v.validate_reader(b"<r><referenceableParamGroup/></r>".as_slice())
            .is_err()
    );
    assert!(
        v.validate_reader(b"<r><referenceableParamGroupRef/></r>".as_slice())
            .is_err()
    );
}
#[test]
fn indexed_root_paths_and_literal_nested_names() {
    let cv = vocabulary();
    let map = mapping(
        &[("/mzML/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    assert_eq!(
        validate(
            &v,
            &format!(
                "<indexedmzML><mzML>{}</mzML></indexedmzML>",
                p("T:A", "alpha")
            )
        ),
        ValidationReport::default()
    );
    // Only the leading indexed wrapper is omitted, never a same-named child.
    let map = mapping(
        &[("/indexedmzML/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    assert_eq!(
        validate(
            &v,
            &format!(
                "<indexedmzML><indexedmzML>{}</indexedmzML></indexedmzML>",
                p("T:A", "alpha")
            )
        ),
        ValidationReport::default()
    );
    let map = mapping(
        &[("//cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::May,
    );
    let v = MzMLValidator::new(&map, &cv);
    assert_eq!(
        validate(&v, &p("T:A", "alpha")),
        ValidationReport::default()
    );
}
#[test]
fn go_bto_skip_after_unknown_obsolete_encounter_and_do_not_fulfill_rules() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/cvParam/@accession", &["GO:1", "BTO:1"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let r = validate(
        &v,
        "<r><cvParam accession='GO:1' name='wrong' value='not integer'/><cvParam accession='BTO:1' name='wrong' unitAccession='absent'/><cvParam accession='GO:missing' name='missing'/></r>",
    );
    assert_eq!(
        r.errors,
        ["Violated mapping rule 'R0' at element '/r', at least one term must be present!"]
    );
    assert_eq!(
        r.warnings,
        ["Unknown CV term: 'GO:missing - missing' at element '/r'"]
    );
}
#[test]
fn binary_pair_order_group_expansion_repeated_errors_and_reset() {
    let cv = vocabulary();
    let map = mapping(
        &[(
            "/r/binaryDataArray/cvParam/@accession",
            &["A:0", "MS:1000523", "MS:1000519", "T:A"],
        )],
        true,
        RequirementLevel::May,
    );
    let v = MzMLValidator::new(&map, &cv);
    let role = p("A:0", "array");
    let good = p("MS:1000523", "float");
    let bad = p("MS:1000519", "integer type");
    let extra = p("T:A", "alpha");
    for pair in [format!("{role}{good}"), format!("{good}{role}")] {
        assert_eq!(
            validate(
                &v,
                &format!("<r><binaryDataArray>{pair}</binaryDataArray></r>")
            ),
            ValidationReport::default()
        );
    }
    let r = validate(
        &v,
        &format!(
            "<r><referenceableParamGroup id='g'>{bad}{role}</referenceableParamGroup><binaryDataArray><referenceableParamGroupRef ref='g'/>{extra}</binaryDataArray><binaryDataArray>{good}</binaryDataArray></r>"
        ),
    );
    let expected = "Binary data array of type 'A:0 ! array' cannot have the value type 'MS:1000519 ! integer type'.";
    assert_eq!(r.errors, [expected, expected]);
    assert!(r.warnings.is_empty());
    let map = mapping(
        &[(
            "/r/binaryDataArray/param/@accession",
            &["A:0", "MS:1000519"],
        )],
        true,
        RequirementLevel::May,
    );
    let mut v = MzMLValidator::new(&map, &cv);
    v.options.tag = "param".into();
    // Source hard-codes cvParam in this one suffix even with a configured tag.
    assert_eq!(
        validate(
            &v,
            &format!("<r><binaryDataArray>{bad}{role}</binaryDataArray></r>")
                .replace("cvParam", "param")
        ),
        ValidationReport::default()
    );
}
#[test]
fn locate_selection_keeps_stable_missing_path_error_and_ignores_mzml_exclusions() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/cvParam/@accession", &["GO:1"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    let term = ParsedCVTerm {
        accession: "GO:1".into(),
        ..Default::default()
    };
    assert!(v.locate_term("/r/cvParam/@accession", &term).unwrap());
    assert!(v.locate_term("/absent/cvParam/@accession", &term).is_err());
    validate(&v, "<absent/>");
    assert!(v.locate_term("/absent/cvParam/@accession", &term).is_err());
}
#[test]
fn shared_limits_count_referenced_uses_and_late_diagnostics() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::Must,
    );
    let mut v = MzMLValidator::new(&map, &cv);
    let xml = format!(
        "<r><referenceableParamGroup id='g'>{}</referenceableParamGroup><target><referenceableParamGroupRef ref='g'/><referenceableParamGroupRef ref='g'/></target></r>",
        p("T:A", "alpha")
    );
    v.options.limits.max_terms = 3;
    assert!(validate(&v, &xml).is_valid());
    v.options.limits.max_terms = 2;
    assert!(v.validate_reader(xml.as_bytes()).is_err());
    v.options.limits.max_terms = 100;
    v.options.limits.max_diagnostics = 1;
    let xml = format!(
        "<r><target>{}{}</target></r>",
        p("T:A", "wrong"),
        p("T:A", "wrong")
    );
    assert!(v.validate_reader(xml.as_bytes()).is_err());
    v.options.limits.max_diagnostics = 2;
    assert_eq!(validate(&v, &xml).errors.len(), 2);
}
#[test]
fn lexical_errors_and_compressed_paths_use_the_shared_parser() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::May,
    );
    let v = MzMLValidator::new(&map, &cv);
    for xml in [
        "",
        "<",
        "<r><cvParam accession='T:A'name='alpha'/></r>",
        "<r>]]></r>",
        "<!DOCTYPE r SYSTEM 'https://example.invalid/a'><r/>",
        "<?xml version='1.0' standalone='yes' encoding='UTF-8'?><r/>",
        "<r/><r/>",
    ] {
        assert!(v.validate_reader(xml.as_bytes()).is_err(), "{xml}");
    }
    let raw = b"<r><cvParam accession='T:A' name='alpha'/></r>";
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    gzip.write_all(raw).unwrap();
    let file = openms::system::file::TempFile::new().unwrap();
    std::fs::write(file.path(), gzip.finish().unwrap()).unwrap();
    assert_eq!(
        v.validate(file.path()).unwrap(),
        ValidationReport::default()
    );
    assert_eq!(
        v.validate_reader(Cursor::new(raw)).unwrap(),
        ValidationReport::default()
    );
}
#[test]
fn default_mapping_convenience_and_native_writer_are_distinct_from_source_goldens() {
    let file = openms::system::file::TempFile::new().unwrap();
    std::fs::write(
        file.path(),
        include_bytes!("data/mzml_validator/MzMLFile_1.mzML"),
    )
    .unwrap();
    assert_eq!(
        mzml::validate_semantics(file.path()).unwrap(),
        ValidationReport::default()
    );
    let options = openms::format::mzml_validator::ValidationOptions {
        limits: openms::format::mzml_validator::ValidationLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(mzml::validate_semantics_with_options(file.path(), &options).is_err());
    assert_eq!(options.limits.max_input_bytes, 1);
    // Native writer regression is not assigned the source writer's warning count.
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &openms::MSExperiment::default()).unwrap();
    let mapping = CVMappingFile::default()
        .read(include_bytes!("data/cv_mapping/ms-mapping.xml").as_slice())
        .unwrap();
    let r = MzMLValidator::new(&mapping, ControlledVocabulary::psi_ms().unwrap())
        .validate_reader(bytes.as_slice())
        .unwrap();
    assert!(r.is_valid(), "{r:#?}");
}

#[test]
fn fixed_special_tags_win_over_configured_cv_tag_and_qnames_are_literal() {
    let cv = vocabulary();
    let empty = CVMappings::default();
    let mut v = MzMLValidator::new(&empty, &cv);
    for (tag, xml) in [
        (
            "referenceableParamGroup",
            "<r><referenceableParamGroup id='g'/></r>",
        ),
        (
            "referenceableParamGroupRef",
            "<r><referenceableParamGroupRef ref='g'/></r>",
        ),
        ("binaryDataArray", "<r><binaryDataArray/></r>"),
    ] {
        v.options.tag = tag.into();
        assert_eq!(validate(&v, xml), ValidationReport::default());
    }
    let map = mapping(
        &[("/p:indexedmzML/mzML/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::Must,
    );
    let v = MzMLValidator::new(&map, &cv);
    assert_eq!(
        validate(
            &v,
            "<p:indexedmzML xmlns:p='urn:test'><mzML><cvParam accession='T:A' name='alpha'/></mzML></p:indexedmzML>"
        ),
        ValidationReport::default()
    );
}

#[test]
fn repeated_borrowed_group_terms_exhaust_shared_work_after_a_successful_use() {
    let cv = vocabulary();
    let map = mapping(
        &[("/r/target/cvParam/@accession", &["T:A"])],
        true,
        RequirementLevel::May,
    );
    let mut v = MzMLValidator::new(&map, &cv);
    v.options.check_term_value_types = false;
    let definition = format!(
        "<referenceableParamGroup id='g'><cvParam accession='T:A' name='alpha' value='{}'/></referenceableParamGroup>",
        "x".repeat(8192)
    );
    let once = format!("<r>{definition}<target><referenceableParamGroupRef ref='g'/></target></r>");
    let many = format!(
        "<r>{definition}<target>{}</target></r>",
        "<referenceableParamGroupRef ref='g'/>".repeat(100)
    );
    // The large value is stored once; every borrowed application must still
    // charge its inspection rather than getting a fresh per-reference budget.
    let (mut low, mut high) = (0, 1_000_000);
    while low < high {
        let middle = (low + high) / 2;
        v.options.limits.max_work = middle;
        if v.validate_reader(once.as_bytes()).is_ok() {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    assert!(low < 1_000_000);
    v.options.limits.max_work = low;
    assert!(validate(&v, &once).is_valid());
    assert!(v.validate_reader(many.as_bytes()).is_err());
    v.options.limits.max_work = 2_000_000;
    assert!(validate(&v, &many).is_valid());
}

#[test]
fn original_source_resources_retain_exact_bytes_on_every_platform() {
    for (bytes, length, hash) in [
        (
            include_bytes!("data/mzml_validator/MzMLFile_1.mzML").as_slice(),
            37187,
            0x865ec686f46133f3_u64,
        ),
        (
            include_bytes!("data/mzml_validator/MzMLFile_3_invalid.mzML").as_slice(),
            22820,
            0xd236d1a690271f3f_u64,
        ),
        (
            include_bytes!("data/mzml_validator/MzMLFile_4_indexed.mzML").as_slice(),
            27526,
            0x4a68d9fbeb7f3888_u64,
        ),
        (
            include_bytes!("../resources/cv/ms-mapping.xml").as_slice(),
            21470,
            0x555ce91d7382c91b_u64,
        ),
    ] {
        assert_eq!(bytes.len(), length);
        let actual = bytes.iter().fold(14695981039346656037_u64, |h, b| {
            (h ^ u64::from(*b)).wrapping_mul(1099511628211)
        });
        assert_eq!(actual, hash);
    }
}

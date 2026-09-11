// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::data_structures::*;

#[test]
fn complete_source_defaults_and_literal_fields() {
    assert_eq!(
        CVReference::default(),
        CVReference {
            name: String::new(),
            identifier: String::new()
        }
    );
    let t = CVMappingTerm::default();
    assert!(!t.use_term_name && !t.use_term && !t.is_repeatable && !t.allow_children);
    assert!(t.accession.is_empty() && t.term_name.is_empty() && t.cv_identifier_ref.is_empty());
    let r = CVMappingRule::default();
    assert_eq!(r.requirement_level, RequirementLevel::Must);
    assert_eq!(r.combinations_logic, CombinationsLogic::Or);
    assert_eq!(
        [
            RequirementLevel::Must as u8,
            RequirementLevel::Should as u8,
            RequirementLevel::May as u8
        ],
        [0, 1, 2]
    );
    assert_eq!(
        [
            CombinationsLogic::Or as u8,
            CombinationsLogic::And as u8,
            CombinationsLogic::Xor as u8
        ],
        [0, 1, 2]
    );
    assert!(
        r.identifier.is_empty()
            && r.element_path.is_empty()
            && r.scope_path.is_empty()
            && r.terms.is_empty()
    );
    let original = CVMappingRule {
        identifier: "my_test_identifier".into(),
        element_path: "my_test_elementpath".into(),
        scope_path: "my_test_scopepath".into(),
        requirement_level: RequirementLevel::Should,
        combinations_logic: CombinationsLogic::Xor,
        terms: vec![CVMappingTerm {
            accession: "BLA:1".into(),
            ..t
        }],
    };
    let mut copy = original.clone();
    assert_eq!(copy, original);
    copy.terms[0].term_name = "my_test_termname".into();
    assert_ne!(copy, original);
    assert!(original.terms[0].term_name.is_empty());
    copy.terms = vec![CVMappingTerm {
        accession: "BLA:2".into(),
        ..Default::default()
    }];
    copy.terms.push(CVMappingTerm {
        accession: "BLA:3".into(),
        ..Default::default()
    });
    assert_eq!(copy.terms.len(), 2);
}

#[test]
fn equality_observes_every_leaf_and_rule_field() {
    let term = CVMappingTerm::default();
    let changes: [fn(&mut CVMappingTerm); 7] = [
        |v| v.accession = "a".into(),
        |v| v.use_term_name = true,
        |v| v.use_term = true,
        |v| v.term_name = "n".into(),
        |v| v.is_repeatable = true,
        |v| v.allow_children = true,
        |v| v.cv_identifier_ref = "c".into(),
    ];
    for change in changes {
        let mut other = term.clone();
        change(&mut other);
        assert_ne!(term, other);
        assert_eq!(other, other.clone());
    }
    let rule = CVMappingRule::default();
    let changes: [fn(&mut CVMappingRule); 6] = [
        |v| v.identifier = "a".into(),
        |v| v.element_path = "p".into(),
        |v| v.requirement_level = RequirementLevel::May,
        |v| v.scope_path = "s".into(),
        |v| v.combinations_logic = CombinationsLogic::And,
        |v| v.terms.push(CVMappingTerm::default()),
    ];
    for change in changes {
        let mut other = rule.clone();
        change(&mut other);
        assert_ne!(rule, other);
    }
    let a = CVReference {
        name: "my_test_name".into(),
        ..Default::default()
    };
    assert_ne!(a, CVReference::default());
    let mut b = a.clone();
    b.identifier = "my_test_identifier".into();
    assert_ne!(a, b);
}

#[test]
fn ordered_bulk_append_duplicate_add_and_explicit_replacement() {
    let a = CVReference {
        name: "first".into(),
        identifier: "Ref1".into(),
    };
    let b = CVReference {
        name: "second".into(),
        identifier: "Ref2".into(),
    };
    let replacement = CVReference {
        name: "last".into(),
        identifier: "Ref1".into(),
    };
    let mut m = CVMappings::default();
    m.set_cv_references(vec![a.clone(), b.clone()]);
    assert!(!m.add_cv_reference(replacement.clone()));
    m.set_cv_references(vec![replacement.clone()]);
    m.set_cv_references(vec![]);
    assert_eq!(m.cv_references(), [a, b, replacement]);
    assert!(m.has_cv_reference("Ref1") && m.has_cv_reference("Ref2"));
    assert!(!m.has_cv_reference("Ref3"));
    // Owned input prevents source CPP-034 iterator invalidation on self alias.
    let own = m.cv_references().to_vec();
    m.set_cv_references(own);
    assert_eq!(m.cv_references().len(), 6);
    let mut copy = m.clone();
    copy.mapping_rules.push(CVMappingRule::default());
    assert_ne!(copy, m);
    copy.mapping_rules.clear();
    assert_eq!(copy, m);
    m.replace_cv_references(vec![]);
    assert!(m.cv_references().is_empty() && !m.has_cv_reference("Ref1"));
    assert_eq!(copy.cv_references().len(), 6);
}

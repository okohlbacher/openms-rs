// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Independent pinned-source projections and synthetic ontology branch oracles.
//! No test reads a C++ checkout. See data/registry_provenance.json.

use openms::chemistry::{
    AASequence, CrossLinksDB, EmpiricalFormula, ModificationsDB, ModifiedPeptideGenerator,
    OboReadOptions, TermSpecificity,
};
use std::collections::BTreeSet;

const XLMOD: &[u8] = include_bytes!("../resources/modifications/XLMOD.obo");
const BASE: &str = include_str!("data/registry_synthetic_base.tsv");
const SYNTHETIC: &str = include_str!("data/registry_synthetic_psi.obo");
const BOUNDARIES: &str = include_str!("data/registry_synthetic_boundaries.obo");

fn full_ids(db: &ModificationsDB) -> Vec<&str> {
    db.entries().iter().map(|m| m.full_id()).collect()
}
fn formula(text: &str) -> EmpiricalFormula {
    EmpiricalFormula::parse(text).unwrap()
}
fn assert_projection(db: &ModificationsDB, mode: &str) {
    let rows: Vec<_> = include_str!("data/registry_xlmod_reference.tsv")
        .lines()
        .filter(|line| line.starts_with(mode))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!(db.len(), rows.len(), "{mode}");
    for (modification, row) in db.entries().iter().zip(rows) {
        assert_eq!(modification.obo_accession(), Some(row[2]));
        assert_eq!(modification.record_id(), None);
        assert_eq!(modification.name(), row[3]);
        assert_eq!(modification.full_name(), row[3]);
        assert_eq!(modification.full_id(), row[4]);
        let origin = if row[5] == "X" && row[6] != "Anywhere" {
            None // Native terminal wildcard, equivalent to source origin X.
        } else {
            row[5].chars().next()
        };
        assert_eq!(modification.origin(), origin, "{}", row[4]);
        assert_eq!(modification.term_specificity().name(), row[6]);
        assert_eq!(
            modification.diff_mono_mass().to_bits(),
            u64::from_str_radix(row[8], 16).unwrap(),
            "{}",
            row[4]
        );
        assert!(modification.diff_formula().is_empty());
        assert!(modification.absolute_formula().is_none());
        assert_eq!(modification.mono_mass(), 0.);
        let synonyms: BTreeSet<String> = row[9]
            .split(';')
            .filter(|s| !s.is_empty())
            .map(str::to_owned)
            .collect();
        assert_eq!(modification.synonyms(), &synonyms, "{}", row[4]);
    }
}

#[test]
fn complete_pinned_xlmod_projection_in_source_accession_and_site_order() {
    let mono = ModificationsDB::from_obo(XLMOD, &OboReadOptions::default()).unwrap();
    assert_projection(&mono, "monolink");
    // CrossLinksDB forces the cross-link projection even with default options.
    let cross = CrossLinksDB::from_obo(XLMOD, &OboReadOptions::default()).unwrap();
    assert_projection(cross.database(), "crosslink");
    assert_projection(CrossLinksDB::global().database(), "crosslink");
    let mut sorted = full_ids(cross.database());
    sorted.sort();
    assert_eq!(cross.all_search_modifications(), sorted);
    let global = ModificationsDB::global();
    assert_eq!(global.len(), 3035 + 92);
    for expected in mono.entries() {
        assert_eq!(
            global
                .get_modification(expected.full_id(), None, None)
                .unwrap(),
            expected.as_ref()
        );
    }
}

#[test]
fn pinned_crosslink_class_test_mass_and_specificity_goldens() {
    let db = CrossLinksDB::global().database();
    for term in [TermSpecificity::Anywhere, TermSpecificity::NTerm] {
        assert_eq!(db.find("DSS", Some('K'), Some(term)).len(), 1);
    }
    let mut matches: Vec<_> = db
        .search_by_mass(138.068_079_61, 0.00001, Some('K'), None)
        .unwrap()
        .into_iter()
        .map(|m| m.full_id())
        .collect();
    matches.sort();
    assert_eq!(
        matches,
        ["BS3 (K)", "BS3 (N-term)", "DSS (K)", "DSS (N-term)"]
    );
    assert_eq!(
        db.find("EDC", None, None)
            .iter()
            .map(|m| m.full_id())
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "EDC (D)",
            "EDC (E)",
            "EDC (K)",
            "EDC (S)",
            "EDC (T)",
            "EDC (Y)",
            "EDC (N-term)",
            "EDC (C-term)",
        ])
    );
}

#[test]
fn psi_aliases_bind_all_existing_unimod_specificities_and_drop_missing_targets() {
    let mut db = ModificationsDB::from_tsv(BASE).unwrap();
    let report = db
        .extend_obo(SYNTHETIC.as_bytes(), &OboReadOptions::default())
        .unwrap();
    assert_eq!(report.records_added, 5);
    assert_eq!(report.unresolved_aliases, 1);
    let aliases = db.find("MOD:9000001", None, None);
    assert_eq!(aliases.len(), 2);
    assert_eq!(
        aliases.iter().map(|m| m.full_id()).collect::<Vec<_>>(),
        ["ReferenceAnchor (M)", "ReferenceAnchor (K)"]
    );
    // The alias's declared N-terminal M specificity is intentionally ignored.
    assert!(
        aliases
            .iter()
            .all(|m| m.term_specificity() == TermSpecificity::Anywhere)
    );
    assert!(db.find("Reference resolved alias", None, None).is_empty());
    assert!(db.find("Reference alias label", None, None).is_empty());
    assert!(db.find("MOD:9000002", None, None).is_empty());
    let without_base =
        ModificationsDB::from_obo(SYNTHETIC.as_bytes(), &OboReadOptions::default()).unwrap();
    assert_eq!(without_base.len(), 5);
    assert!(without_base.find("MOD:9000001", None, None).is_empty());
    let no_change = db.get_modification("MOD:9000010", None, None).unwrap();
    assert_eq!(no_change.name(), "MOD:9000010");
    assert_eq!(
        no_change.full_id(),
        "Reference unchanged absolute record (M)"
    );
    assert_eq!(db.find("Reference no change", None, None), vec![no_change]);
    assert_eq!(db.find("Reference final synonym", None, None).len(), 1);
}

#[test]
fn source_origin_filter_terminal_expansion_and_reaction_site_complement() {
    let mono =
        ModificationsDB::from_obo(BOUNDARIES.as_bytes(), &OboReadOptions::default()).unwrap();
    assert_eq!(
        full_ids(&mono),
        [
            "Reference union of sites (K)",
            "Reference union of sites (M)",
            "Reference union of sites (Protein N-term)",
            "Reference union of sites (Protein C-term)",
            "Reference mono only (K)",
        ]
    );
    let cross = CrossLinksDB::from_obo(BOUNDARIES.as_bytes(), &OboReadOptions::default()).unwrap();
    assert_eq!(
        full_ids(cross.database()),
        [
            "Reference union of sites (K)",
            "Reference union of sites (M)",
            "Reference union of sites (N-term)",
            "Reference union of sites (C-term)",
            "Reference cross only (K)",
            "Reference cross only (M)",
        ]
    );
    // A count other than 1/2 is not silently reduced to either chemistry class.
    assert_eq!(mono.find("Reference site alias", None, None).len(), 4);
    assert_eq!(
        cross
            .database()
            .find("Reference site alias", None, None)
            .len(),
        4
    );
}

#[test]
fn eof_flush_and_unknown_stanza_isolation() {
    let no_newline = SYNTHETIC.trim_end();
    let db = ModificationsDB::from_obo(no_newline.as_bytes(), &OboReadOptions::default()).unwrap();
    assert_eq!(
        db.get_modification("MOD:9000014", None, None)
            .unwrap()
            .diff_mono_mass(),
        4.25
    );
    // Native correction: a non-Term stanza cannot overwrite a preceding term.
    let text = "[Term]\nid: MOD:9000090\nname: Before typedef\nproperty_value: Origin: \"K\" xsd:string\n[Typedef]\nid: ignored_relation\nname: Not a modification\nproperty_value: Origin: \"M\" xsd:string\n[Term]\nid: MOD:9000091\nname: After typedef\nproperty_value: Origin: \"M\" xsd:string";
    let db = ModificationsDB::from_obo(text.as_bytes(), &OboReadOptions::default()).unwrap();
    assert_eq!(full_ids(&db), ["Before typedef (K)", "After typedef (M)"]);
    assert!(db.find("ignored_relation", None, None).is_empty());
}

#[test]
fn malformed_records_and_all_parser_budgets_leave_existing_registry_unchanged() {
    let valid =
        "[Term]\nid: MOD:9000080\nname: Valid prefix\nproperty_value: Origin: \"M\" xsd:string\n";
    for invalid in [
        "property_value: DiffMono: \"NaN\" xsd:float",
        "property_value: MassMono: \"inf\" xsd:float",
        "property_value: DiffAvg: \"invalid\" xsd:float",
        "property_value: DiffFormula: \"Qq\" xsd:string",
        "property_value: TermSpec: \"not-a-terminus\" xsd:string",
        "property_value: DiffMono: missing-quotes xsd:float",
        "synonym: missing-quotes EXACT []",
    ] {
        let mut db = ModificationsDB::from_tsv(BASE).unwrap();
        let before = db.entries().to_vec();
        let text = format!(
            "{valid}[Term]\nid: MOD:9000081\nname: Malformed suffix\nproperty_value: Origin: \"M\" xsd:string\n{invalid}\n"
        );
        assert!(
            db.extend_obo(text.as_bytes(), &OboReadOptions::default())
                .is_err(),
            "{invalid}"
        );
        assert_eq!(db.entries(), before);
        assert!(db.find("MOD:9000080", None, None).is_empty());
        assert_eq!(db.find("UniMod:990001", None, None).len(), 2);
    }
    let limits = [
        OboReadOptions {
            max_input_bytes: SYNTHETIC.len() - 1,
            ..Default::default()
        },
        OboReadOptions {
            max_line_bytes: 10,
            ..Default::default()
        },
        OboReadOptions {
            max_terms: 1,
            ..Default::default()
        },
        OboReadOptions {
            max_records: 1,
            ..Default::default()
        },
        OboReadOptions {
            max_aliases: 0,
            ..Default::default()
        },
        OboReadOptions {
            max_registry_bytes: 1,
            ..Default::default()
        },
    ];
    for options in limits {
        let mut db = ModificationsDB::from_tsv(BASE).unwrap();
        let before = db.entries().to_vec();
        assert!(
            db.extend_obo(SYNTHETIC.as_bytes(), &options).is_err(),
            "{options:?}"
        );
        assert_eq!(db.entries(), before);
        assert!(db.find("MOD:9000001", None, None).is_empty());
    }
}

#[test]
fn repeated_alias_target_expansion_is_bounded_even_when_results_deduplicate() {
    // Ten stored specificity records have forty existing lookup associations.
    // One alias adds ten associations; subsequent occurrences add none, but
    // still visit all ten targets and must consume the expansion-work budget.
    let base = BASE.repeat(5);
    let alias = "[Term]\nid: MOD:9000070\nname: Repeated alias\ndef: \"Synthetic.\" [UniMod:990001]\nproperty_value: Origin: \"M\" xsd:string\n";
    let options = OboReadOptions {
        max_aliases: 50,
        ..Default::default()
    };
    let mut at_limit = ModificationsDB::from_tsv(&base).unwrap();
    let report = at_limit
        .extend_obo(alias.repeat(5).as_bytes(), &options)
        .unwrap();
    assert_eq!(report.records_added, 0);
    assert_eq!(report.aliases_added, 10);
    assert_eq!(at_limit.find("MOD:9000070", None, None).len(), 10);

    let mut beyond = ModificationsDB::from_tsv(&base).unwrap();
    let before = beyond.entries().to_vec();
    assert!(
        beyond
            .extend_obo(alias.repeat(6).as_bytes(), &options)
            .is_err()
    );
    assert_eq!(beyond.entries(), before);
    assert!(beyond.find("MOD:9000070", None, None).is_empty());
}

#[test]
fn absolute_formula_no_change_delta_precedence_and_retained_registry_ownership() {
    let plain = AASequence::parse("M").unwrap();
    let db = ModificationsDB::from_obo(SYNTHETIC.as_bytes(), &OboReadOptions::default()).unwrap();
    let unchanged = AASequence::parse_with_registry("M(MOD:9000010)", &db).unwrap();
    assert_eq!(unchanged.formula().unwrap(), plain.formula().unwrap());
    assert_eq!(
        unchanged.mono_mass().unwrap().to_bits(),
        plain.mono_mass().unwrap().to_bits()
    );
    let replaced = AASequence::parse_with_registry("M(MOD:9000011)", &db).unwrap();
    assert_eq!(replaced.formula().unwrap(), formula("C4H9NO3"));
    assert_eq!(
        replaced.mono_mass().unwrap(),
        formula("C4H9NO3").mono_mass()
    );
    let delta = AASequence::parse_with_registry("M(MOD:9000012)", &db).unwrap();
    assert_eq!(delta.formula().unwrap(), formula("C5H11NO3S"));
    let mut terminal = plain.clone();
    terminal
        .set_n_terminal_modification_with_registry("MOD:9000013", &db)
        .unwrap();
    assert_eq!(
        terminal.mono_mass().unwrap(),
        plain.mono_mass().unwrap() + 12.5
    );
    assert!(terminal.formula().is_err());
    let handles =
        ModifiedPeptideGenerator::get_modifications_with_registry(&["MOD:9000011"], &db).unwrap();
    drop(db);
    let mut generated = plain;
    ModifiedPeptideGenerator::default()
        .apply_fixed_modifications(&handles, &mut generated)
        .unwrap();
    assert_eq!(generated, replaced);
    assert_eq!(
        unchanged
            .residue_modification(0)
            .unwrap()
            .unwrap()
            .full_id(),
        "Reference unchanged absolute record (M)"
    );
}

#[test]
fn empty_annotated_unimod_export_and_unrepresentable_terminal_mass() {
    let acetyl = ModifiedPeptideGenerator::get_modifications(&["Acetyl (N-term)"]).unwrap();
    let mut empty = AASequence::default();
    ModifiedPeptideGenerator::default()
        .apply_fixed_modifications(&acetyl, &mut empty)
        .unwrap();
    assert!(empty.n_terminal_modification().is_some());
    // Literal source toUniModString returns immediately for an empty sequence.
    assert_eq!(empty.to_unimod_string().unwrap(), "");

    let text = "[Term]\nid: MOD:9000099\nname: Reference negative terminal\nproperty_value: Origin: \"X\" xsd:string\nproperty_value: TermSpec: \"N-term\" xsd:string\nproperty_value: DiffMono: \"-2\" xsd:float";
    let db = ModificationsDB::from_obo(text.as_bytes(), &OboReadOptions::default()).unwrap();
    let sequence = AASequence::parse_with_registry(".(MOD:9000099)M", &db).unwrap();
    assert!(sequence.mono_mass().unwrap() > 0.);
    // Native bracket signs mean delta. A negative absolute terminal component
    // cannot be emitted faithfully by this representation and is checked.
    assert!(sequence.to_unimod_string().is_err());
}

#[test]
fn literal_empty_absolute_formula_is_absent_but_explicit_zero_composition_is_retained() {
    let prefix = "[Term]\nid: MOD:9000098\nname: Reference empty formula\nproperty_value: Origin: \"M\" xsd:string\nproperty_value: DiffMono: \"12.5\" xsd:float\n";
    let empty = format!("{prefix}property_value: Formula: \"\" xsd:string");
    let db = ModificationsDB::from_obo(empty.as_bytes(), &OboReadOptions::default()).unwrap();
    assert!(db.entries()[0].absolute_formula().is_none());
    let peptide = AASequence::parse_with_registry("M(MOD:9000098)", &db).unwrap();
    let water = formula("H2O").mono_mass();
    let expected_internal = (formula("C5H11NO2S").mono_mass() + 12.5) - water;
    assert_eq!(
        peptide.mono_mass().unwrap().to_bits(),
        (water + expected_internal).to_bits()
    );
    assert!(peptide.formula().is_err());
    // Source retains a nonempty absolute formula string even when its parsed
    // atom counts are zero. This must not silently become the absent case.
    for (literal, charge) in [("H0", 0), ("+", 1)] {
        let text = format!("{prefix}property_value: Formula: \"{literal}\" xsd:string");
        let db = ModificationsDB::from_obo(text.as_bytes(), &OboReadOptions::default()).unwrap();
        let absolute = db.entries()[0].absolute_formula().unwrap();
        assert!(absolute.is_empty());
        assert_eq!(absolute.charge(), charge);
        assert!(AASequence::parse_with_registry("M(MOD:9000098)", &db).is_err());
    }
}

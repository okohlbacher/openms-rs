// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "unimod-xml")]

//! `FORMAT/UnimodXMLFile.h`, its handler and `CHEMISTRY/UnimodXMLDataProvider.h`.
//!
//! Tier 1: `tests/data/unimod_provider_release.tsv` is the output of
//! `UnimodXMLDataProvider::loadModifications` in the executed Release build
//! (`../oracle/A5a-unimod`, ibminode06, Linux x86_64 Release
//! bc9cc12/c19e494/174b576) on the installed `unimod.xml` and
//! `custom_mods.xml`, which are byte-identical to
//! `resources/modifications/` (git blobs `27c7c09e` and `b5dac074` at core
//! `bc9cc12`). Every record is compared field by field. Tier 4: the records
//! agree with the compiled table `tools/generate_modifications.py` derives
//! from the same files, except where the source's classification rule is
//! documented to differ.

use openms::chemistry::{ModificationsDB, ResidueModification, TermSpecificity};
use openms::format::unimod_xml::{UnimodXmlDataProvider, UnimodXmlFile, records_to_modifications};

const ORACLE: &str = include_str!("data/unimod_provider_release.tsv");

fn load(file: &str) -> Vec<ResidueModification> {
    let path = format!(
        "{}/resources/modifications/{file}",
        env!("CARGO_MANIFEST_DIR")
    );
    UnimodXmlDataProvider::new(path)
        .load_modifications()
        .unwrap()
}

fn losses(modification: &ResidueModification) -> String {
    modification
        .neutral_losses()
        .iter()
        .map(|loss| {
            format!(
                "{}@{}@{}",
                loss.formula().to_source_string(),
                loss.mono_mass(),
                loss.average_mass()
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn oracle_losses(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    text.split(';')
        .map(|loss| {
            let parts: Vec<&str> = loss.split('@').collect();
            let mono: f64 = parts[1].parse().unwrap();
            let average: f64 = parts[2].parse().unwrap();
            format!("{}@{mono}@{average}", parts[0])
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn provider_matches_the_release_build_record_for_record() {
    for file in ["unimod.xml", "custom_mods.xml"] {
        let rows: Vec<Vec<&str>> = ORACLE
            .lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|fields| fields[0] == file)
            .collect();
        let count: usize = ORACLE
            .lines()
            .find_map(|line| line.strip_prefix(&format!("count\t{file}\t")))
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(rows.len(), count, "{file}");
        let loaded = load(file);
        assert_eq!(loaded.len(), count, "{file}");
        for (row, modification) in rows.iter().zip(&loaded) {
            let context = format!("{file}: {}", row[2]);
            assert_eq!(modification.name(), row[1], "{context}");
            assert_eq!(modification.full_id(), row[2], "{context}");
            assert_eq!(
                modification.record_id().unwrap().to_string(),
                row[3],
                "{context}"
            );
            assert_eq!(
                modification.origin().unwrap_or('X').to_string(),
                row[4],
                "{context}"
            );
            // getTermSpecificityName() spells ANYWHERE "none".
            let term = match modification.term_specificity() {
                TermSpecificity::Anywhere => "none",
                other => other.name(),
            };
            assert_eq!(term, row[5], "{context}");
            assert_eq!(modification.full_name(), row[6], "{context}");
            assert_eq!(
                modification.diff_mono_mass(),
                row[7].parse::<f64>().unwrap(),
                "{context}"
            );
            assert_eq!(
                modification.diff_average_mass(),
                row[8].parse::<f64>().unwrap(),
                "{context}"
            );
            assert_eq!(
                modification.diff_formula().to_source_string(),
                row[9],
                "{context}"
            );
            assert_eq!(
                modification.source_classification_name(),
                row[10],
                "{context}"
            );
            assert_eq!(losses(modification), oracle_losses(row[11]), "{context}");
            assert_eq!(modification.mono_mass(), 0.0, "{context}");
        }
    }
}

/// The compiled table and the provider read the same files. They agree on
/// every field but `classification`, where the table keeps each specificity's
/// own value and the source provider gives every site of a modification its
/// last specificity's value.
#[test]
fn provider_agrees_with_the_compiled_table() {
    let mut loaded = load("unimod.xml");
    loaded.extend(load("custom_mods.xml"));
    let table: Vec<_> = ModificationsDB::global()
        .entries()
        .iter()
        .filter(|entry| entry.obo_accession().is_none())
        .collect();
    assert_eq!(loaded.len(), table.len());
    let mut differing_classifications = 0;
    for (index, (xml, compiled)) in loaded.iter().zip(&table).enumerate() {
        assert_eq!(xml.full_id(), compiled.full_id(), "{index}");
        assert_eq!(xml.record_id(), compiled.record_id());
        assert_eq!(xml.origin(), compiled.origin());
        assert_eq!(xml.term_specificity(), compiled.term_specificity());
        assert_eq!(xml.full_name(), compiled.full_name());
        assert_eq!(xml.diff_formula(), compiled.diff_formula());
        assert_eq!(xml.diff_mono_mass(), compiled.diff_mono_mass());
        assert_eq!(xml.diff_average_mass(), compiled.diff_average_mass());
        assert_eq!(xml.is_hidden(), compiled.is_hidden());
        assert_eq!(xml.neutral_losses(), compiled.neutral_losses());
        if xml.classification() != compiled.classification() {
            differing_classifications += 1;
            // The provider's value is that of the modification's last site.
            let last = table
                .iter()
                .rposition(|entry| {
                    entry.record_id() == xml.record_id() && entry.name() == xml.name()
                })
                .unwrap();
            assert_eq!(xml.classification(), table[last].classification());
        }
    }
    assert!(differing_classifications > 0);
}

#[test]
fn handler_rules_on_small_documents() {
    let xml = r#"<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">
 <umod:modifications>
  <umod:mod title="Test" full_name="Test mod" record_id="7">
   <umod:specificity hidden="1" site="S" position="Anywhere" classification="First">
    <umod:NeutralLoss mono_mass="97.976896" avge_mass="97.9952" flag="false">
     <umod:element symbol="H" number="3"/>
     <umod:element symbol="O" number="4"/>
     <umod:element symbol="P" number="1"/>
    </umod:NeutralLoss>
    <umod:NeutralLoss mono_mass="0" avge_mass="0" flag="false"/>
   </umod:specificity>
   <umod:specificity hidden="0" site="N-term" position="Protein N-term" classification="Last"/>
   <umod:specificity hidden="0" site="K" position="Somewhere" classification="Last"/>
   <umod:delta mono_mass="1.5" avge_mass="1.25">
    <umod:element symbol="13C" number="2"/>
    <umod:element symbol="H" number="-1"/>
   </umod:delta>
   <umod:Ignore><umod:element symbol="C" number="9"/></umod:Ignore>
  </umod:mod>
 </umod:modifications>
</umod:unimod>"#;
    let load = UnimodXmlFile.load_from_reader(xml.as_bytes()).unwrap();
    assert_eq!(load.records.len(), 3);
    assert_eq!(load.warnings.len(), 1);
    let mods = records_to_modifications(load.records).unwrap();
    let ids: Vec<_> = mods.iter().map(|m| m.full_id().to_owned()).collect();
    assert_eq!(ids, ["Test (S)", "Test (Protein N-term)", "Test (K)"]);
    for m in &mods {
        assert_eq!(m.diff_formula().to_source_string(), "(13)C2H-1");
        assert_eq!(m.diff_mono_mass(), 1.5);
        assert_eq!(m.classification(), "Last");
    }
    assert!(mods[0].is_hidden());
    assert_eq!(mods[0].neutral_losses().len(), 1);
    assert_eq!(
        mods[0].neutral_losses()[0].formula().to_source_string(),
        "H3O4P1"
    );
    assert!(mods[1].neutral_losses().is_empty());
    assert_eq!(mods[1].origin(), None);

    let missing = r#"<umod:mod title="T" record_id="1"><umod:specificity site="A"/></umod:mod>"#;
    assert!(UnimodXmlFile.load_from_reader(missing.as_bytes()).is_err());
    let unbalanced = r#"<umod:mod title="T" record_id="1"></umod:delta>"#;
    assert!(
        UnimodXmlFile
            .load_from_reader(unbalanced.as_bytes())
            .is_err()
    );
    let bad_id = r#"<umod:mod title="T" record_id="-1"/>"#;
    assert!(UnimodXmlFile.load_from_reader(bad_id.as_bytes()).is_err());
    assert!(
        UnimodXmlDataProvider::new("/nonexistent/unimod.xml")
            .load_modifications()
            .is_err()
    );
}

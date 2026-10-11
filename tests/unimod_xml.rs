// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/UnimodXMLFile.h`, `FORMAT/HANDLERS/UnimodXMLHandler.h` and
//! `CHEMISTRY/UnimodXMLDataProvider.h`.
//!
//! The class-test expectations come from `UnimodXMLFile_test.cpp`, run over
//! `resources/modifications/unimod.xml`, which is byte-identical to the
//! pinned `share/OpenMS/CHEMISTRY/unimod.xml` the class test loads (sha256
//! `22231d32...91091` for both). The registry cross-check compares the
//! runtime reader with the table `tools/generate_modifications.py` derives
//! from the same two XML files independently, in Python. See
//! `tests/data/unimod_xml_provenance.json`.
//!
//! The tier-1 tests (`provider_matches_the_release_build_record_for_record`,
//! `provider_agrees_with_the_compiled_table`,
//! `handler_rules_on_small_documents`) were group A's, written for its own
//! reader and replayed here against this one: `tests/data/unimod_provider_release.tsv`
//! is `UnimodXMLDataProvider::loadModifications` of the executed Release
//! build (`../oracle/A5a-unimod`, ibminode06, Linux x86_64 Release
//! bc9cc12/c19e494/174b576, run twice, identical) on the installed
//! `unimod.xml` (2,791 records) and `custom_mods.xml` (244).

#![cfg(feature = "kb-xml")]

use openms::Error;
use openms::chemistry::{EmpiricalFormula, ModificationsDB, ResidueModification, TermSpecificity};
use openms::format::unimod_xml::{
    self, ReadLimits, ReadOptions, SpecificityAttributes, UnimodXMLDataProvider, read_with_warnings,
};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn resource(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("resources/modifications")
        .join(name)
}

fn unimod() -> Vec<ResidueModification> {
    unimod_xml::load(resource("unimod.xml")).unwrap()
}

/// OpenMS `TEST_REAL_SIMILAR` with the class-test default tolerance.
fn similar(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 1e-5 * expected.abs().max(actual.abs()) || actual == expected
}

/// `UnimodXMLFile_test.cpp`, section `load`: more than one modification.
#[test]
fn load_reads_the_pinned_unimod_file() {
    let modifications = unimod();
    assert!(modifications.len() > 1);
    // One record per specificity: the embedded table holds 3035 records for
    // 1543 UniMod and 14 custom modifications, and the custom file has 22
    // specificities (tools/generate_modifications.py, provenance.json).
    let custom = unimod_xml::load(resource("custom_mods.xml")).unwrap();
    assert_eq!(modifications.len() + custom.len(), 3035);
    let titles: std::collections::BTreeSet<_> =
        modifications.iter().map(|m| m.record_id()).collect();
    assert_eq!(titles.len(), 1543);
}

/// `UnimodXMLFile_test.cpp`, section `[EXTRA] diff formula and diff mono
/// mass agree for every entry` (OpenMS issue 10030): only `<delta>` elements
/// make up the formula, so it weighs what `mono_mass` says, within 1e-3 Da.
#[test]
fn diff_formula_and_diff_mono_mass_agree_for_every_entry() {
    let modifications = unimod();
    let mismatches: Vec<_> = modifications
        .iter()
        .filter(|m| (m.diff_formula().mono_mass() - m.diff_mono_mass()).abs() > 1e-3)
        .map(|m| m.full_id().to_owned())
        .collect();
    assert_eq!(mismatches, Vec::<String>::new());

    let mut seen = 0;
    for modification in &modifications {
        if modification.name() == "ICAT-D" {
            seen += 1;
            assert_eq!(
                modification.diff_formula(),
                &EmpiricalFormula::parse("C20H34N4O5S").unwrap()
            );
            assert!(similar(modification.diff_mono_mass(), 442.224991));
        } else if modification.name() == "ICAT-D:2H(8)" {
            seen += 1;
            assert_eq!(
                modification.diff_formula(),
                &EmpiricalFormula::parse("C20H26(2)H8N4O5S").unwrap()
            );
            assert!(similar(modification.diff_mono_mass(), 450.275205));
        }
    }
    assert!(seen >= 2, "both ICAT-D entries are present");
}

/// `UnimodXMLFile_test.cpp`, section `[EXTRA] neutral loss formulas and
/// neutral loss masses agree for every entry`. The formula and mass vectors
/// cannot differ in length here, because each neutral loss is one value
/// holding all three; the mass checks carry over with the class test's
/// tolerances (1e-3 Da mono, 0.05 Da average).
#[test]
fn neutral_loss_formulas_and_masses_agree_for_every_entry() {
    let modifications = unimod();
    let with_neutral_loss = modifications
        .iter()
        .filter(|m| !m.neutral_losses().is_empty())
        .count();
    assert!(with_neutral_loss > 0);
    let mut mass_mismatches = Vec::new();
    for modification in &modifications {
        for loss in modification.neutral_losses() {
            if (loss.formula().mono_mass() - loss.mono_mass()).abs() > 1e-3
                || (loss.formula().average_mass() - loss.average_mass()).abs() > 0.05
            {
                mass_mismatches.push(modification.full_id().to_owned());
            }
        }
    }
    assert_eq!(mass_mismatches, Vec::<String>::new());
}

/// The runtime reader, keeping each specificity's attributes, reproduces the
/// registry the default build embeds, record for record and field for field
/// (`ResidueModification`'s equality compares every stored field).
#[test]
fn runtime_reader_reproduces_the_embedded_registry() {
    let mut records = unimod();
    records.extend(unimod_xml::load(resource("custom_mods.xml")).unwrap());
    let global = ModificationsDB::global().entries();
    assert!(
        global.len() > records.len(),
        "XLMOD follows the Unimod records"
    );
    for (index, (read, embedded)) in records.iter().zip(global).enumerate() {
        assert_eq!(read, embedded.as_ref(), "record {index}");
    }
    let db = ModificationsDB::from_records(records.clone()).unwrap();
    assert_eq!(db.len(), records.len());
    let oxidation = db.get_modification("Oxidation", Some('M'), None).unwrap();
    assert_eq!(oxidation.full_id(), "Oxidation (M)");
}

/// Source mode: `hidden` is dropped, and the last specificity's
/// classification, canonicalised by `setSourceClassification`, is every
/// specificity's (`UnimodXMLHandler.cpp:60-63`).
#[test]
fn source_mode_reproduces_the_source_classification_rule() {
    let retained = unimod();
    let source =
        unimod_xml::load_with_options(resource("unimod.xml"), &ReadOptions::source()).unwrap();
    assert_eq!(source.len(), retained.len());
    let mut last: BTreeMap<Option<u32>, &str> = BTreeMap::new();
    for record in &retained {
        last.insert(record.record_id(), record.classification());
    }
    let mut differs = 0;
    for (source, retained) in source.iter().zip(&retained) {
        assert!(!source.is_hidden());
        assert_eq!(source.full_id(), retained.full_id());
        assert_eq!(source.diff_formula(), retained.diff_formula());
        assert_eq!(source.neutral_losses(), retained.neutral_losses());
        // Unimod spells every classification but one exactly as the source's
        // getSourceClassificationName does, so canonicalisation is the
        // identity; "Synth. pep. protect. gp." is not in the source's table
        // (ResidueModification.cpp:306-387) and becomes UNKNOWN, named "".
        let expected = match last[&retained.record_id()] {
            "Synth. pep. protect. gp." => "",
            other => other,
        };
        assert_eq!(source.classification(), expected);
        if source.classification() != retained.classification() {
            differs += 1;
        }
    }
    // Acetyl alone has specificities classified "Post-translational",
    // "Chemical derivative", "Multiple" and "Artefact".
    assert!(differs > 0);
}

/// `UnimodXMLDataProvider::loadModifications`: the file's records, each with
/// its full identifier.
#[test]
fn data_provider_loads_the_file_records() {
    let provider = UnimodXMLDataProvider::new(resource("custom_mods.xml"));
    assert_eq!(provider.filename(), resource("custom_mods.xml"));
    assert_eq!(
        provider.options().specificity_attributes,
        SpecificityAttributes::Retain
    );
    let records = provider.load_modifications().unwrap();
    assert_eq!(
        records,
        unimod_xml::load(resource("custom_mods.xml")).unwrap()
    );
    assert!(records.iter().all(|m| !m.full_id().is_empty()));

    let source =
        UnimodXMLDataProvider::with_options(resource("custom_mods.xml"), ReadOptions::source());
    assert!(
        source
            .load_modifications()
            .unwrap()
            .iter()
            .all(|m| !m.is_hidden())
    );

    let missing = UnimodXMLDataProvider::new(resource("no-such-unimod.xml"));
    assert!(matches!(missing.load_modifications(), Err(Error::Io(_))));
}

#[test]
fn limits_refuse_instead_of_truncating() {
    let options = ReadOptions {
        limits: ReadLimits {
            max_records: 100,
            ..ReadLimits::default()
        },
        ..ReadOptions::default()
    };
    assert!(matches!(
        unimod_xml::load_with_options(resource("unimod.xml"), &options),
        Err(Error::InvalidValue(_))
    ));
    let options = ReadOptions {
        limits: ReadLimits {
            max_input_bytes: 1000,
            ..ReadLimits::default()
        },
        ..ReadOptions::default()
    };
    assert!(matches!(
        unimod_xml::load_with_options(resource("unimod.xml"), &options),
        Err(Error::InvalidValue(_))
    ));
}

const ORACLE: &str = include_str!("data/unimod_provider_release.tsv");

/// The source's records: `hidden` false and the last specificity's
/// classification at every site.
fn load_source(file: &str) -> Vec<ResidueModification> {
    UnimodXMLDataProvider::with_options(resource(file), ReadOptions::source())
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
        let loaded = load_source(file);
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

/// The compiled table and the reader read the same files. In the default
/// (`Retain`) mode they agree on every field; in `Source` mode, which the
/// Release build's records follow, they differ in the classification, where
/// the table keeps each specificity's own value and the source gives every
/// site of a modification its last specificity's value.
#[test]
fn provider_agrees_with_the_compiled_table() {
    let mut retained = unimod();
    retained.extend(unimod_xml::load(resource("custom_mods.xml")).unwrap());
    let mut source = load_source("unimod.xml");
    source.extend(load_source("custom_mods.xml"));
    let table: Vec<_> = ModificationsDB::global()
        .entries()
        .iter()
        .filter(|entry| entry.obo_accession().is_none())
        .collect();
    assert_eq!(retained.len(), table.len());
    assert_eq!(source.len(), table.len());
    for (index, (xml, compiled)) in retained.iter().zip(&table).enumerate() {
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
        assert_eq!(xml.classification(), compiled.classification());
    }
    let mut differing_classifications = 0;
    for (index, (xml, compiled)) in source.iter().zip(&table).enumerate() {
        assert_eq!(xml.full_id(), compiled.full_id(), "{index}");
        if xml.classification() != compiled.classification() {
            differing_classifications += 1;
            // The source's value is that of the modification's last site,
            // canonicalised by `getSourceClassificationName`.
            let last = table
                .iter()
                .rposition(|entry| {
                    entry.record_id() == xml.record_id() && entry.name() == xml.name()
                })
                .unwrap();
            assert_eq!(
                xml.classification(),
                table[last].source_classification_name()
            );
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
    // The source's records: every site carries the last classification.
    let (mods, warnings) = read_with_warnings(xml.as_bytes(), &ReadOptions::source()).unwrap();
    assert_eq!(mods.len(), 3);
    assert_eq!(warnings.len(), 1);
    let ids: Vec<_> = mods.iter().map(|m| m.full_id().to_owned()).collect();
    assert_eq!(ids, ["Test (S)", "Test (Protein N-term)", "Test (K)"]);
    for m in &mods {
        assert_eq!(m.diff_formula().to_source_string(), "(13)C2H-1");
        assert_eq!(m.diff_mono_mass(), 1.5);
        // "Last" is no name `getSourceClassificationName` knows.
        assert_eq!(m.classification(), "");
        assert_eq!(m.source_classification_name(), "");
    }
    assert_eq!(mods[0].neutral_losses().len(), 1);
    assert_eq!(
        mods[0].neutral_losses()[0].formula().to_source_string(),
        "H3O4P1"
    );
    assert!(mods[1].neutral_losses().is_empty());
    assert_eq!(mods[1].origin(), None);
    // Retained: each site keeps its own hidden flag and classification.
    let (kept, _) = read_with_warnings(xml.as_bytes(), &ReadOptions::default()).unwrap();
    assert!(kept[0].is_hidden());
    assert!(!kept[1].is_hidden());
    assert_eq!(kept[0].classification(), "First");
    assert_eq!(kept[1].classification(), "Last");

    let missing = r#"<umod:mod title="T" record_id="1"><umod:specificity site="A"/></umod:mod>"#;
    assert!(unimod_xml::read(missing.as_bytes()).is_err());
    let unbalanced = r#"<umod:mod title="T" record_id="1"></umod:delta>"#;
    assert!(unimod_xml::read(unbalanced.as_bytes()).is_err());
    let bad_id = r#"<umod:mod title="T" record_id="-1"/>"#;
    assert!(unimod_xml::read(bad_id.as_bytes()).is_err());
    assert!(
        UnimodXMLDataProvider::new("/nonexistent/unimod.xml")
            .load_modifications()
            .is_err()
    );
}

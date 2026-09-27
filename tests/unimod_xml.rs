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

#![cfg(feature = "kb-xml")]

use openms::Error;
use openms::chemistry::{EmpiricalFormula, ModificationsDB, ResidueModification};
use openms::format::unimod_xml::{
    self, ReadLimits, ReadOptions, SpecificityAttributes, UnimodXMLDataProvider,
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
    assert!(global.len() > records.len(), "XLMOD follows the Unimod records");
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
/// specificity's (`UnimodXMLHandler.cpp:57-60`).
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

    let source = UnimodXMLDataProvider::with_options(resource("custom_mods.xml"), ReadOptions::source());
    assert!(source.load_modifications().unwrap().iter().all(|m| !m.is_hidden()));

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

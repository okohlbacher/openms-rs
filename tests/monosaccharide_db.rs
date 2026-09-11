// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::MonosaccharideDB;

#[test]
fn source_public_lookup_sections_singleton_and_unknown_behavior() {
    let db = MonosaccharideDB::global();
    assert!(std::ptr::eq(db, MonosaccharideDB::global()));
    for name in ["Hex", "HexNAc", "Fuc", "Neu5Ac"] {
        assert!(db.has_symbol(name));
        assert!(db.get(name).is_some());
    }
    assert!(!db.has_symbol("NotAGlycanSymbol"));
    assert!(db.get("NotAGlycanSymbol").is_none());
    assert!(db.get_or_error("NotAGlycanSymbol").is_err());
    assert!(std::ptr::eq(
        db.get("HexNAc").unwrap(),
        db.get_or_error("HexNAc").unwrap()
    ));
    assert!(!db.is_empty());
    assert_eq!(db.len(), db.all_symbols().len());
    assert_eq!(db.len(), 24);
}

#[test]
fn all_source_data_mass_literals_formulas_and_primary_symbol_order() {
    // Independent transcription of all 24 JSON scientific values, not a formula
    // recomputation or an assertion against the generated Rust projection.
    let expected = [
        ("Dec", 282.095082153_f64, "C10H18O9"),
        ("Fuc", 146.057908799, "C6H10O4"),
        ("Hep", 192.063388102, "C7H12O6"),
        ("Hex", 162.052823418, "C6H10O5"),
        ("HexN", 161.068807836, "C6H11N1O4"),
        ("HexNAc", 203.07937252, "C8H13N1O5"),
        ("HexNAc(S)", 283.036187378, "C8H13N1O8S1"),
        ("HexNS", 241.025622694, "C6H11N1O7S1"),
        ("HexP", 242.019153939, "C6H11O8P1"),
        ("HexS", 242.009638277, "C6H10O8S1"),
        ("Neu", 249.084851823, "C9H15N1O7"),
        ("Neu5Ac", 291.095416506, "C11H17N1O8"),
        ("Neu5Gc", 307.090331126, "C11H17N1O9"),
        ("Non", 252.08451747, "C9H16O8"),
        ("Oct", 222.073952786, "C8H14O7"),
        ("Pen", 132.042258735, "C5H8O4"),
        ("Sug", 42.0105646837, "C2H2O1"),
        ("Tet", 102.031694051, "C4H6O3"),
        ("Tri", 72.0211293674, "C3H4O2"),
        ("a-Hex", 176.032087974, "C6H8O6"),
        ("d-Hex", 146.057908799, "C6H10O4"),
        ("en,a-Hex", 158.02152329, "C6H6O5"),
        ("phosphate", 79.9663305208, "H1O3P1"),
        ("sulfate", 79.9568148587, "H0O3S1"),
    ];
    let db = MonosaccharideDB::global();
    assert_eq!(
        db.all_symbols(),
        expected.iter().map(|row| row.0).collect::<Vec<_>>()
    );
    for (symbol, mass, formula) in expected {
        let record = db.get_or_error(symbol).unwrap();
        assert_eq!(record.symbol, symbol);
        assert_eq!(record.mass.to_bits(), mass.to_bits(), "{symbol}");
        assert_eq!(record.formula, formula, "{symbol}");
    }
    assert_eq!(
        db.get("Hex").unwrap().name,
        "A generic monosaccharide with 6 backbone carbons"
    );
    assert_eq!(db.get("Neu5Ac").unwrap().name, "Neu5Ac");
}

#[test]
fn every_alias_returns_the_same_record_and_does_not_add_a_primary_symbol() {
    let db = MonosaccharideDB::global();
    let expected = [
        ("S", "sulfate"),
        ("P", "phosphate"),
        ("dHex", "d-Hex"),
        ("Fucose", "Fuc"),
        ("en,aHex", "en,a-Hex"),
        ("enHexA", "en,a-Hex"),
        ("aHex", "a-Hex"),
        ("HexA", "a-Hex"),
        ("Neuraminic acid", "Neu"),
        ("Sialic Acid", "Neu5Ac"),
        ("NeuAc", "Neu5Ac"),
        ("NeuGc", "Neu5Gc"),
    ];
    let mut synonyms = 0;
    for symbol in db.all_symbols() {
        synonyms += db.get(symbol).unwrap().synonyms.len();
    }
    assert_eq!(synonyms, expected.len());
    for (alias, primary) in expected {
        assert!(db.has_symbol(alias));
        assert!(std::ptr::eq(
            db.get(alias).unwrap(),
            db.get(primary).unwrap()
        ));
        assert!(db.get(primary).unwrap().synonyms.iter().any(|s| s == alias));
        assert!(!db.all_symbols().contains(&alias));
    }
}

#[test]
fn exact_queries_and_owned_record_clones_do_not_modify_the_registry() {
    let db = MonosaccharideDB::global();
    for unknown in ["", "hex", " Hex", "Hex ", "Hex\n", "Ｈex", "Neu5Ac\0"] {
        assert!(!db.has_symbol(unknown));
        assert!(db.get(unknown).is_none());
    }
    let long = "Neu5Ac".repeat(200_000);
    assert!(db.get(&long).is_none());
    assert_eq!(
        db.get_or_error(&long).unwrap_err().to_string(),
        "invalid value: unknown monosaccharide symbol"
    );
    let mut copy = db.get("NeuAc").unwrap().clone();
    copy.symbol = "local".into();
    copy.mass = -1.;
    copy.formula.clear();
    copy.synonyms.push("caller alias".into());
    assert_eq!(db.get("NeuAc").unwrap().symbol, "Neu5Ac");
    assert_eq!(db.get("NeuAc").unwrap().synonyms.len(), 2);
    assert_eq!(db.get("NeuAc").unwrap().mass, 291.095416506);
    assert_eq!(db.get("NeuAc").unwrap().formula, "C11H17N1O8");
    assert!(!db.has_symbol("caller alias"));
}

#[test]
fn concurrent_global_access_retains_one_immutable_identity() {
    let expected = MonosaccharideDB::global() as *const MonosaccharideDB as usize;
    let workers: Vec<_> = (0..8)
        .map(|_| {
            std::thread::spawn(|| {
                let db = MonosaccharideDB::global();
                assert_eq!(db.get("Fucose").unwrap().symbol, "Fuc");
                db as *const MonosaccharideDB as usize
            })
        })
        .collect();
    for worker in workers {
        assert_eq!(worker.join().unwrap(), expected);
    }
}

#[cfg(feature = "rna-json")]
#[test]
fn independent_json_decoder_matches_every_original_field_and_mass_bit() {
    let raw: serde_json::Value = serde_json::from_str(include_str!(
        "../resources/monosaccharides/monosaccharides.json"
    ))
    .unwrap();
    let rows = raw["monosaccharides"].as_object().unwrap();
    let db = MonosaccharideDB::global();
    assert_eq!(db.len(), rows.len());
    for (symbol, row) in rows {
        let record = db.get_or_error(symbol).unwrap();
        assert_eq!(record.symbol, *symbol);
        assert_eq!(
            record.mass.to_bits(),
            row["mass"].as_f64().unwrap().to_bits()
        );
        assert_eq!(record.name, row["name"].as_str().unwrap());
        assert_eq!(record.formula, row["formula"].as_str().unwrap());
        assert_eq!(
            record.synonyms,
            row["synonyms"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap())
                .collect::<Vec<_>>()
        );
    }
}

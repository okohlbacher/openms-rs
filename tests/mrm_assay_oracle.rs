// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential tests (evidence tier 1) of `MRMAssay` and `MRMDecoy` against
//! the executed C++ Release build.
//!
//! The oracle `../oracle/O1c-assays/driver.cpp` (not committed; hashes in
//! `tests/data/mrm_assay_provenance.json`) ran on ibminode06 against the
//! Linux x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576`,
//! twice per case with byte-identical output. It loads the class-test TraML
//! inputs with `TraMLFile`, which is not ported, and prints a canonical dump
//! of every experiment: `input_*` files dump the loaded inputs, from which
//! this test rebuilds each experiment, and the other files dump the result of
//! the same call the class test (or a documented extra case) makes.
//!
//! Compared per record: protein id and sequence; peptide id, sequence,
//! charge, first RT, group label, protein references and modifications;
//! compound id, charge and RT; transition id, references, precursor and
//! product m/z, library intensity, decoy type, the three flags, product
//! charge, every interpretation (ordinal, rank, ion type, CV-term accessions
//! and values) and every meta value. Numbers are compared exactly.

#![allow(clippy::field_reassign_with_default)]

use openms::analysis::openswath::mrm_assay::{MRMAssay, UisOptions};
use openms::analysis::openswath::mrm_decoy::{DecoyMethod, DecoyOptions, MRMDecoy};
use openms::analysis::openswath::mrm_ion_series::{FormulaLossParsing, IonSeriesOptions};
use openms::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{
    Compound, Interpretation, Peptide, PeptideCompound, PeptideModification, Protein,
    RetentionTime, TraMLProduct, UnresolvedModificationPolicy,
};
use openms::chemistry::PeptideFragmentType;
use openms::metadata::{CVTerm, MetaValue, MetaValueData};
use std::path::PathBuf;

fn data(name: &str) -> String {
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "tests/data/mrm_assay/oracle",
        name,
    ]
    .iter()
    .collect();
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some(other) => out.push(other),
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn num(text: &str) -> f64 {
    text.parse()
        .unwrap_or_else(|_| panic!("not a number: {text}"))
}

fn opt_int(text: &str) -> Option<i32> {
    (text != "NA").then(|| text.parse().unwrap())
}

fn opt_num(text: &str) -> Option<f64> {
    (text != "NA").then(|| num(text))
}

/// A typed value as the oracle prints it: `s:`, `i:`, `d:` or `e:`.
#[derive(Debug, PartialEq)]
enum Value {
    Text(String),
    Int(i64),
    Float(f64),
    Empty,
}

fn parse_value(text: &str) -> Value {
    let (kind, rest) = text.split_at(2);
    match kind {
        "s:" | "o:" => Value::Text(unescape(rest)),
        "i:" => Value::Int(rest.parse().unwrap()),
        "d:" => Value::Float(num(rest)),
        "e:" => Value::Empty,
        _ => panic!("bad value {text}"),
    }
}

fn value_of(meta: &MetaValue) -> Value {
    match meta.data() {
        MetaValueData::String(s) => Value::Text(s.clone()),
        MetaValueData::Integer(i) => Value::Int(*i),
        MetaValueData::Float(f) => Value::Float(*f),
        MetaValueData::Empty => Value::Empty,
        other => Value::Text(format!("{other:?}")),
    }
}

fn meta_value(value: &Value) -> MetaValue {
    match value {
        Value::Text(s) => MetaValue::from(s.as_str()),
        Value::Int(i) => MetaValue::from(*i),
        Value::Float(f) => MetaValue::try_from(*f).unwrap(),
        Value::Empty => MetaValue::default(),
    }
}

/// `key=value;key=value` pairs; keys and accessions never contain `=` or `;`.
fn parse_pairs(text: &str) -> Vec<(String, Value)> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|pair| {
            let (key, value) = pair.split_once('=').unwrap();
            (unescape(key), parse_value(value))
        })
        .collect()
}

fn iontype(text: &str) -> PeptideFragmentType {
    match text {
        "a" => PeptideFragmentType::AIon,
        "b" => PeptideFragmentType::BIon,
        "c" => PeptideFragmentType::CIon,
        "x" => PeptideFragmentType::XIon,
        "y" => PeptideFragmentType::YIon,
        "z" => PeptideFragmentType::ZIon,
        "NonIdentified" => PeptideFragmentType::NonIdentified,
        "Unannotated" => PeptideFragmentType::Unannotated,
        other => panic!("unsupported ion type {other}"),
    }
}

fn iontype_name(t: PeptideFragmentType) -> String {
    match t {
        PeptideFragmentType::AIon => "a".into(),
        PeptideFragmentType::BIon => "b".into(),
        PeptideFragmentType::CIon => "c".into(),
        PeptideFragmentType::XIon => "x".into(),
        PeptideFragmentType::YIon => "y".into(),
        PeptideFragmentType::ZIon => "z".into(),
        PeptideFragmentType::NonIdentified => "NonIdentified".into(),
        PeptideFragmentType::Unannotated => "Unannotated".into(),
        other => format!("{other:?}"),
    }
}

/// One comparable record, with numbers kept as numbers.
#[derive(Debug, PartialEq)]
enum Record {
    Protein(String, String),
    Peptide {
        id: String,
        sequence: String,
        charge: Option<i32>,
        rt: Option<f64>,
        group: String,
        refs: Vec<String>,
        mods: Vec<(i32, i32, f64, f64)>,
    },
    Compound(String, Option<i32>, Option<f64>),
    Transition {
        id: String,
        peptide_ref: String,
        compound_ref: String,
        precursor_mz: f64,
        product_mz: f64,
        intensity: f64,
        decoy: i32,
        flags: String,
        product_charge: Option<i32>,
        interpretations: Vec<(i32, i32, String, Vec<(String, Value)>)>,
        metas: Vec<(String, Value)>,
    },
}

fn parse_dump(text: &str) -> Vec<Record> {
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "protein" => Record::Protein(unescape(f[1]), unescape(f[2])),
                "peptide" => Record::Peptide {
                    id: unescape(f[1]),
                    sequence: unescape(f[2]),
                    charge: opt_int(f[3]),
                    rt: opt_num(f[4]),
                    group: unescape(f[5]),
                    refs: if f[6].is_empty() {
                        Vec::new()
                    } else {
                        f[6].split(',').map(unescape).collect()
                    },
                    mods: if f[7].is_empty() {
                        Vec::new()
                    } else {
                        f[7].split(';')
                            .map(|m| {
                                let p: Vec<&str> = m.split(':').collect();
                                (
                                    p[0].parse().unwrap(),
                                    p[1].parse().unwrap(),
                                    num(p[2]),
                                    num(p[3]),
                                )
                            })
                            .collect()
                    },
                },
                "compound" => Record::Compound(unescape(f[1]), opt_int(f[2]), opt_num(f[3])),
                "transition" => Record::Transition {
                    id: unescape(f[1]),
                    peptide_ref: unescape(f[2]),
                    compound_ref: unescape(f[3]),
                    precursor_mz: num(f[4]),
                    product_mz: num(f[5]),
                    intensity: num(f[6]),
                    decoy: f[7].parse().unwrap(),
                    flags: f[8].into(),
                    product_charge: opt_int(f[9]),
                    interpretations: if f[10].is_empty() {
                        Vec::new()
                    } else {
                        f[10]
                            .split('|')
                            .map(|i| {
                                let p: Vec<&str> = i.splitn(4, ':').collect();
                                (
                                    p[0].parse().unwrap(),
                                    p[1].parse().unwrap(),
                                    p[2].to_owned(),
                                    parse_pairs(p[3]),
                                )
                            })
                            .collect()
                    },
                    metas: parse_pairs(f.get(11).copied().unwrap_or_default()),
                },
                other => panic!("unexpected record {other}: {line}"),
            }
        })
        .collect()
}

fn rt_of(p: &PeptideCompound) -> Option<f64> {
    p.retention_time().ok()
}

fn records(exp: &TargetedExperiment) -> Vec<Record> {
    let mut out = Vec::new();
    for p in exp.proteins() {
        out.push(Record::Protein(p.id.clone(), p.sequence.clone()));
    }
    for p in exp.peptides() {
        out.push(Record::Peptide {
            id: p.id.clone(),
            sequence: p.sequence.clone(),
            charge: p.charge_state(),
            rt: rt_of(p),
            group: p.peptide_group_label().to_owned(),
            refs: p.protein_refs.clone(),
            mods: p
                .mods
                .iter()
                .map(|m| (m.location, m.unimod_id, m.mono_mass_delta, m.avg_mass_delta))
                .collect(),
        });
    }
    for c in exp.compounds() {
        out.push(Record::Compound(c.id.clone(), c.charge_state(), rt_of(c)));
    }
    for t in exp.transitions() {
        let terms = |i: &Interpretation| {
            let mut v: Vec<(String, Value)> = i
                .cv_terms
                .terms()
                .iter()
                .flat_map(|(acc, list)| list.iter().map(move |t| (acc.clone(), value_of(&t.value))))
                .collect();
            v.sort_by(|a, b| format!("{}={:?}", a.0, a.1).cmp(&format!("{}={:?}", b.0, b.1)));
            v
        };
        out.push(Record::Transition {
            id: t.name().to_owned(),
            peptide_ref: t.peptide_ref().to_owned(),
            compound_ref: t.compound_ref().to_owned(),
            precursor_mz: t.precursor_mz(),
            product_mz: t.product_mz(),
            intensity: t.library_intensity(),
            decoy: match t.decoy_transition_type() {
                DecoyTransitionType::Unknown => 0,
                DecoyTransitionType::Target => 1,
                DecoyTransitionType::Decoy => 2,
            },
            flags: format!(
                "{}{}{}",
                u8::from(t.is_detecting_transition()),
                u8::from(t.is_identifying_transition()),
                u8::from(t.is_quantifying_transition())
            ),
            product_charge: t.product().charge_state(),
            interpretations: t
                .product()
                .interpretation_list()
                .iter()
                .map(|i| {
                    (
                        i32::from(i.ordinal),
                        i32::from(i.rank),
                        iontype_name(i.iontype),
                        terms(i),
                    )
                })
                .collect(),
            metas: t
                .cv_terms
                .metadata
                .iter()
                .map(|(k, v)| (k.clone(), value_of(v)))
                .collect(),
        });
    }
    out
}

fn rebuild(dump: &str) -> TargetedExperiment {
    let mut exp = TargetedExperiment::new();
    for record in parse_dump(dump) {
        match record {
            Record::Protein(id, sequence) => {
                let mut p = Protein::default();
                p.id = id;
                p.sequence = sequence;
                exp.add_protein(p);
            }
            Record::Peptide {
                id,
                sequence,
                charge,
                rt,
                group,
                refs,
                mods,
            } => {
                let mut p = Peptide::default();
                p.id = id;
                p.sequence = sequence;
                if let Some(c) = charge {
                    p.set_charge_state(c);
                }
                if let Some(rt) = rt {
                    let mut r = RetentionTime::new();
                    r.set_rt(rt);
                    p.rts.push(r);
                }
                p.set_peptide_group_label(group);
                p.protein_refs = refs;
                for (location, unimod_id, mono, avg) in mods {
                    let mut m = PeptideModification::default();
                    m.location = location;
                    m.unimod_id = unimod_id;
                    m.mono_mass_delta = mono;
                    m.avg_mass_delta = avg;
                    p.mods.push(m);
                }
                exp.add_peptide(p);
            }
            Record::Compound(id, charge, rt) => {
                let mut c = Compound::default();
                c.id = id;
                if let Some(ch) = charge {
                    c.set_charge_state(ch);
                }
                if let Some(rt) = rt {
                    let mut r = RetentionTime::new();
                    r.set_rt(rt);
                    c.rts.push(r);
                }
                exp.add_compound(c);
            }
            Record::Transition {
                id,
                peptide_ref,
                compound_ref,
                precursor_mz,
                product_mz,
                intensity,
                decoy,
                flags,
                product_charge,
                interpretations,
                metas,
            } => {
                let mut t = ReactionMonitoringTransition::new();
                t.set_name(id);
                t.set_peptide_ref(peptide_ref);
                t.set_compound_ref(compound_ref);
                t.set_precursor_mz(precursor_mz);
                t.set_library_intensity(intensity);
                t.set_decoy_transition_type(match decoy {
                    1 => DecoyTransitionType::Target,
                    2 => DecoyTransitionType::Decoy,
                    _ => DecoyTransitionType::Unknown,
                });
                let flag = |i: usize| flags.as_bytes()[i] == b'1';
                t.set_detecting_transition(flag(0));
                t.set_identifying_transition(flag(1));
                t.set_quantifying_transition(flag(2));
                let mut product = TraMLProduct::new();
                if let Some(c) = product_charge {
                    product.set_charge_state(c);
                }
                for (ordinal, rank, kind, terms) in interpretations {
                    let mut i = Interpretation::default();
                    i.ordinal = u8::try_from(ordinal).unwrap();
                    i.rank = u8::try_from(rank).unwrap();
                    i.iontype = iontype(&kind);
                    for (accession, value) in terms {
                        let mut term = CVTerm::new(accession, "", "MS");
                        term.value = meta_value(&value);
                        i.cv_terms.add(term).unwrap();
                    }
                    product.add_interpretation(i);
                }
                t.set_product(product);
                // The product m/z lives on the product; set it after.
                t.set_product_mz(product_mz);
                for (key, value) in metas {
                    t.cv_terms.metadata.insert(key, meta_value(&value));
                }
                exp.add_transition(t);
            }
        }
    }
    exp
}

#[track_caller]
fn assert_matches(case: &str, actual: &TargetedExperiment) {
    let expected = parse_dump(&data(&format!("{case}.tsv")));
    let actual = records(actual);
    assert_eq!(actual.len(), expected.len(), "{case}: record count differs");
    for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
        assert_eq!(a, e, "{case}: record {i} differs");
    }
}

fn input(file: &str) -> TargetedExperiment {
    rebuild(&data(&format!("input_{file}.tsv")))
}

fn source_options<S: Into<String>>(
    types: impl IntoIterator<Item = S>,
    charges: impl IntoIterator<Item = usize>,
    specific: bool,
    unspecific: bool,
) -> IonSeriesOptions {
    let mut options = IonSeriesOptions::new(types, charges).with_losses(specific, unspecific);
    // The source skips unresolvable modifications with a warning.
    options.unresolved_modifications = UnresolvedModificationPolicy::SourceSkip;
    // Each oracle case runs in a fresh process, so the source's static
    // formula-loss cache starts empty per call here too.
    options.formula_loss = FormulaLossParsing::source_static();
    options
}

#[test]
fn inputs_round_trip_through_the_dump() {
    for file in [
        "MRMAssay_reannotateTransitions_input.TraML",
        "MRMAssay_restrictTransitions_input.TraML",
        "MRMAssay_detectingTransitions_input.TraML",
        "MRMAssay_detectingTransistionCompound_input.TraML",
        "MRMAssay_uisTransitions_input_1.TraML",
        "MRMAssay_uisTransitions_input_3.TraML",
        "MRMDecoyGenerator_input.TraML",
    ] {
        let exp = input(file);
        assert_matches(&format!("input_{file}"), &exp);
    }
}

// MRMAssay_test.cpp reannotateTransitions, both parameter sets.
#[test]
fn reannotate_transitions() {
    let mut exp = input("MRMAssay_reannotateTransitions_input.TraML");
    MRMAssay::reannotate_transitions(
        &mut exp,
        0.05,
        0.05,
        &source_options(["y"], [2], false, false),
    )
    .unwrap();
    assert_matches("reannotate_1", &exp);
    let mut exp = input("MRMAssay_reannotateTransitions_input.TraML");
    MRMAssay::reannotate_transitions(
        &mut exp,
        0.05,
        0.05,
        &source_options(["y", "b"], [2, 3], true, true),
    )
    .unwrap();
    assert_matches("reannotate_2", &exp);
}

// MRMAssay_test.cpp restrictTransitions, plus the class test's swath list.
#[test]
fn restrict_transitions() {
    let options = source_options(["y"], [1], false, false);
    let mut exp = input("MRMAssay_restrictTransitions_input.TraML");
    MRMAssay::restrict_transitions(&mut exp, 400.0, 2000.0, &[], &options).unwrap();
    assert_matches("restrict", &exp);
    let mut swathes = vec![(400.0, 425.0)];
    let mut lower = 424.0;
    while lower < 1174.5 {
        swathes.push((lower, lower + 26.0));
        lower += 25.0;
    }
    let mut exp = input("MRMAssay_restrictTransitions_input.TraML");
    MRMAssay::restrict_transitions(&mut exp, 400.0, 2000.0, &swathes, &options).unwrap();
    assert_matches("restrict_swathes", &exp);
}

// MRMAssay_test.cpp detectingTransitions.
#[test]
fn detecting_transitions() {
    let mut exp = input("MRMAssay_detectingTransitions_input.TraML");
    MRMAssay::detecting_transitions(&mut exp, 4, 6);
    assert_matches("detecting", &exp);
}

// MRMAssay_test.cpp filterMinMaxTransitionsCompound, plus
// filterUnreferencedDecoysCompound on the same input.
#[test]
fn compound_filters() {
    let mut exp = input("MRMAssay_detectingTransistionCompound_input.TraML");
    MRMAssay::filter_min_max_transitions_compound(&mut exp, 3, 6);
    assert_matches("compound_minmax", &exp);
    let mut exp = input("MRMAssay_detectingTransistionCompound_input.TraML");
    MRMAssay::filter_unreferenced_decoys_compound(&mut exp).unwrap();
    assert_matches("compound_unreferenced", &exp);
}

fn uis(file: &str, options: IonSeriesOptions, ms2: bool) -> TargetedExperiment {
    let mut exp = input(file);
    let mut uis = UisOptions::new(options, 0.05);
    uis.enable_ms2_precursors = ms2;
    uis.shuffle_seed = Some(42);
    MRMAssay::uis_transitions(&mut exp, &uis).unwrap();
    exp
}

// MRMAssay_test.cpp uisTransitions, all five parameter sets (seed 42).
#[test]
fn uis_transitions() {
    let input_1 = "MRMAssay_uisTransitions_input_1.TraML";
    let input_3 = "MRMAssay_uisTransitions_input_3.TraML";
    assert_matches(
        "uis_1",
        &uis(input_1, source_options(["y"], [2], true, false), false),
    );
    assert_matches(
        "uis_2",
        &uis(input_1, source_options(["y"], [2], true, true), false),
    );
    assert_matches(
        "uis_3",
        &uis(input_3, source_options(["b"], [3], true, true), false),
    );
    assert_matches(
        "uis_4",
        &uis(
            input_3,
            source_options(["y", "b"], [2, 3], true, true),
            false,
        ),
    );
    assert_matches(
        "uis_5",
        &uis(
            input_3,
            source_options(["y", "b"], [2, 3], true, true),
            true,
        ),
    );
}

fn decoys(method: DecoyMethod, fraction: f64, switch_kr: bool) -> TargetedExperiment {
    let exp = input("MRMDecoyGenerator_input.TraML");
    let mut options = DecoyOptions::new(source_options(["b", "y", "a"], 1..=5, true, false));
    options.method = method;
    options.aim_decoy_fraction = fraction;
    options.switch_kr = switch_kr;
    options.max_attempts = 5;
    options.identity_threshold = 0.7;
    options.precursor_mz_shift = 0.1;
    options.product_mz_shift = 20.0;
    options.product_mz_threshold = 0.8;
    let mut dec = TargetedExperiment::new();
    MRMDecoy::new()
        .generate_decoys(&exp, &mut dec, &options)
        .unwrap();
    dec
}

// MRMDecoy_test.cpp generateDecoys (pseudo-reverse), plus reverse, a 0.5
// fraction (the default-seeded RandomShuffler) and switchKR.
#[test]
fn generate_decoys() {
    assert_matches(
        "decoys_pseudo",
        &decoys(DecoyMethod::PseudoReverse, 1.0, false),
    );
    assert_matches("decoys_reverse", &decoys(DecoyMethod::Reverse, 1.0, false));
    assert_matches(
        "decoys_fraction",
        &decoys(DecoyMethod::PseudoReverse, 0.5, false),
    );
    assert_matches(
        "decoys_switchkr",
        &decoys(DecoyMethod::PseudoReverse, 1.0, true),
    );
}

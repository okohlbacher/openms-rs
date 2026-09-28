// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/PEFFFile.h`: the sequence-producing `PEFFEntry` members
//! (`getSequence`, `getModifiedSequence`, `getVariantSequences`,
//! `getProcessedSequence`, `digestWithVariants`, `generatePeptides`,
//! `enumeratePEFFModifications_`).
//!
//! Every expectation is transcribed from `PEFFFile_test.cpp` at the pinned
//! core revision `bc9cc12`; each test quotes its section. The proteoform
//! tables are that file's `REF_TABLE`, `K16R_TABLE`, `S14A_TABLE` and
//! `CV_TABLE`, row for row. The source's `ModificationsDB` loads PSI-MOD by
//! default; the global registry here holds Unimod only, so the two PSI-MOD
//! assertions resolve against a registry extended with the one term they
//! use, authored below from its identifier, residue and Unimod link.

use openms::chemistry::{
    AASequence, ModificationsDB, OboReadOptions, ProteaseDigestion, TermSpecificity,
};
use openms::format::peff::{
    PEFFEntry, PEFFFile, PEFFModification, PEFFProcessedRegion, PEFFVariantComplex,
    PEFFVariantSimple, PeptideGenerationOptions, VariantDigestOptions,
};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn trypsin() -> ProteaseDigestion {
    let mut digestor = ProteaseDigestion::default();
    digestor.set_enzyme("Trypsin").unwrap();
    digestor.missed_cleavages = 0;
    digestor
}

fn variant(position: usize, residue: char, tag: &str) -> PEFFVariantSimple {
    PEFFVariantSimple::new(
        position,
        residue,
        tag,
        openms::format::peff::ANNOTATION_ID_NOT_SET,
    )
}

fn digest_options(
    min_length: usize,
    max_length: usize,
    include_reference: bool,
    include_variants: bool,
    include_modifications: bool,
) -> VariantDigestOptions {
    VariantDigestOptions {
        min_length,
        max_length,
        include_reference,
        include_variants,
        include_modifications,
        ..VariantDigestOptions::default()
    }
}

#[allow(clippy::too_many_arguments)]
fn generation_options(
    fixed: &[&str],
    variable: &[&str],
    max_variable: usize,
    min_length: usize,
    max_length: usize,
    include_reference: bool,
    include_variants: bool,
    include_modifications: bool,
) -> PeptideGenerationOptions {
    PeptideGenerationOptions {
        fixed_mods: fixed.iter().map(|s| (*s).to_owned()).collect(),
        variable_mods: variable.iter().map(|s| (*s).to_owned()).collect(),
        max_variable_mods_per_peptide: max_variable,
        min_length,
        max_length,
        include_reference,
        include_peff_variants: include_variants,
        include_peff_modifications: include_modifications,
        ..PeptideGenerationOptions::default()
    }
}

fn modified_residues(sequence: &AASequence) -> usize {
    (0..sequence.len())
        .filter(|&i| sequence.residue_modification(i).unwrap().is_some())
        .count()
}

/// The global registry plus PSI-MOD `MOD:00046` (O-phospho-L-serine, origin
/// S, linked to Unimod 21), the one PSI-MOD term the source tests use.
fn with_psi_mod_00046() -> ModificationsDB {
    const TERM: &str = "format-version: 1.2\n\n[Term]\nid: MOD:00046\n\
                        name: O-phospho-L-serine\n\
                        def: \"Phosphorylated L-serine.\" [UniMod:21]\n\
                        property_value: Origin: \"S\" xsd:string\n\
                        property_value: TermSpec: \"none\" xsd:string\n";
    let mut db = ModificationsDB::global().clone();
    db.extend_obo(TERM.as_bytes(), &OboReadOptions::default())
        .unwrap();
    db
}

/// Section `[PEFFEntry] AASequence getSequence() const`.
#[test]
fn aa_sequence() {
    let entry = PEFFEntry {
        sequence: "PEPTIDE".into(),
        ..PEFFEntry::default()
    };
    let sequence = entry.aa_sequence().unwrap();
    assert_eq!(sequence.to_string(), "PEPTIDE");
    assert_eq!(sequence.len(), 7);
}

/// Section `[PEFFEntry] AASequence getModifiedSequence() const`.
#[test]
fn modified_sequence_without_modifications() {
    let entry = PEFFEntry {
        sequence: "PEPTIDE".into(),
        ..PEFFEntry::default()
    };
    assert_eq!(entry.modified_sequence().unwrap().len(), 7);
}

/// Section `[PEFFEntry] getVariantSequences`.
#[test]
fn variant_sequences() {
    let mut entry = PEFFEntry {
        sequence: "PEPTIDE".into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(3, 'A', "dbSNP:test"));
    entry.complex_variants.push(PEFFVariantComplex::new(
        2,
        4,
        "XXX",
        "ClinVar:123",
        openms::format::peff::ANNOTATION_ID_NOT_SET,
    ));

    let simple = entry.variant_sequences(false);
    assert_eq!(simple.sequences.len(), 1);
    assert_eq!(simple.descriptions.len(), 1);
    assert_eq!(simple.sequences[0].to_string(), "PEATIDE");
    assert!(simple.descriptions[0].contains("P3A"));

    let all = entry.variant_sequences(true);
    assert_eq!(all.sequences.len(), 2);
    assert_eq!(all.sequences[0].to_string(), "PEATIDE");
    assert_eq!(all.sequences[1].to_string(), "PXXXIDE");
    assert!(all.descriptions[1].contains("2-4>XXX"));
}

/// Section `[PEFFEntry] AASequence getProcessedSequence(const std::string&
/// region_type) const`.
#[test]
fn processed_sequence_after_a_signal_peptide() {
    let mut entry = PEFFEntry {
        sequence: "SIGNALPEPTIDERESTOFPROTEIN".into(),
        ..PEFFEntry::default()
    };
    entry.processed_regions.push(PEFFProcessedRegion::new(
        1,
        13,
        "PEFF:0001021",
        "signal peptide",
        "",
        openms::format::peff::ANNOTATION_ID_NOT_SET,
    ));
    assert_eq!(
        entry
            .processed_sequence("PEFF:0001021")
            .unwrap()
            .to_string(),
        "RESTOFPROTEIN"
    );
}

/// Section `[PEFFEntry] digestWithVariants`.
#[test]
fn digest_with_variants() {
    let mut entry = PEFFEntry {
        sequence: "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAK"
            .into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(5, 'A', "dbSNP:rs1"));
    entry.simple_variants.push(variant(10, 'M', "dbSNP:rs2"));
    entry.simple_variants.push(variant(20, 'F', "dbSNP:rs3"));
    entry.simple_variants.push(variant(30, '*', "dbSNP:rs4"));
    let digestor = trypsin();

    let all = entry
        .digest_with_variants(&digestor, &digest_options(4, 30, true, true, false))
        .unwrap();
    assert!(!all.sequences.is_empty());
    assert_eq!(all.descriptions.len(), all.sequences.len());
    let references = all.descriptions.iter().filter(|d| d.is_empty()).count();
    let variants = all.descriptions.len() - references;
    for sequence in &all.sequences {
        assert!(!sequence.to_string().contains('*'));
    }
    assert!(references > 0);
    assert!(variants > 0);

    let variants_only = entry
        .digest_with_variants(&digestor, &digest_options(4, 30, false, true, false))
        .unwrap();
    for description in &variants_only.descriptions {
        assert!(!description.is_empty());
    }
    assert!(variants_only.sequences.len() < all.sequences.len());

    let mut chymotrypsin = ProteaseDigestion::default();
    chymotrypsin.set_enzyme("Chymotrypsin").unwrap();
    chymotrypsin.missed_cleavages = 0;
    let chymo = entry
        .digest_with_variants(&chymotrypsin, &digest_options(4, 30, true, true, false))
        .unwrap();
    assert!(!chymo.sequences.is_empty());
    assert_ne!(chymo.sequences.len(), all.sequences.len());
}

/// Section `[PEFFEntry] digestWithVariants_combinatorics`.
#[test]
fn digest_with_variants_combinatorics() {
    let mut entry = PEFFEntry {
        sequence: "MKAEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(4, 'L', "var1"));
    entry.simple_variants.push(variant(6, 'L', "var2"));
    let peptides = entry
        .digest_with_variants(&trypsin(), &digest_options(2, 20, true, true, false))
        .unwrap();
    let strings: Vec<String> = peptides.sequences.iter().map(ToString::to_string).collect();
    let found = |s: &str| strings.iter().any(|x| x == s);
    assert!(found("AEPTIDER"));
    assert!(found("ALPTIDER"));
    assert!(found("AEPLIDER"));
    assert!(found("ALPLIDER"));
    let versions = strings
        .iter()
        .filter(|s| ["AEPTIDER", "ALPTIDER", "AEPLIDER", "ALPLIDER"].contains(&s.as_str()))
        .count();
    assert_eq!(versions, 4);
}

/// Section `[PEFFEntry] digestWithVariants_modifications`.
#[test]
fn digest_with_variants_and_modifications() {
    let mut entry = PEFFEntry {
        sequence: "MKAMEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry
        .modifications
        .push(PEFFModification::new(4, "UNIMOD:35", "Oxidation"));
    let digestor = trypsin();

    let modifications = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, false, true))
        .unwrap();
    let unmodified = modifications
        .sequences
        .iter()
        .zip(&modifications.descriptions)
        .any(|(s, d)| s.to_string() == "AMEPTIDER" && d.is_empty());
    assert!(unmodified);

    entry.simple_variants.push(variant(5, 'L', "var1"));
    let combined = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, true, true))
        .unwrap();
    let variants_only = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, true, false))
        .unwrap();
    let modifications_only = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, false, true))
        .unwrap();
    assert!(combined.sequences.len() >= variants_only.sequences.len());
    assert!(combined.sequences.len() >= modifications_only.sequences.len());
}

/// Section `[PEFFEntry] digestWithVariants_no_variants_flag`.
#[test]
fn digest_without_the_variants_flag() {
    let mut entry = PEFFEntry {
        sequence: "MKAEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(4, 'L', "var1"));
    let digestor = trypsin();
    let with = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, true, false))
        .unwrap();
    let without = entry
        .digest_with_variants(&digestor, &digest_options(2, 20, true, false, false))
        .unwrap();
    assert!(without.sequences.len() < with.sequences.len());
    for description in &without.descriptions {
        assert!(description.is_empty());
    }
}

/// Section `[PEFFEntry] generatePeptides_basic`.
#[test]
fn generate_peptides_basic() {
    let mut entry = PEFFEntry {
        sequence: "MKAMEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(4, 'L', "dbSNP:test"));
    entry
        .modifications
        .push(PEFFModification::new(5, "UNIMOD:21", "Phospho"));
    let peptides = entry
        .generate_peptides(
            &trypsin(),
            &generation_options(&[], &[], 2, 2, 20, true, true, true),
        )
        .unwrap();
    assert!(!peptides.sequences.is_empty());
    assert!(
        peptides
            .descriptions
            .iter()
            .any(|d| d.is_empty() || !d.contains("M4L"))
    );
    assert!(peptides.descriptions.iter().any(|d| d.contains("M4L")));
}

/// Section `[PEFFEntry] generatePeptides_with_fixed_mods`.
#[test]
fn generate_peptides_with_fixed_modifications() {
    let entry = PEFFEntry {
        sequence: "MKACEPTIDER".into(),
        ..PEFFEntry::default()
    };
    let peptides = entry
        .generate_peptides(
            &trypsin(),
            &generation_options(&["Carbamidomethyl (C)"], &[], 0, 2, 20, true, false, false),
        )
        .unwrap();
    assert!(!peptides.sequences.is_empty());
    for sequence in &peptides.sequences {
        if sequence.as_str().contains('C') {
            assert!(sequence.is_modified());
        }
    }
}

/// Section `[PEFFEntry] generatePeptides_with_variable_mods`.
#[test]
fn generate_peptides_with_variable_modifications() {
    let entry = PEFFEntry {
        sequence: "MKAMEPTIDER".into(),
        ..PEFFEntry::default()
    };
    let digestor = trypsin();
    let with = entry
        .generate_peptides(
            &digestor,
            &generation_options(&[], &["Oxidation (M)"], 2, 2, 20, true, false, false),
        )
        .unwrap();
    assert!(!with.sequences.is_empty());
    let without = entry
        .generate_peptides(
            &digestor,
            &generation_options(&[], &[], 2, 2, 20, true, false, false),
        )
        .unwrap();
    assert!(with.sequences.len() >= without.sequences.len());
}

/// Section `[PEFFEntry] generatePeptides_full_workflow`.
#[test]
fn generate_peptides_full_workflow() {
    let mut entry = PEFFEntry {
        sequence: "MKACMEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry.simple_variants.push(variant(5, 'L', "var1"));
    let digestor = trypsin();
    let annotations_only = entry
        .generate_peptides(
            &digestor,
            &generation_options(&[], &[], 2, 2, 20, true, true, false),
        )
        .unwrap();
    let with_sample_mods = entry
        .generate_peptides(
            &digestor,
            &generation_options(
                &["Carbamidomethyl (C)"],
                &["Oxidation (M)"],
                1,
                2,
                20,
                true,
                true,
                false,
            ),
        )
        .unwrap();
    assert!(with_sample_mods.sequences.len() >= annotations_only.sequences.len());
}

/// Section `[PEFFEntry] generatePeptides_no_reference_with_peff_mods`.
#[test]
fn generate_peptides_without_reference_drops_the_unmodified_peptide() {
    let mut entry = PEFFEntry {
        sequence: "MKPEPTIDER".into(),
        ..PEFFEntry::default()
    };
    entry
        .modifications
        .push(PEFFModification::new(4, "UNIMOD:21", "Phospho"));
    let peptides = entry
        .generate_peptides(
            &trypsin(),
            &generation_options(&[], &[], 2, 2, 20, false, false, true),
        )
        .unwrap();
    let reference = AASequence::parse("PEPTIDER").unwrap().to_string();
    assert!(
        !peptides
            .sequences
            .iter()
            .any(|s| s.to_string() == reference)
    );
    assert!(!peptides.sequences.is_empty());
    for sequence in &peptides.sequences {
        assert!(sequence.is_modified());
    }
}

/// Sections `[EXTRA] round-trip: PEFFFile must consume the byte-exact output
/// of the UniPEFF TOPP tool` and `[EXTRA] AASequence materialization from
/// UniPEFF output -- what a top-down search engine sees`.
#[test]
fn unipeff_output_materialises() {
    let document = PEFFFile::new()
        .load(data("PEFFFile_UniPEFF_test.peff"))
        .unwrap();
    assert_eq!(document.entries.len(), 2);
    let e1 = &document.entries[0];

    let raw = e1.aa_sequence().unwrap();
    assert_eq!(raw.len(), 100);
    assert_eq!(raw.as_str(), e1.sequence);
    assert!(!raw.is_modified());

    let modified = e1.modified_sequence().unwrap();
    assert_eq!(modified.len(), 100);
    assert_eq!(modified.as_str(), e1.sequence);
    assert!(modified.is_modified());
    assert!(modified_residues(&modified) >= 1);
    let methyl = ModificationsDB::global()
        .get_modification("UNIMOD:34", Some('L'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let uplift = modified.mono_mass().unwrap() - raw.mono_mass().unwrap();
    assert!(uplift + 1e-6 >= methyl.diff_mono_mass());

    assert_eq!(
        e1.modifications.iter().filter(|m| m.position == 0).count(),
        3
    );

    let mature = e1.processed_sequence("PEFF:0001020").unwrap();
    assert_eq!(mature.len(), 60);
    assert_eq!(mature.as_str(), &e1.sequence[40..]);

    let e2 = &document.entries[1];
    assert!(e2.modifications.is_empty());
    let raw2 = e2.aa_sequence().unwrap();
    let modified2 = e2.modified_sequence().unwrap();
    assert_eq!(raw2.len(), 31);
    assert_eq!(modified2.len(), 31);
    assert!(!modified2.is_modified());
    assert!((raw2.mono_mass().unwrap() - modified2.mono_mass().unwrap()).abs() < 1e-5);
}

/// Section `[EXTRA] AASequence materialization from UniPEFF output`, the
/// controlled PSI-MOD regression: `MOD:00046` on S8 aliases Unimod Phospho.
#[test]
fn psi_mod_accession_resolves_through_its_unimod_alias() {
    let db = with_psi_mod_00046();
    let mut crafted = PEFFEntry {
        sequence: "PEPTIDES".into(),
        ..PEFFEntry::default()
    };
    crafted
        .modifications
        .push(PEFFModification::new(8, "MOD:00046", "O-phospho-L-serine"));
    let reference = crafted.aa_sequence().unwrap();
    let (phospho, warnings) = crafted.modified_sequence_with_registry(&db).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(phospho.len(), 8);
    assert!(phospho.residue_modification(7).unwrap().is_some());
    let phospho_on_s = ModificationsDB::global()
        .get_modification("UNIMOD:21", Some('S'), Some(TermSpecificity::Anywhere))
        .unwrap();
    let delta = phospho.mono_mass().unwrap() - reference.mono_mass().unwrap();
    assert!((delta - phospho_on_s.diff_mono_mass()).abs() < 1e-5);
}

/// Without PSI-MOD, the same accession does not resolve: the source's
/// warning, and the sequence unmodified.
#[test]
fn unresolved_modifications_warn_and_are_skipped() {
    let mut crafted = PEFFEntry {
        sequence: "PEPTIDES".into(),
        ..PEFFEntry::default()
    };
    crafted
        .modifications
        .push(PEFFModification::new(8, "MOD:00046", "O-phospho-L-serine"));
    let (sequence, warnings) = crafted
        .modified_sequence_with_registry(ModificationsDB::global())
        .unwrap();
    assert!(!sequence.is_modified());
    assert_eq!(
        warnings,
        ["Could not resolve modification: MOD:00046 (O-phospho-L-serine)"]
    );
}

const REF_TABLE: &[(u32, &str, &str)] = &[
    (
        0x00,
        "MKLSLGLLPKARGSPKERLAPTDVMGNECS",
        "MKLSLGLLPKARGSPKERLAPTDVMGNECS",
    ),
    (
        0x01,
        "MKLSLGLLPKARGS[UNIMOD:21]PKERLAPTDVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PKERLAPTDVMGNECS",
    ),
    (
        0x02,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ERLAPTDVMGNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ERLAPTDVMGNECS",
    ),
    (
        0x03,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ERLAPTDVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ERLAPTDVMGNECS",
    ),
    (
        0x04,
        "MKLSLGLLPKARGSPKER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGSPKER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x05,
        "MKLSLGLLPKARGS[UNIMOD:21]PKER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PKER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x06,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x07,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x08,
        "MKLSLGLLPKARGSPKERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGSPKERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x09,
        "MKLSLGLLPKARGS[UNIMOD:21]PKERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PKERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x0A,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x0B,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x0C,
        "MKLSLGLLPKARGSPKER[UNIMOD:34]LAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGSPKER(Methyl)LAPT(Phospho)DVMGNECS",
    ),
    (
        0x0D,
        "MKLSLGLLPKARGS[UNIMOD:21]PKER[UNIMOD:34]LAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PKER(Methyl)LAPT(Phospho)DVMGNECS",
    ),
    (
        0x0E,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ER(Methyl)LAPT(Phospho)DVMGNECS",
    ),
    (
        0x0F,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ER(Methyl)LAPT(Phospho)DVMGNECS",
    ),
    (
        0x10,
        "MKLSLGLLPKARGSPKERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPKERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x11,
        "MKLSLGLLPKARGS[UNIMOD:21]PKERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PKERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x12,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x13,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x14,
        "MKLSLGLLPKARGSPKER[UNIMOD:34]LAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPKER(Methyl)LAPTDVM(Oxidation)GNECS",
    ),
    (
        0x15,
        "MKLSLGLLPKARGS[UNIMOD:21]PKER[UNIMOD:34]LAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PKER(Methyl)LAPTDVM(Oxidation)GNECS",
    ),
    (
        0x16,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ER[UNIMOD:34]LAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ER(Methyl)LAPTDVM(Oxidation)GNECS",
    ),
    (
        0x17,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ER[UNIMOD:34]LAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ER(Methyl)LAPTDVM(Oxidation)GNECS",
    ),
    (
        0x18,
        "MKLSLGLLPKARGSPKERLAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPKERLAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x19,
        "MKLSLGLLPKARGS[UNIMOD:21]PKERLAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PKERLAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1A,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ERLAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ERLAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1B,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ERLAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ERLAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1C,
        "MKLSLGLLPKARGSPKER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPKER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1D,
        "MKLSLGLLPKARGS[UNIMOD:21]PKER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PKER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1E,
        "MKLSLGLLPKARGSPK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPK(Acetyl)ER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0x1F,
        "MKLSLGLLPKARGS[UNIMOD:21]PK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PK(Acetyl)ER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
];

const K16R_TABLE: &[(u32, &str, &str)] = &[
    (
        0x0,
        "MKLSLGLLPKARGSPRERLAPTDVMGNECS",
        "MKLSLGLLPKARGSPRERLAPTDVMGNECS",
    ),
    (
        0x1,
        "MKLSLGLLPKARGS[UNIMOD:21]PRERLAPTDVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PRERLAPTDVMGNECS",
    ),
    (
        0x2,
        "MKLSLGLLPKARGSPRER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGSPRER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x4,
        "MKLSLGLLPKARGSPRERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGSPRERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x8,
        "MKLSLGLLPKARGSPRERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPRERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x5,
        "MKLSLGLLPKARGS[UNIMOD:21]PRERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGS(Phospho)PRERLAPT(Phospho)DVMGNECS",
    ),
    (
        0xA,
        "MKLSLGLLPKARGSPRER[UNIMOD:34]LAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGSPRER(Methyl)LAPTDVM(Oxidation)GNECS",
    ),
    (
        0xF,
        "MKLSLGLLPKARGS[UNIMOD:21]PRER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGS(Phospho)PRER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
];

const S14A_TABLE: &[(u32, &str, &str)] = &[
    (
        0x0,
        "MKLSLGLLPKARGAPKERLAPTDVMGNECS",
        "MKLSLGLLPKARGAPKERLAPTDVMGNECS",
    ),
    (
        0x1,
        "MKLSLGLLPKARGAPK[UNIMOD:1]ERLAPTDVMGNECS",
        "MKLSLGLLPKARGAPK(Acetyl)ERLAPTDVMGNECS",
    ),
    (
        0x2,
        "MKLSLGLLPKARGAPKER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGAPKER(Methyl)LAPTDVMGNECS",
    ),
    (
        0x4,
        "MKLSLGLLPKARGAPKERLAPT[UNIMOD:21]DVMGNECS",
        "MKLSLGLLPKARGAPKERLAPT(Phospho)DVMGNECS",
    ),
    (
        0x8,
        "MKLSLGLLPKARGAPKERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGAPKERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x3,
        "MKLSLGLLPKARGAPK[UNIMOD:1]ER[UNIMOD:34]LAPTDVMGNECS",
        "MKLSLGLLPKARGAPK(Acetyl)ER(Methyl)LAPTDVMGNECS",
    ),
    (
        0xC,
        "MKLSLGLLPKARGAPKERLAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGAPKERLAPT(Phospho)DVM(Oxidation)GNECS",
    ),
    (
        0xF,
        "MKLSLGLLPKARGAPK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSLGLLPKARGAPK(Acetyl)ER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
];

const CV_TABLE: &[(u32, &str, &str)] = &[
    (
        0x00,
        "MKLSABLPKARGSPKERLAPTDVMGNECS",
        "MKLSABLPKARGSPKERLAPTDVMGNECS",
    ),
    (
        0x01,
        "MKLSABLPKARGS[UNIMOD:21]PKERLAPTDVMGNECS",
        "MKLSABLPKARGS(Phospho)PKERLAPTDVMGNECS",
    ),
    (
        0x10,
        "MKLSABLPKARGSPKERLAPTDVM[UNIMOD:35]GNECS",
        "MKLSABLPKARGSPKERLAPTDVM(Oxidation)GNECS",
    ),
    (
        0x1F,
        "MKLSABLPKARGS[UNIMOD:21]PK[UNIMOD:1]ER[UNIMOD:34]LAPT[UNIMOD:21]DVM[UNIMOD:35]GNECS",
        "MKLSABLPKARGS(Phospho)PK(Acetyl)ER(Methyl)LAPT(Phospho)DVM(Oxidation)GNECS",
    ),
];

const REF_SEQ: &str = "MKLSLGLLPKARGSPKERLAPTDVMGNECS";

type ModDef = (usize, &'static str, &'static str);

const REF_MODS: &[ModDef] = &[
    (14, "UNIMOD:21", "Phospho"),
    (16, "UNIMOD:1", "Acetyl"),
    (18, "UNIMOD:34", "Methyl"),
    (22, "UNIMOD:21", "Phospho"),
    (25, "UNIMOD:35", "Oxidation"),
];
const K16R_MODS: &[ModDef] = &[
    (14, "UNIMOD:21", "Phospho"),
    (18, "UNIMOD:34", "Methyl"),
    (22, "UNIMOD:21", "Phospho"),
    (25, "UNIMOD:35", "Oxidation"),
];
const S14A_MODS: &[ModDef] = &[
    (16, "UNIMOD:1", "Acetyl"),
    (18, "UNIMOD:34", "Methyl"),
    (22, "UNIMOD:21", "Phospho"),
    (25, "UNIMOD:35", "Oxidation"),
];
const CV_MODS: &[ModDef] = &[
    (13, "UNIMOD:21", "Phospho"),
    (15, "UNIMOD:1", "Acetyl"),
    (17, "UNIMOD:34", "Methyl"),
    (21, "UNIMOD:21", "Phospho"),
    (24, "UNIMOD:35", "Oxidation"),
];

fn build(sequence: &str, definitions: &[ModDef], mask: u32) -> PEFFEntry {
    let mut entry = PEFFEntry {
        identifier: "sp:Q9TEST1".into(),
        prefix: "sp".into(),
        sequence: sequence.into(),
        sequence_length: sequence.len(),
        ..PEFFEntry::default()
    };
    for (bit, (position, accession, name)) in definitions.iter().enumerate() {
        if mask & (1 << bit) != 0 {
            entry
                .modifications
                .push(PEFFModification::new(*position, *accession, *name));
        }
    }
    entry
}

fn check_table(sequence: &str, definitions: &[ModDef], table: &[(u32, &str, &str)]) {
    for &(mask, proforma, aaseq) in table {
        let entry = build(sequence, definitions, mask);
        assert_eq!(PEFFFile::to_proforma(&entry), proforma, "mask {mask:#x}");
        let materialised = entry.modified_sequence().unwrap();
        assert_eq!(materialised.to_string(), aaseq, "mask {mask:#x}");
        assert_eq!(materialised.as_str(), sequence);
    }
}

/// Section `[EXTRA] proteoform enumeration -- exact ProForma strings for a
/// realistic UniProt-style entry (50+ proteoforms)`, blocks A to F.
#[test]
fn proteoform_tables() {
    assert_eq!(REF_SEQ.len(), 30);
    check_table(REF_SEQ, REF_MODS, REF_TABLE);

    // Block B-pre: the variant backbones come from the variant API.
    let mut variants_only = PEFFEntry {
        identifier: "sp:Q9TEST1".into(),
        prefix: "sp".into(),
        sequence: REF_SEQ.into(),
        sequence_length: REF_SEQ.len(),
        ..PEFFEntry::default()
    };
    variants_only.simple_variants = vec![
        variant(16, 'R', "dbSNP:rs_test_K16R"),
        variant(14, 'A', "dbSNP:rs_test_S14A"),
    ];
    variants_only.complex_variants = vec![PEFFVariantComplex::new(
        5,
        7,
        "AB",
        "ClinVar:test_SLG_AB",
        openms::format::peff::ANNOTATION_ID_NOT_SET,
    )];
    let simple = variants_only.variant_sequences(false);
    assert_eq!(simple.sequences.len(), 2);
    assert_eq!(
        simple.sequences[0].as_str(),
        "MKLSLGLLPKARGSPRERLAPTDVMGNECS"
    );
    assert_eq!(
        simple.sequences[1].as_str(),
        "MKLSLGLLPKARGAPKERLAPTDVMGNECS"
    );
    assert!(simple.descriptions[0].contains("K16R"));
    assert!(simple.descriptions[0].contains("dbSNP:rs_test_K16R"));
    assert!(simple.descriptions[1].contains("S14A"));
    assert!(simple.descriptions[1].contains("dbSNP:rs_test_S14A"));
    let all = variants_only.variant_sequences(true);
    assert_eq!(all.sequences.len(), 3);
    assert_eq!(all.sequences[2].as_str(), "MKLSABLPKARGSPKERLAPTDVMGNECS");
    assert!(all.descriptions[2].contains("5-7>AB"));
    assert!(all.descriptions[2].contains("ClinVar:test_SLG_AB"));

    // Blocks B, C and D.
    let k16r = all.sequences[0].as_str().to_owned();
    assert_eq!(k16r.len(), 30);
    assert_eq!(&k16r[15..16], "R");
    check_table(&k16r, K16R_MODS, K16R_TABLE);
    let s14a = all.sequences[1].as_str().to_owned();
    assert_eq!(&s14a[13..14], "A");
    check_table(&s14a, S14A_MODS, S14A_TABLE);
    let cv = all.sequences[2].as_str().to_owned();
    assert_eq!(cv, "MKLSABLPKARGSPKERLAPTDVMGNECS");
    assert_eq!(cv.len(), 29);
    check_table(&cv, CV_MODS, CV_TABLE);

    // Block E: processed regions and a PSI-MOD + Unimod mix.
    let mut full = PEFFEntry {
        identifier: "sp:Q9TEST1".into(),
        prefix: "sp".into(),
        sequence: REF_SEQ.into(),
        sequence_length: REF_SEQ.len(),
        ..PEFFEntry::default()
    };
    full.modifications = vec![
        PEFFModification::new(14, "MOD:00046", "O-phospho-L-serine"),
        PEFFModification::new(25, "UNIMOD:35", "Oxidation"),
    ];
    let region = |start, end, accession: &str, name: &str| {
        PEFFProcessedRegion::new(
            start,
            end,
            accession,
            name,
            "",
            openms::format::peff::ANNOTATION_ID_NOT_SET,
        )
    };
    full.processed_regions = vec![
        region(1, 1, "PEFF:0001035", "initiator methionine"),
        region(1, 9, "PEFF:0001021", "signal peptide"),
        region(10, 30, "PEFF:0001020", "mature protein"),
        region(11, 15, "PEFF:0001022", "transit peptide"),
        region(21, 25, "PEFF:0001034", "propeptide"),
    ];
    let (full_aa, _) = full
        .modified_sequence_with_registry(&with_psi_mod_00046())
        .unwrap();
    assert_eq!(
        full_aa.to_string(),
        "MKLSLGLLPKARGS(Phospho)PKERLAPTDVM(Oxidation)GNECS"
    );
    assert_eq!(full_aa.as_str(), REF_SEQ);
    assert_eq!(
        PEFFFile::to_proforma(&full),
        "MKLSLGLLPKARGS[MOD:00046]PKERLAPTDVM[UNIMOD:35]GNECS"
    );
    let chain = full.processed_sequence("PEFF:0001020").unwrap();
    assert_eq!(chain.len(), 21);
    assert_eq!(chain.to_string(), "KARGSPKERLAPTDVMGNECS");
    assert!(!chain.is_modified());
    assert_eq!(&chain.as_str()[4..5], "S");
    assert!(chain.residue_modification(4).unwrap().is_none());
    assert_eq!(
        full.processed_sequence("PEFF:0001035").unwrap().to_string(),
        "M"
    );
    assert_eq!(
        full.processed_sequence("PEFF:0001022").unwrap().to_string(),
        "ARGSP"
    );
    assert_eq!(
        full.processed_sequence("PEFF:0001034").unwrap().to_string(),
        "PTDVM"
    );

    // Block F: unlocalised modifications and disulfide bonds.
    let mut unlocalised = PEFFEntry {
        identifier: "sp:Q9TEST1-UNLOC".into(),
        prefix: "sp".into(),
        sequence: "PEPTIDES".into(),
        sequence_length: 8,
        ..PEFFEntry::default()
    };
    unlocalised.modifications = vec![
        PEFFModification::new(0, "UNIMOD:21", "Phospho"),
        PEFFModification::new(0, "UNIMOD:35", "Oxidation"),
    ];
    assert_eq!(
        PEFFFile::to_proforma(&unlocalised),
        "<?[UNIMOD:21][UNIMOD:35]>PEPTIDES"
    );
    let unlocalised_aa = unlocalised.modified_sequence().unwrap();
    assert_eq!(unlocalised_aa.to_string(), "PEPTIDES");
    assert!(!unlocalised_aa.is_modified());

    let mut mixed = unlocalised.clone();
    mixed.identifier = "sp:Q9TEST1-MIXED".into();
    mixed.modifications = vec![
        PEFFModification::new(0, "UNIMOD:35", "Oxidation"),
        PEFFModification::new(8, "UNIMOD:21", "Phospho"),
    ];
    assert_eq!(
        PEFFFile::to_proforma(&mixed),
        "<?[UNIMOD:35]>PEPTIDES[UNIMOD:21]"
    );
    assert_eq!(
        mixed.modified_sequence().unwrap().to_string(),
        "PEPTIDES(Phospho)"
    );

    let mut with_bonds = PEFFEntry {
        identifier: "sp:Q9TEST1-DSB".into(),
        prefix: "sp".into(),
        sequence: "PEPCIDESCPEP".into(),
        sequence_length: 12,
        ..PEFFEntry::default()
    };
    with_bonds.disulfide_bonds = vec![
        openms::format::peff::PEFFDisulfideBond::new(
            "4",
            "9",
            "intra",
            openms::format::peff::ANNOTATION_ID_NOT_SET,
        ),
        openms::format::peff::PEFFDisulfideBond::new(
            "?",
            "9",
            "endpoint-unknown",
            openms::format::peff::ANNOTATION_ID_NOT_SET,
        ),
    ];
    assert_eq!(PEFFFile::to_proforma(&with_bonds), "PEPCIDESCPEP");
    let bonded = with_bonds.modified_sequence().unwrap();
    assert_eq!(bonded.to_string(), "PEPCIDESCPEP");
    assert!(!bonded.is_modified());
}

/// Section `[EXTRA] proteoform enumeration -- per-proteoform PTM cap
/// (generatePeptides max_peff_mods_per_proteoform)`.
#[test]
fn proteoform_cap() {
    const SEQ: &str = "ASASASASASASASASAS";
    let mut entry = PEFFEntry {
        identifier: "sp:P02662_TEST".into(),
        prefix: "sp".into(),
        sequence: SEQ.into(),
        sequence_length: SEQ.len(),
        ..PEFFEntry::default()
    };
    for position in [2, 4, 6, 8, 10, 12, 14, 16, 18] {
        assert_eq!(&SEQ[position - 1..position], "S");
        entry
            .modifications
            .push(PEFFModification::new(position, "UNIMOD:21", "Phospho"));
    }
    let mut no_cleavage = ProteaseDigestion::default();
    no_cleavage.set_enzyme("no cleavage").unwrap();
    let generate = |cap: usize| {
        let options = PeptideGenerationOptions {
            max_variable_mods_per_peptide: 0,
            min_length: 1,
            max_length: 0,
            include_reference: true,
            include_peff_variants: false,
            include_peff_modifications: true,
            max_peff_mods_per_proteoform: cap,
            ..PeptideGenerationOptions::default()
        };
        entry.generate_peptides(&no_cleavage, &options).unwrap()
    };
    let count = |cap: usize| {
        let peptides = generate(cap);
        assert_eq!(peptides.descriptions.len(), peptides.sequences.len());
        if cap > 0 {
            for sequence in &peptides.sequences {
                assert_eq!(sequence.as_str(), SEQ);
                assert!(modified_residues(sequence) <= cap);
            }
        }
        peptides.sequences.len()
    };
    assert_eq!(count(1), 10);
    assert_eq!(count(2), 46);
    assert_eq!(count(3), 130);
    assert_eq!(count(0), 512);

    let capped = generate(3);
    assert_eq!(capped.sequences.len(), 130);
    assert!(
        capped
            .sequences
            .iter()
            .any(|s| !s.is_modified() && s.to_string() == SEQ)
    );
    assert!(
        capped
            .sequences
            .iter()
            .any(|s| s.to_string() == "AS(Phospho)AS(Phospho)AS(Phospho)ASASASASASAS")
    );
}

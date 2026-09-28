// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FragmentIndex` against the upstream class test
//! `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` at core
//! `bc9cc12` (evidence tier 3: its expected peptide sets, counts, ranges and
//! self-hit invariants, transcribed section by section with source line
//! numbers), plus native checks (tier 4).

use openms::analysis::id::fragment_index::{
    FragmentIndex, FragmentIndexParameters, Peptide, SnesAnchor, SpectrumMatchesTopN,
};
use openms::chemistry::{
    AASequence, C13C12_MASSDIFF_U, DigestionSpecificity, ModifiedPeptideGenerator, PROTON_MASS_U,
    TheoreticalSpectrumGenerator,
};
use openms::format::fasta::FASTAEntry;
use openms::{MSSpectrum, Precursor, Result};

fn entry(id: &str, sequence: &str) -> FASTAEntry {
    FASTAEntry::new(id, id, sequence)
}

fn index(params: FragmentIndexParameters) -> FragmentIndex {
    FragmentIndex::new(params).expect("valid parameters")
}

/// `FragmentIndex_test::testDigestion`: same set by `{sequence_, mod_bitmask_}`.
fn test_digestion(fi: &FragmentIndex, expected: &[(u32, u32, (u16, u16))]) -> bool {
    expected.len() == fi.peptides().len()
        && expected.iter().all(|&(_, mask, seq)| {
            fi.peptides()
                .iter()
                .any(|p| p.sequence == seq && p.mod_bitmask == mask)
        })
}

/// `FragmentIndex_test::peptidesSorted`.
fn peptides_sorted(fi: &FragmentIndex) -> bool {
    fi.peptides()
        .windows(2)
        .all(|w| w[1].precursor_mz >= w[0].precursor_mz)
}

/// `FragmentIndex_test::fragmentsSorted`.
fn fragments_sorted(fi: &FragmentIndex) -> bool {
    fi.fragments()
        .chunks(fi.bucket_size())
        .all(|b| b.windows(2).all(|w| w[1].peptide_idx >= w[0].peptide_idx))
}

/// `FragmentIndex_test::fragmentCountForPeptide`.
fn fragment_count(fi: &FragmentIndex, peptide_idx: u32) -> usize {
    fi.fragments()
        .iter()
        .filter(|f| f.peptide_idx == peptide_idx)
        .count()
}

fn spectrum(
    peptide: &AASequence,
    charge: u8,
    precursor_mz: f64,
    known: bool,
) -> Result<MSSpectrum> {
    let ions = TheoreticalSpectrumGenerator::default().generate(peptide, charge, charge, None)?;
    Ok(MSSpectrum {
        ms_level: 2,
        peaks: ions.peaks,
        precursors: vec![Precursor {
            mz: precursor_mz,
            charge: if known { i32::from(charge) } else { 0 },
            ..Precursor::default()
        }],
        ..MSSpectrum::default()
    })
}

/// `FragmentIndex_test::testQuery`: every peptide's own theoretical spectrum
/// at `charge` finds it at that charge with every peak matched.
fn test_query(fi: &FragmentIndex, charge: u8, known: bool, entries: &[FASTAEntry]) -> Result<bool> {
    let mut all = true;
    for (idx, pep) in fi.peptides().iter().enumerate() {
        let modified = fi.reconstruct_modified_sequence(pep, entries)?;
        let spec = spectrum(&modified, charge, modified.mz(i32::from(charge))?, known)?;
        let mut sms = SpectrumMatchesTopN::default();
        fi.query_spectrum(&spec, &mut sms)?;
        let mut found = false;
        for s in &sms.hits {
            if s.peptide_idx == idx && s.precursor_charge == u16::from(charge) {
                found = s.num_matched as usize >= spec.peaks.len();
            }
        }
        all = all && found;
    }
    Ok(all)
}

fn no_mods(mut p: FragmentIndexParameters) -> FragmentIndexParameters {
    p.modifications_fixed.clear();
    p.modifications_variable.clear();
    p
}

// FragmentIndex_test.cpp:178-262 (build())
#[test]
fn build_digests_as_class_test() -> Result<()> {
    let entries = [
        entry("t", "ARGEPADSSRKDFDMDMDM"),
        entry("t2", "HALLORTSCHSM"),
    ];
    let with_mods: [(u32, u32, (u16, u16)); 11] = [
        (0, 0, (2, 8)),
        (0, 0, (11, 8)),
        (0, 1, (11, 8)),
        (0, 2, (11, 8)),
        (0, 3, (11, 8)),
        (0, 4, (11, 8)),
        (0, 5, (11, 8)),
        (0, 6, (11, 8)),
        (1, 0, (0, 6)),
        (1, 0, (6, 6)),
        (1, 1, (6, 6)),
    ];
    let unmod_no_minmax = [
        (0, 0, (0, 2)),
        (0, 0, (2, 8)),
        (0, 0, (10, 1)),
        (0, 0, (11, 8)),
        (1, 0, (0, 6)),
        (1, 0, (6, 6)),
    ];
    let unmod_minmax = [(0, 0, (0, 2)), (1, 0, (0, 6)), (1, 0, (6, 6))];
    let unmod_missed = [
        (0, 0, (0, 2)),
        (0, 0, (2, 8)),
        (0, 0, (11, 8)),
        (0, 0, (0, 10)),
        (0, 0, (2, 9)),
        (0, 0, (10, 9)),
        (1, 0, (0, 6)),
        (1, 0, (6, 6)),
        (1, 0, (0, 12)),
    ];
    let mut p = no_mods(FragmentIndexParameters::default());
    p.enzyme = "Trypsin".into();
    p.missed_cleavages = 0;
    p.peptide_min_mass = 0;
    p.peptide_min_length = 0;
    p.peptide_max_mass = 5000;
    let mut fi = index(p.clone());
    fi.build(&entries)?;
    assert!(test_digestion(&fi, &unmod_no_minmax));
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));

    fi.clear();
    p.peptide_min_length = 2;
    p.peptide_max_length = 6;
    fi.set_parameters(p.clone())?;
    fi.build(&entries)?;
    assert!(test_digestion(&fi, &unmod_minmax));
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));

    fi.clear();
    p.peptide_max_length = 100;
    p.missed_cleavages = 1;
    fi.set_parameters(p.clone())?;
    fi.build(&entries)?;
    assert!(test_digestion(&fi, &unmod_missed));
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));

    fi.clear();
    p.missed_cleavages = 0;
    p.peptide_min_mass = 0;
    p.peptide_min_length = 6;
    p.modifications_variable = vec!["Oxidation (M)".into()];
    p.modifications_fixed = vec!["Carbamidomethyl (C)".into()];
    fi.set_parameters(p)?;
    fi.build(&entries)?;
    assert!(test_digestion(&fi, &with_mods));
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));
    Ok(())
}

fn specificity_params(
    specificity: DigestionSpecificity,
    min: usize,
    max: usize,
) -> FragmentIndexParameters {
    let mut p = no_mods(FragmentIndexParameters::default());
    p.enzyme = "Trypsin".into();
    p.missed_cleavages = 0;
    p.enzyme_specificity = specificity;
    p.peptide_min_length = min;
    p.peptide_max_length = max;
    p.peptide_min_mass = 0;
    p.peptide_max_mass = 50000;
    p
}

// FragmentIndex_test.cpp:270-375 ([EXTRA] peptide:enzyme_specificity)
#[test]
fn enzyme_specificity_changes_digestion() -> Result<()> {
    let entries = [entry("t", "AKACDEFGRHILMNPQSTV")];
    let count = |p: FragmentIndexParameters, e: &[FASTAEntry]| -> Result<usize> {
        let mut fi = index(p);
        fi.build(e)?;
        Ok(fi.peptides().len())
    };
    assert_eq!(
        count(
            specificity_params(DigestionSpecificity::Full, 2, 100),
            &entries
        )?,
        3
    );
    assert!(
        count(
            specificity_params(DigestionSpecificity::Semi, 2, 100),
            &entries
        )? > 3
    );
    assert_eq!(
        count(
            specificity_params(DigestionSpecificity::None, 8, 12),
            &entries
        )?,
        50
    );
    assert!(
        50 > count(
            specificity_params(DigestionSpecificity::Semi, 8, 12),
            &entries
        )?
    );
    assert_eq!(
        count(
            specificity_params(DigestionSpecificity::None, 8, 12),
            &[entry("t", "ABC")]
        )?,
        0
    );
    Ok(())
}

// FragmentIndex_test.cpp:379-392
#[test]
fn build_rejects_entries_longer_than_65535() -> Result<()> {
    let mut fi = index(FragmentIndexParameters::default());
    assert!(fi.build(&[entry("contig1", &"A".repeat(70000))]).is_err());
    let mut ok = index(FragmentIndexParameters::default());
    ok.build(&[entry("contig_ok", &"A".repeat(65535))])?;
    assert!(ok.is_build());
    Ok(())
}

// FragmentIndex_test.cpp:395-406
#[test]
fn clear_resets_the_index() -> Result<()> {
    let mut fi = index(FragmentIndexParameters::default());
    fi.build(&[
        entry("t", "ARGEPADSSRKDFDMDMDM"),
        entry("t2", "HALLORTSCHS"),
    ])?;
    fi.clear();
    assert!(fi.peptides().is_empty());
    assert!(!fi.is_build());
    Ok(())
}

const PROTEIN: &str = "MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYDKVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVLFDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDYHTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMALARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDSVWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN";

// FragmentIndex_test.cpp:410-439
#[test]
fn query_spectrum_self_hits_every_charge() -> Result<()> {
    let entries = [entry("test1", PROTEIN)];
    let mut p = FragmentIndexParameters::default();
    p.max_fragment_charge = 4;
    p.min_precursor_charge = 1;
    p.max_precursor_charge = 4;
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 5_000_000;
    p.min_ion_index = 0;
    let mut fi = index(p);
    fi.build(&entries)?;
    assert!(!fi.peptides().is_empty());
    for charge in 1..=4 {
        assert!(
            test_query(&fi, charge, false, &entries)?,
            "charge {charge} unknown"
        );
        assert!(
            test_query(&fi, charge, true, &entries)?,
            "charge {charge} known"
        );
    }
    Ok(())
}

// FragmentIndex_test.cpp:443-499
#[test]
fn isotope_errors_map_to_the_same_peptide() -> Result<()> {
    let entries = [entry("test1", PROTEIN)];
    let mut p = no_mods(FragmentIndexParameters::default());
    p.min_isotope_error = -3;
    p.max_isotope_error = 3;
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 90000;
    let mut fi = index(p);
    fi.build(&entries)?;
    let peptide = AASequence::parse("EVAEAATGEDASSPPPK")?;
    for iso in -3..=3 {
        let spec = spectrum(
            &peptide,
            1,
            peptide.mz(1)? + f64::from(iso) * C13C12_MASSDIFF_U,
            true,
        )?;
        let mut sms = SpectrumMatchesTopN::default();
        fi.query_spectrum(&spec, &mut sms)?;
        assert!(!sms.hits.is_empty(), "isotope {iso}");
        for hit in &sms.hits {
            let result = fi.peptides()[hit.peptide_idx];
            assert_eq!(result.sequence.0, 5);
            assert_eq!(usize::from(result.sequence.1), peptide.len());
        }
    }
    Ok(())
}

// FragmentIndex_test.cpp:503-569
#[test]
fn tolerances_absorb_jitter_and_precursor_offset() -> Result<()> {
    let entries = [entry("test1", PROTEIN)];
    let mut p = no_mods(FragmentIndexParameters::default());
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 90000;
    p.min_ion_index = 0;
    p.fragment_mz_tolerance = 0.05;
    p.fragment_mz_tolerance_unit_ppm = false;
    p.precursor_mass_tolerance_lower = 2.0;
    p.precursor_mass_tolerance_upper = 2.0;
    p.precursor_mass_tolerance_unit_ppm = false;
    let mut fi = index(p);
    fi.build(&entries)?;
    let peptide = AASequence::parse("EVAEAATGEDASSPPPK")?;
    let mut spec = spectrum(&peptide, 1, peptide.mz(1)? + 1.9, true)?;
    for (i, peak) in spec.peaks.iter_mut().enumerate() {
        let factor = ((i % 91) as i32 - 45) as f32 * 0.001_f32;
        peak.mz += f64::from(factor);
    }
    let mut sms = SpectrumMatchesTopN::default();
    fi.query_spectrum(&spec, &mut sms)?;
    let mut found = false;
    for hit in &sms.hits {
        let sequence = fi.peptides()[hit.peptide_idx].sequence;
        if sequence.0 == 5 && usize::from(sequence.1) == peptide.len() && hit.isotope_error == 0 {
            found = true;
            assert!(hit.num_matched as usize >= spec.peaks.len());
        }
    }
    assert!(found);
    Ok(())
}

fn no_cleavage() -> FragmentIndexParameters {
    let mut p = no_mods(FragmentIndexParameters::default());
    p.enzyme = "no cleavage".into();
    p.peptide_min_length = 0;
    p.peptide_max_length = 100;
    p.peptide_min_mass = 0;
    p.peptide_max_mass = 50000;
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 50000;
    p.min_ion_index = 0;
    p
}

// FragmentIndex_test.cpp:573-611
#[test]
fn lightweight_fragment_count() -> Result<()> {
    let entries = [entry("p", "PEPTIDER")];
    let mut fi = index(no_cleavage());
    fi.build(&entries)?;
    assert_eq!(fi.peptides().len(), 1);
    assert_eq!(fragment_count(&fi, 0), 2 * 7);
    fi.clear();
    let mut p = no_cleavage();
    p.min_ion_index = 2;
    fi.set_parameters(p)?;
    fi.build(&entries)?;
    assert_eq!(fi.peptides().len(), 1);
    assert_eq!(fragment_count(&fi, 0), 2 * 5);
    Ok(())
}

// FragmentIndex_test.cpp:615-645
#[test]
fn multi_mod_per_site() -> Result<()> {
    let mut p = no_cleavage();
    p.modifications_max_variable_mods_per_peptide = 2;
    p.modifications_variable = vec!["Oxidation (M)".into()];
    p.modifications_fixed = vec!["Carbamidomethyl (C)".into()];
    let mut fi = index(p);
    fi.build(&[entry("p", "ACACK")])?;
    assert_eq!(fi.peptides().len(), 1);
    assert_eq!(fi.peptides()[0].mod_bitmask, 0);
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));
    assert_eq!(fragment_count(&fi, 0), 8);
    Ok(())
}

// FragmentIndex_test.cpp:649-695
#[test]
fn fixed_plus_variable_mods() -> Result<()> {
    let entries = [entry("p", "ACMK")];
    let mut p = no_cleavage();
    p.modifications_max_variable_mods_per_peptide = 2;
    p.modifications_variable = vec!["Oxidation (M)".into()];
    p.modifications_fixed = vec!["Carbamidomethyl (C)".into()];
    let mut fi = index(p);
    fi.build(&entries)?;
    assert_eq!(fi.peptides().len(), 2);
    assert!(peptides_sorted(&fi));
    assert!(fragments_sorted(&fi));
    for i in 0..fi.peptides().len() {
        assert_eq!(fragment_count(&fi, i as u32), 6);
    }
    for pep in fi.peptides() {
        let reconstructed = fi.reconstruct_modified_sequence(pep, &entries)?;
        assert_eq!(reconstructed.len(), 4);
        assert!(reconstructed.residue_modification(1)?.is_some());
        assert_eq!(
            reconstructed.residue_modification(2)?.is_some(),
            pep.mod_bitmask & 1 != 0
        );
    }
    Ok(())
}

// FragmentIndex_test.cpp:700-845 (cross_validate_vs_ModifiedPeptideGenerator)
#[test]
fn cross_validates_against_modified_peptide_generator() -> Result<()> {
    let run_mpg =
        |seq: &str, fixed: &[&str], variable: &[&str], max: usize| -> Result<Vec<(f32, String)>> {
            let generator = ModifiedPeptideGenerator::default();
            let mut peptide = AASequence::parse(seq)?;
            if !fixed.is_empty() {
                generator.apply_fixed_modifications(
                    &ModifiedPeptideGenerator::get_modifications(fixed)?,
                    &mut peptide,
                )?;
            }
            let variants = if variable.is_empty() {
                vec![peptide]
            } else {
                generator.variable_modifications(
                    &ModifiedPeptideGenerator::get_modifications(variable)?,
                    &peptide,
                    max,
                    true,
                )?
            };
            let mut out = variants
                .iter()
                .map(|v| Ok((v.mz(1)? as f32, v.to_string())))
                .collect::<Result<Vec<_>>>()?;
            out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            Ok(out)
        };
    let run_fi =
        |seq: &str, fixed: &[&str], variable: &[&str], max: usize| -> Result<Vec<(f32, String)>> {
            let entries = [entry("p", seq)];
            let mut p = no_cleavage();
            p.min_ion_index = 2;
            p.modifications_max_variable_mods_per_peptide = max;
            p.modifications_variable = variable.iter().map(|s| s.to_string()).collect();
            p.modifications_fixed = fixed.iter().map(|s| s.to_string()).collect();
            let mut fi = index(p);
            fi.build(&entries)?;
            let mut out = fi
                .peptides()
                .iter()
                .map(|pep| {
                    Ok((
                        pep.precursor_mz,
                        fi.reconstruct_modified_sequence(pep, &entries)?.to_string(),
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            Ok(out)
        };
    let cases: [(&str, &[&str], &[&str], usize); 6] = [
        ("ACMACK", &["Carbamidomethyl (C)"], &["Oxidation (M)"], 2),
        ("DFDMDMDM", &[], &["Oxidation (M)"], 2),
        (
            "KAAAAAAAMA",
            &[],
            &["Carbamyl (N-term)", "Oxidation (M)"],
            2,
        ),
        (
            "ACAACAACA",
            &[],
            &["Glutathione (C)", "Carbamidomethyl (C)"],
            1,
        ),
        ("AAAAAAAAA", &["Carbamidomethyl (C)"], &["Oxidation (M)"], 2),
        ("ACMACMACA", &["Carbamidomethyl (C)"], &["Oxidation (M)"], 3),
    ];
    for (seq, fixed, variable, max) in cases {
        let mpg = run_mpg(seq, fixed, variable, max)?;
        let fi = run_fi(seq, fixed, variable, max)?;
        assert_eq!(fi.len(), mpg.len(), "{seq}");
        for (a, b) in fi.iter().zip(&mpg) {
            // TEST_REAL_SIMILAR: relative 1e-5.
            assert!(
                ((a.0 - b.0) / b.0).abs() <= 1e-5,
                "{seq}: {} vs {}",
                a.0,
                b.0
            );
            assert_eq!(a.1, b.1, "{seq}");
        }
    }
    Ok(())
}

// FragmentIndex_test.cpp:849-875
#[test]
fn symmetric_window_self_hit() -> Result<()> {
    let entries = [entry("TEST1", "PEPTIDER")];
    let mut p = FragmentIndexParameters::default();
    p.min_ion_index = 0;
    p.min_isotope_error = 0;
    p.max_isotope_error = 0;
    let mut fi = index(p);
    fi.build(&entries)?;
    assert!(test_query(&fi, 2, true, &entries)?);
    Ok(())
}

fn peptider_spectrum(shift_ppm: f64) -> Result<MSSpectrum> {
    let seq = AASequence::parse("PEPTIDER")?;
    let ions = TheoreticalSpectrumGenerator::default().generate(&seq, 1, 1, None)?;
    Ok(MSSpectrum {
        ms_level: 2,
        peaks: ions.peaks,
        precursors: vec![Precursor {
            mz: seq.mz(2)? * (1.0 + shift_ppm * 1e-6),
            charge: 2,
            ..Precursor::default()
        }],
        ..MSSpectrum::default()
    })
}

// FragmentIndex_test.cpp:878-934
#[test]
fn asymmetric_window_compensates_calibration_offset() -> Result<()> {
    let entries = [entry("TEST", "PEPTIDER")];
    let mut p = FragmentIndexParameters::default();
    p.min_ion_index = 0;
    p.min_isotope_error = 0;
    p.max_isotope_error = 0;
    p.precursor_mass_tolerance_lower = 5.0;
    p.precursor_mass_tolerance_upper = 5.0;
    let mut fi = index(p.clone());
    fi.build(&entries)?;
    let spec = peptider_spectrum(8.0)?;
    let mut tight = SpectrumMatchesTopN::default();
    fi.query_spectrum(&spec, &mut tight)?;
    assert!(tight.hits.is_empty());
    p.precursor_mass_tolerance_lower = 15.0;
    fi.set_parameters(p)?;
    let mut asym = SpectrumMatchesTopN::default();
    fi.query_spectrum(&spec, &mut asym)?;
    assert!(!asym.hits.is_empty());
    Ok(())
}

// FragmentIndex_test.cpp:938-951
#[test]
fn open_search_mode_threshold() {
    assert!(FragmentIndex::is_open_search_mode(500.0, 1500.0, true));
    assert!(!FragmentIndex::is_open_search_mode(999.0, 999.0, true));
    assert!(!FragmentIndex::is_open_search_mode(1000.0, 1000.0, true));
    assert!(FragmentIndex::is_open_search_mode(1000.0001, 1000.0, true));
    assert!(!FragmentIndex::is_open_search_mode(0.9, 0.9, false));
    assert!(!FragmentIndex::is_open_search_mode(1.0, 1.0, false));
    assert!(FragmentIndex::is_open_search_mode(1.1, 0.5, false));
}

// FragmentIndex_test.cpp:955-1007
#[test]
fn open_mode_collapses_isotope_iteration() -> Result<()> {
    let entries = [entry("TEST", "PEPTIDER")];
    let run = |lo: i16, hi: i16| -> Result<usize> {
        let mut p = FragmentIndexParameters::default();
        p.precursor_mass_tolerance_lower = 0.5;
        p.precursor_mass_tolerance_upper = 1.5;
        p.precursor_mass_tolerance_unit_ppm = false;
        p.min_isotope_error = lo;
        p.max_isotope_error = hi;
        let mut fi = index(p);
        fi.build(&entries)?;
        let mut sms = SpectrumMatchesTopN::default();
        fi.query_spectrum(&peptider_spectrum(0.0)?, &mut sms)?;
        Ok(sms.hits.len())
    };
    assert_eq!(run(-2, 2)?, run(0, 0)?);
    Ok(())
}

// FragmentIndex_test.cpp:1010-1041
#[test]
fn asymmetric_closed_window_with_isotope_errors() -> Result<()> {
    let entries = [
        entry("P1", "PEPTIDER"),
        entry("P2", "PEPTIDERG"),
        entry("P3", "PEPTIDERA"),
    ];
    let mut p = FragmentIndexParameters::default();
    p.enzyme = "no cleavage".into();
    p.precursor_mass_tolerance_lower = 5.0;
    p.precursor_mass_tolerance_upper = 15.0;
    p.min_isotope_error = -1;
    p.max_isotope_error = 2;
    p.min_ion_index = 0;
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 50000;
    let mut fi = index(p);
    fi.build(&entries)?;
    assert!(test_query(&fi, 2, true, &entries)?);
    Ok(())
}

// FragmentIndex_test.cpp:1045-1086 (validation)
#[test]
fn invalid_windows_are_rejected() -> Result<()> {
    let with = |lower: f64, upper: f64| {
        let mut p = FragmentIndexParameters::default();
        p.precursor_mass_tolerance_lower = lower;
        p.precursor_mass_tolerance_upper = upper;
        p
    };
    assert!(FragmentIndex::new(with(-5.0, 10.0)).is_err());
    assert!(FragmentIndex::new(with(0.0, 0.0)).is_err());
    assert!(FragmentIndex::new(with(10.0, f64::NAN)).is_err());
    // Through the Param path, as the source's setParameters.
    let mut param = FragmentIndexParameters::defaults()?;
    param.set_value(
        "precursor:mass_tolerance_lower",
        openms::param::ParamValue::Float(-5.0),
        "",
        &[],
    )?;
    assert!(FragmentIndex::from_param(&param).is_err());
    Ok(())
}

// FragmentIndex_test.cpp:1089-1172
#[test]
fn peptides_in_mass_window_is_half_open() -> Result<()> {
    let mut p = no_mods(FragmentIndexParameters::default());
    p.precursor_mass_tolerance_unit_ppm = false;
    p.precursor_mass_tolerance_lower = 1000.0;
    p.precursor_mass_tolerance_upper = 1000.0;
    p.enzyme = "no cleavage".into();
    p.peptide_min_length = 1;
    let mut fi = index(p);
    let entries = [
        entry("P1", "AAAK"),
        entry("P2", "AAAAK"),
        entry("P3", "AAAAAK"),
        entry("P4", "AAAAAAK"),
        entry("P5", "AAAAAAAK"),
    ];
    fi.build(&entries)?;
    let peptides: Vec<Peptide> = fi.peptides().to_vec();
    assert_eq!(peptides.len(), 5);
    assert!(peptides[0].precursor_mz < peptides[1].precursor_mz);
    assert!(peptides[3].precursor_mz < peptides[4].precursor_mz);
    let p2 = peptides[2].precursor_mz;
    let r1 = fi.get_peptides_in_mass_window(p2, (-10.0, 10.0));
    assert_eq!(r1, (2, 3));
    assert!(peptides[r1.1].precursor_mz > p2 + 10.0);
    assert_eq!(
        fi.get_peptides_in_mass_window(peptides[4].precursor_mz, (-10.0, 10.0)),
        (4, 5)
    );
    assert_eq!(
        fi.get_peptides_in_mass_window(peptides[0].precursor_mz, (-10.0, 10.0)),
        (0, 1)
    );
    let gap = (peptides[1].precursor_mz + peptides[2].precursor_mz) * 0.5;
    let gap_tol = (peptides[2].precursor_mz - peptides[1].precursor_mz) * 0.25;
    let r4 = fi.get_peptides_in_mass_window(gap, (-gap_tol, gap_tol));
    assert_eq!(r4.0, r4.1);
    assert_eq!(fi.get_peptides_in_mass_window(100.0, (-10.0, 10.0)), (0, 0));
    assert_eq!(
        fi.get_peptides_in_mass_window(10000.0, (-10.0, 10.0)),
        (5, 5)
    );
    assert_eq!(
        fi.get_peptides_in_mass_window(p2, (-10000.0, 10000.0)),
        (0, 5)
    );
    let (lo, hi) = (peptides[1].precursor_mz, peptides[2].precursor_mz);
    assert_eq!(fi.get_peptides_in_mass_window(lo, (0.0, hi - lo)), (1, 3));
    // Reversed window (native check of the source's defensive branch).
    assert_eq!(fi.get_peptides_in_mass_window(p2, (10.0, -10.0)), (0, 0));
    Ok(())
}

// FragmentIndex_test.cpp:1175-1261
#[test]
fn range_second_is_never_scored() -> Result<()> {
    let mut p = FragmentIndexParameters::default();
    p.precursor_mass_tolerance_unit_ppm = false;
    p.precursor_mass_tolerance_lower = 5.0;
    p.precursor_mass_tolerance_upper = 5.0;
    p.min_isotope_error = 0;
    p.max_isotope_error = 0;
    p.enzyme = "no cleavage".into();
    p.peptide_min_length = 1;
    p.min_ion_index = 0;
    p.fragment_min_mz = 0;
    p.fragment_max_mz = 50000;
    let mut fi = index(p);
    let entries = [entry("A", "PEPTIDER"), entry("B", "PEPTIDERA")];
    fi.build(&entries)?;
    let peptides = fi.peptides().to_vec();
    assert_eq!(peptides.len(), 2);
    assert!(peptides[0].precursor_mz < peptides[1].precursor_mz);
    assert!(peptides[1].precursor_mz > peptides[0].precursor_mz + 5.0);
    assert_eq!(
        fi.get_peptides_in_mass_window(peptides[0].precursor_mz, (-5.0, 5.0)),
        (0, 1)
    );
    let generator = TheoreticalSpectrumGenerator {
        add_first_prefix_ion: true,
        ..TheoreticalSpectrumGenerator::default()
    };
    let theo_b = generator.generate(&AASequence::parse("PEPTIDERA")?, 1, 1, None)?;
    let spec = MSSpectrum {
        ms_level: 2,
        peaks: theo_b.peaks,
        precursors: vec![Precursor {
            mz: f64::from(peptides[0].precursor_mz),
            charge: 1,
            ..Precursor::default()
        }],
        ..MSSpectrum::default()
    };
    let mut sms = SpectrumMatchesTopN::default();
    fi.query_spectrum(&spec, &mut sms)?;
    assert!(
        !sms.hits
            .iter()
            .any(|h| h.peptide_idx == 1 && h.num_matched > 0)
    );
    Ok(())
}

fn snes_params(fixed: &[&str]) -> FragmentIndexParameters {
    let mut p = no_mods(FragmentIndexParameters::default());
    p.enzyme_specificity = DigestionSpecificity::None;
    p.peptide_min_length = 8;
    p.peptide_max_length = 12;
    p.peptide_min_mass = 0;
    p.peptide_max_mass = 50000;
    p.modifications_fixed = fixed.iter().map(|s| s.to_string()).collect();
    p.snes_enabled = true;
    p
}

const SNES_PROTEIN: &str = "AKACDEFGRHILMNPQSTV";

// FragmentIndex_test.cpp:1267-1304
#[test]
fn snes_mother_enumeration() -> Result<()> {
    let mut fi = index(snes_params(&[]));
    fi.build(&[entry("p", SNES_PROTEIN)])?;
    assert!(fi.is_snes_mode());
    assert_eq!(fi.peptides().len(), 24);
    let c = fi
        .peptides()
        .iter()
        .filter(|p| FragmentIndex::is_single_c_mother(p.mod_bitmask))
        .count();
    assert_eq!(c, 12);
    assert_eq!(fi.peptides().len() - c, 12);
    Ok(())
}

// FragmentIndex_test.cpp:1307-1330
#[test]
fn snes_is_skipped_for_specific_digestion() -> Result<()> {
    let mut p = specificity_params(DigestionSpecificity::Full, 2, 100);
    p.snes_enabled = true;
    let mut fi = index(p);
    fi.build(&[entry("p", SNES_PROTEIN)])?;
    assert!(!fi.is_snes_mode());
    assert_eq!(fi.peptides().len(), 3);
    Ok(())
}

fn find_mother(fi: &FragmentIndex, single_c: bool, pred: impl Fn(&Peptide) -> bool) -> Peptide {
    *fi.peptides()
        .iter()
        .find(|p| {
            FragmentIndex::is_single_c_mother(p.mod_bitmask) == single_c
                && p.protein_idx == 0
                && pred(p)
        })
        .expect("mother present")
}

// FragmentIndex_test.cpp:1333-1379
#[test]
fn realize_snes_length_single_n() -> Result<()> {
    let entries = [entry("p", SNES_PROTEIN)];
    let mut fi = index(snes_params(&[]));
    fi.build(&entries)?;
    let target = AASequence::parse("AKACDEFGRH")?.mono_mass()? + PROTON_MASS_U;
    let mother = find_mother(&fi, false, |p| p.sequence.0 == 0);
    let realized = fi.realize_snes_length(&mother, &entries, target, 10.0, 10.0, true)?;
    assert_eq!(realized, Some(10));
    let seq = fi.reconstruct_realized_sub_sequence(&mother, &entries, 10, 0)?;
    assert_eq!(seq.as_str(), "AKACDEFGRH");
    Ok(())
}

// FragmentIndex_test.cpp:1382-1442
#[test]
fn realize_snes_length_single_c() -> Result<()> {
    let entries = [entry("p", SNES_PROTEIN)];
    let mut fi = index(snes_params(&[]));
    fi.build(&entries)?;
    let target = AASequence::parse("ILMNPQSTV")?.mono_mass()? + PROTON_MASS_U;
    let mother = find_mother(&fi, true, |p| {
        usize::from(p.sequence.0) + usize::from(p.sequence.1) == 19
    });
    assert_eq!(
        fi.realize_snes_length(&mother, &entries, target, 10.0, 10.0, true)?,
        Some(9)
    );
    let shifted = target + 0.05;
    assert_eq!(
        fi.realize_snes_length(&mother, &entries, shifted, 500.0, 10.0, true)?,
        Some(9)
    );
    assert_eq!(
        fi.realize_snes_length(&mother, &entries, shifted, 10.0, 500.0, true)?,
        None
    );
    let seq = fi.reconstruct_realized_sub_sequence(&mother, &entries, 9, 0)?;
    assert_eq!(seq.as_str(), "ILMNPQSTV");
    Ok(())
}

// FragmentIndex_test.cpp:1445-1494
#[test]
fn snes_index_is_smaller_than_naive() -> Result<()> {
    let entries = [entry("p", SNES_PROTEIN)];
    let mut naive_params = snes_params(&[]);
    naive_params.snes_enabled = false;
    let mut naive = index(naive_params);
    naive.build(&entries)?;
    let mut snes = index(snes_params(&[]));
    snes.build(&entries)?;
    assert!(!naive.is_snes_mode());
    assert!(snes.is_snes_mode());
    assert!(snes.peptides().len() < naive.peptides().len());
    assert!(!naive.peptides().is_empty() && !snes.peptides().is_empty());
    // Native: fewer fragments too, the point of one series per mother.
    assert!(snes.num_fragments() < naive.num_fragments());
    Ok(())
}

// FragmentIndex_test.cpp:1497-1541
#[test]
fn realized_sub_sequence_applies_fixed_modifications() -> Result<()> {
    let entries = [entry("p", SNES_PROTEIN)];
    let mut fi = index(snes_params(&["Carbamidomethyl (C)"]));
    fi.build(&entries)?;
    let mother = find_mother(&fi, false, |p| p.sequence.0 == 0);
    let seq = fi.reconstruct_realized_sub_sequence(&mother, &entries, 8, 0)?;
    assert_eq!(seq.as_str(), "AKACDEFG");
    assert_eq!(seq.len(), 8);
    assert_eq!(seq.to_string(), "AKAC(Carbamidomethyl)DEFG");
    Ok(())
}

// Tier 4: the SNES spectrum query is a documented gap and says so.
#[test]
fn snes_query_is_reported_unsupported() -> Result<()> {
    let entries = [entry("p", SNES_PROTEIN)];
    let mut fi = index(snes_params(&[]));
    fi.build(&entries)?;
    let spec = spectrum(&AASequence::parse("AKACDEFGRH")?, 1, 1100.0, true)?;
    let mut sms = SpectrumMatchesTopN::default();
    assert!(fi.query_spectrum(&spec, &mut sms).is_err());
    Ok(())
}

// Tier 4: parameters round-trip and keep the source defaults.
#[test]
fn parameters_round_trip() -> Result<()> {
    let defaults = FragmentIndexParameters::defaults()?;
    let parsed = FragmentIndexParameters::from_param(&defaults)?;
    assert_eq!(parsed, FragmentIndexParameters::default());
    assert_eq!(
        defaults.description("fragment:mass_tolerance_unit")?,
        "Unit of fragment m"
    );
    let mut changed = parsed.clone();
    changed.snes_enabled = true;
    changed.enzyme_specificity = DigestionSpecificity::None;
    changed.min_isotope_error = -1;
    changed.fragment_mz_tolerance_unit_ppm = false;
    assert_eq!(
        FragmentIndexParameters::from_param(&changed.to_param()?)?,
        changed
    );
    // SNES without y ions is refused (source updateMembers_).
    changed.add_y_ions = false;
    assert!(FragmentIndex::new(changed).is_err());
    Ok(())
}

// Tier 4: Σ sets and slot order follow the source rules.
#[test]
fn sigma_delta_sets_and_mod_slots() -> Result<()> {
    let mut p = FragmentIndexParameters::default();
    p.modifications_variable = vec!["Oxidation (M)".into(), "Acetyl (Protein N-term)".into()];
    p.modifications_max_variable_mods_per_peptide = 2;
    let fi = index(p);
    let ox = 15.994_915;
    let base = fi.snes_sigma_delta_set(SnesAnchor::None);
    assert_eq!(base.len(), 3);
    assert!((base[1] - ox).abs() < 1e-5 && (base[2] - 2.0 * ox).abs() < 1e-5);
    assert!(fi.snes_sigma_delta_set(SnesAnchor::ProtNterm).len() > base.len());
    assert_eq!(fi.snes_sigma_delta_set(SnesAnchor::ProtCterm), base);
    let slots = fi.build_mod_slots(b"MAMK", true, false);
    assert_eq!(slots.len(), 3);
    assert_eq!(
        slots[0].position,
        openms::analysis::id::fragment_index::ModSlot::NTERM_SLOT
    );
    assert_eq!((slots[1].position, slots[2].position), (0, 2));
    assert_eq!(fi.build_mod_slots(b"MAMK", false, false).len(), 2);
    Ok(())
}

// Tier 4: the build is identical for every thread count.
#[test]
fn build_is_identical_across_thread_counts() -> Result<()> {
    let entries = [entry("test1", PROTEIN), entry("t", "ARGEPADSSRKDFDMDMDM")];
    let build = |threads| -> Result<(Vec<Peptide>, Vec<(u32, f32)>)> {
        let mut fi = index(FragmentIndexParameters::default());
        fi.threads = threads;
        fi.build(&entries)?;
        Ok((
            fi.peptides().to_vec(),
            fi.fragments()
                .iter()
                .map(|f| (f.peptide_idx, f.fragment_mz))
                .collect(),
        ))
    };
    let serial = build(openms::concept::parallel::Threads::serial())?;
    assert_eq!(
        build(openms::concept::parallel::Threads::from_cli(4))?,
        serial
    );
    Ok(())
}

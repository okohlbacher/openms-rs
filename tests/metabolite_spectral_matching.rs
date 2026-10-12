// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ANALYSIS/ID/MetaboliteSpectralMatching` against the upstream class test
//! (`MetaboliteSpectralMatching_test.cpp` at core `bc9cc12`), whose library,
//! parameters and literals are transcribed here, plus hyperscore values
//! derived by hand from the source formula and native error cases.

use openms::analysis::id::metabolite_spectral_matching::{
    IonizationMode, MassErrorUnit, MetaboliteSpectralMatching, ReportMode, SpectralMatch,
};
use openms::format::mztab::{MzTab, MzTabSmallMoleculeSectionRow};
use openms::identification::PeakAnnotation;
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D, Precursor};
use openms::metadata::{DriftTimeUnit, ImTypes, MetaValue};
use openms::param::{Param, ParamValue};
use openms::{Error, Result};

fn peaks(points: &[(f64, f32)]) -> Vec<Peak1D> {
    points
        .iter()
        .map(|&(mz, intensity)| Peak1D { mz, intensity })
        .collect()
}

fn precursor(mz: f64, charge: i32) -> Precursor {
    Precursor {
        mz,
        charge,
        ..Precursor::default()
    }
}

/// The class test's `makeLibrary`.
fn library() -> MSExperiment {
    let entry = |name: &str,
                 ccs: &str,
                 hmdb: &str,
                 formula: &str,
                 inchi: &str,
                 smiles: &str,
                 mz: f64,
                 rt: f64| {
        let mut s = MSSpectrum {
            name: name.into(),
            rt,
            precursors: vec![precursor(mz, 1)],
            peaks: peaks(&[
                (100.0, 999.0),
                (150.0, 500.0),
                (200.0, 300.0),
                (250.0, 100.0),
            ]),
            ..MSSpectrum::default()
        };
        for (key, value) in [
            ("Metabolite_Name", name),
            ("CCS", ccs),
            ("HMDB_ID", hmdb),
            ("Sum_Formula", formula),
            ("Inchi_String", inchi),
            ("SMILES_String", smiles),
            ("Precursor_Ion", "[M+H]+"),
        ] {
            s.metadata.insert(key.into(), MetaValue::from(value));
        }
        s
    };
    MSExperiment {
        spectra: vec![
            entry(
                "compound_A",
                "150.0",
                "HMDB0000001",
                "C10H15N",
                "InChI=1S/compoundA",
                "CCO",
                300.0,
                0.0,
            ),
            entry(
                "compound_B",
                "200.0",
                "HMDB0000002",
                "C12H17N",
                "InChI=1S/compoundB",
                "CCCO",
                300.1,
                1.0,
            ),
        ],
        ..MSExperiment::default()
    }
}

/// The class test's `makeParams`.
fn params(ccs_error_percent: f64) -> Param {
    let mut p = MetaboliteSpectralMatching::defaults().unwrap();
    for (key, value) in [
        ("prec_mass_error_value", ParamValue::Float(1.0)),
        ("frag_mass_error_value", ParamValue::Float(0.5)),
        ("mass_error_unit", ParamValue::String("Da".into())),
        ("ionization_mode", ParamValue::String("positive".into())),
        ("report_mode", ParamValue::String("top3".into())),
        ("merge_spectra", ParamValue::String("false".into())),
        ("ccs_error_percent", ParamValue::Float(ccs_error_percent)),
    ] {
        p.set_value(key, value, "", &[]).unwrap();
    }
    p
}

/// The class test's experimental spectrum, without CCS information.
fn experiment() -> MSSpectrum {
    MSSpectrum {
        ms_level: 2,
        rt: 10.0,
        precursors: vec![precursor(300.05, 1)],
        peaks: peaks(&[
            (100.0, 800.0),
            (150.0, 400.0),
            (200.0, 250.0),
            (250.0, 80.0),
        ]),
        ..MSSpectrum::default()
    }
}

fn run(spectrum: MSSpectrum, ccs_error_percent: f64) -> Result<(MzTab, Vec<String>)> {
    let msm = MetaboliteSpectralMatching::from_param(&params(ccs_error_percent))?;
    let mut msexp = MSExperiment {
        spectra: vec![spectrum],
        ..MSExperiment::default()
    };
    let mut spec_db = library();
    let mut mztab = MzTab::default();
    let warnings = msm.run(&mut msexp, &mut spec_db, &mut mztab, "")?;
    Ok((mztab, warnings))
}

fn opt(row: &MzTabSmallMoleculeSectionRow, name: &str) -> String {
    row.opt
        .iter()
        .find(|entry| entry.name == name)
        .map(|entry| entry.value.to_cell_string())
        .unwrap_or_else(|| panic!("missing optional column {name}"))
}

/// Class test: `SpectralMatch CCS getters and setters`.
#[test]
fn spectral_match_ccs_fields() {
    let mut sm = SpectralMatch::default();
    assert_eq!(sm.observed_ccs, ImTypes::DRIFTTIME_NOT_SET);
    assert_eq!(sm.found_ccs, ImTypes::DRIFTTIME_NOT_SET);
    sm.observed_ccs = 167.3;
    sm.found_ccs = 170.1;
    let copy = sm.clone();
    assert_eq!(copy.observed_ccs, 167.3);
    assert_eq!(copy.found_ccs, 170.1);
}

/// Class test: `MetaboliteSpectralMatching default parameters`, plus the full
/// default set from the source constructor.
#[test]
fn default_parameters() -> Result<()> {
    let p = MetaboliteSpectralMatching::defaults()?;
    assert_eq!(p.value("ccs_error_percent")?, &ParamValue::Float(0.0));
    assert_eq!(p.value("prec_mass_error_value")?, &ParamValue::Float(100.0));
    assert_eq!(p.value("frag_mass_error_value")?, &ParamValue::Float(500.0));
    assert_eq!(p.value("mass_error_unit")?.as_str()?, "ppm");
    assert_eq!(p.value("report_mode")?.as_str()?, "top3");
    assert_eq!(p.value("ionization_mode")?.as_str()?, "positive");
    assert_eq!(p.value("merge_spectra")?.as_str()?, "true");
    assert_eq!(p.entry("ccs_error_percent")?.min_float, 0.0);
    let msm = MetaboliteSpectralMatching::from_param(&params(0.0))?;
    assert_eq!(msm.mass_error_unit, MassErrorUnit::Da);
    assert_eq!(msm.report_mode, ReportMode::Top3);
    assert_eq!(msm.ionization_mode, IonizationMode::Positive);
    assert!(!msm.merge_spectra);
    Ok(())
}

/// Class test: `CCS filtering disabled reports all matches with CCS columns`.
#[test]
fn ccs_filtering_disabled_reports_all_matches() -> Result<()> {
    let mut spectrum = experiment();
    spectrum
        .metadata
        .insert("CCS".into(), MetaValue::from("152.0"));
    let (mztab, _) = run(spectrum, 0.0)?;
    let rows = &mztab.small_molecule_data;
    assert_eq!(rows.len(), 2);
    for row in rows {
        assert_eq!(opt(row, "opt_observed_ccs"), "152.0");
        let library = opt(row, "opt_library_ccs");
        assert!(library == "150.0" || library == "200.0");
        assert_ne!(opt(row, "opt_ccs_error_percent"), "null");
    }
    Ok(())
}

/// Class test: `CCS filtering enabled keeps only within-tolerance matches`.
#[test]
fn ccs_filtering_enabled_keeps_close_matches() -> Result<()> {
    let mut spectrum = experiment();
    spectrum
        .metadata
        .insert("CCS".into(), MetaValue::from("152.0"));
    let (mztab, warnings) = run(spectrum, 5.0)?;
    let rows = &mztab.small_molecule_data;
    assert_eq!(rows.len(), 1);
    assert_eq!(opt(&rows[0], "opt_library_ccs"), "150.0");
    assert_eq!(opt(&rows[0], "opt_observed_ccs"), "152.0");
    assert!(warnings.is_empty());
    Ok(())
}

/// Class test: `CCS filtering converts 1/K0 (VSSC) drift time to CCS`.
#[test]
fn ccs_filtering_converts_inverse_mobility() -> Result<()> {
    let mut spectrum = experiment();
    spectrum.drift_time = 0.726;
    spectrum.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    let (mztab, _) = run(spectrum, 5.0)?;
    let rows = &mztab.small_molecule_data;
    assert_eq!(rows.len(), 1);
    assert_eq!(opt(&rows[0], "opt_library_ccs"), "150.0");
    assert_ne!(opt(&rows[0], "opt_observed_ccs"), "null");
    Ok(())
}

/// Class test: `CCS filtering keeps matches when CCS is missing`; the warning
/// text is the source's `OPENMS_LOG_WARN` line.
#[test]
fn ccs_filtering_keeps_matches_without_ccs() -> Result<()> {
    let (mztab, warnings) = run(experiment(), 5.0)?;
    let rows = &mztab.small_molecule_data;
    assert_eq!(rows.len(), 2);
    for row in rows {
        assert_eq!(opt(row, "opt_observed_ccs"), "null");
        assert_eq!(opt(row, "opt_ccs_error_percent"), "null");
    }
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].starts_with("Warning: 'ccs_error_percent' is set to 5 but"));
    Ok(())
}

/// The remaining small-molecule cells, derived from `exportMzTab_` for the
/// class test's first scenario: identifier (Metabolite_Name, since no GNPS or
/// MassBank ID), formula, SMILES, InChI, library precursor m/z and charge, RT
/// floored to 0.1 s, the fixed database fields, and the ppm error
/// `floor((300 - 300.05) / 300 * 1e6 * 100) / 100`, the double nearest -166.67.
#[test]
fn exported_rows_follow_the_source() -> Result<()> {
    let (mztab, _) = run(experiment(), 0.0)?;
    let rows = &mztab.small_molecule_data;
    let a = rows
        .iter()
        .find(|r| r.description.get() == "compound_A")
        .expect("compound_A row");
    assert_eq!(a.identifier.get()[0].get(), "compound_A");
    assert_eq!(a.chemical_formula.get(), "C10H15N");
    assert_eq!(a.smiles.get(), "CCO");
    assert_eq!(a.inchi_key.get(), "InChI=1S/compoundA");
    assert_eq!(a.exp_mass_to_charge.to_cell_string(), "300.0");
    assert_eq!(a.charge.to_cell_string(), "1");
    assert_eq!(a.retention_time.get()[0].to_cell_string(), "10.0");
    assert_eq!(a.database.get(), "MassBank");
    assert_eq!(a.database_version.get(), "Sep 27, 2013");
    // StringUtils::toStr(double): fixed with 15 fractional digits, trailing
    // zeros trimmed (NumericFormatting.h appendNumeric).
    assert_eq!(opt(a, "opt_ppm_error"), "-166.669999999999987");
    assert_eq!(opt(a, "opt_adduct_ion"), "[M+H]+");
    assert_eq!(opt(a, "opt_sec_id"), "HMDB0000001");
    assert_eq!(opt(a, "opt_source_idx"), "0");
    let names: Vec<_> = a.opt.iter().map(|e| e.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "opt_ppm_error",
            "opt_adduct_ion",
            "opt_match_score",
            "opt_sec_id",
            "opt_source_idx",
            "opt_spec_native_id",
            "opt_observed_ccs",
            "opt_library_ccs",
            "opt_ccs_error_percent"
        ]
    );
    // Both library entries share their peaks, so both score
    // ln(800*999 + 400*500 + 250*300 + 80*100) + ln(4!) = ln(1_082_200) + ln(24).
    let expected = (1_082_200.0_f64).ln() + (24.0_f64).ln();
    for row in rows {
        // std::ostream default precision 6 (%g): 17.07256... prints as 17.0726.
        assert!((expected - 17.072_560_394_538_72).abs() < 1e-12);
        assert_eq!(opt(row, "opt_match_score"), "17.0726");
    }
    Ok(())
}

/// Class test: `computeHyperScore`, plus the value from the source formula.
#[test]
fn hyper_score_matches_class_test() -> Result<()> {
    let exp = MSSpectrum {
        peaks: peaks(&[
            (100.0, 100.0),
            (200.0, 200.0),
            (300.0, 300.0),
            (400.0, 50.0),
        ]),
        ..MSSpectrum::default()
    };
    let db = MSSpectrum {
        peaks: peaks(&[
            (100.0, 999.0),
            (200.0, 999.0),
            (300.0, 999.0),
            (400.0, 999.0),
        ]),
        ..MSSpectrum::default()
    };
    let score = MetaboliteSpectralMatching::compute_hyper_score(0.5, false, &exp, &db, 0.0)?;
    assert!(score > 0.0);
    let expected = (999.0_f64 * 650.0).ln() + (24.0_f64).ln();
    assert!((score - expected).abs() <= 1e-12 * expected);
    let empty = MSSpectrum::default();
    assert_eq!(
        MetaboliteSpectralMatching::compute_hyper_score(0.5, false, &empty, &db, 0.0)?,
        0.0
    );
    // Fewer than three matched peaks score 0.
    let two = MSSpectrum {
        peaks: peaks(&[(100.0, 100.0), (200.0, 200.0)]),
        ..MSSpectrum::default()
    };
    assert_eq!(
        MetaboliteSpectralMatching::compute_hyper_score(0.5, false, &two, &db, 0.0)?,
        0.0
    );
    Ok(())
}

#[test]
fn hyper_score_annotations() -> Result<()> {
    let exp = MSSpectrum {
        peaks: peaks(&[(100.0, 100.0), (200.0, 200.0), (300.0, 300.0)]),
        ..MSSpectrum::default()
    };
    let mut db = MSSpectrum {
        peaks: peaks(&[(100.1, 5.0), (200.0, 6.0), (300.0, 7.0)]),
        ..MSSpectrum::default()
    };
    let mut annotations: Vec<PeakAnnotation> = Vec::new();
    // Without data arrays no annotations are produced.
    MetaboliteSpectralMatching::compute_hyper_score_with_annotations(
        0.5,
        false,
        &exp,
        &db,
        &mut annotations,
        0.0,
    )?;
    assert!(annotations.is_empty());
    db.string_data_arrays.push(DataArray {
        data: vec!["a".into(), "b".into(), "c".into()],
        ..DataArray::default()
    });
    db.integer_data_arrays.push(DataArray {
        data: vec![1, 2, 3],
        ..DataArray::default()
    });
    MetaboliteSpectralMatching::compute_hyper_score_with_annotations(
        0.5,
        false,
        &exp,
        &db,
        &mut annotations,
        0.0,
    )?;
    assert_eq!(annotations.len(), 3);
    assert_eq!(annotations[0].annotation, "a");
    assert_eq!(annotations[0].charge, 1);
    assert_eq!(annotations[0].mz, 100.0); // the experimental peak's m/z
    assert_eq!(annotations[2].intensity, 300.0);
    db.integer_data_arrays[0].data.pop();
    assert!(matches!(
        MetaboliteSpectralMatching::compute_hyper_score_with_annotations(
            0.5,
            false,
            &exp,
            &db,
            &mut annotations,
            0.0
        ),
        Err(Error::InvalidValue(_))
    ));
    Ok(())
}

#[test]
fn run_refuses_merging_and_bad_libraries() -> Result<()> {
    let mut msexp = MSExperiment {
        spectra: vec![experiment()],
        ..MSExperiment::default()
    };
    let mut mztab = MzTab::default();
    let merging = MetaboliteSpectralMatching::default();
    assert!(matches!(
        merging.run(&mut msexp, &mut library(), &mut mztab, ""),
        Err(Error::Unsupported(_))
    ));
    let msm = MetaboliteSpectralMatching::from_param(&params(0.0))?;
    let mut no_precursor = library();
    no_precursor.spectra[1].precursors.clear();
    assert!(matches!(
        msm.run(&mut msexp, &mut no_precursor, &mut mztab, ""),
        Err(Error::InvalidValue(_))
    ));
    let mut no_hmdb = library();
    no_hmdb.spectra[0].metadata.remove("HMDB_ID");
    assert!(matches!(
        msm.run(&mut msexp, &mut no_hmdb, &mut mztab, ""),
        Err(Error::InvalidValue(_))
    ));
    assert!(mztab.small_molecule_data.is_empty());
    Ok(())
}

/// `best` keeps one row; negative mode skips the positively charged library.
#[test]
fn report_and_ionization_modes() -> Result<()> {
    let mut p = params(0.0);
    p.set_value("report_mode", ParamValue::String("best".into()), "", &[])?;
    let msm = MetaboliteSpectralMatching::from_param(&p)?;
    let mut msexp = MSExperiment {
        spectra: vec![experiment()],
        ..MSExperiment::default()
    };
    let mut mztab = MzTab::default();
    msm.run(&mut msexp, &mut library(), &mut mztab, "")?;
    assert_eq!(mztab.small_molecule_data.len(), 1);
    // Equal scores keep library (precursor m/z) order: compound_A first.
    assert_eq!(mztab.small_molecule_data[0].description.get(), "compound_A");
    p.set_value(
        "ionization_mode",
        ParamValue::String("negative".into()),
        "",
        &[],
    )?;
    let msm = MetaboliteSpectralMatching::from_param(&p)?;
    msm.run(&mut msexp, &mut library(), &mut mztab, "")?;
    assert!(mztab.small_molecule_data.is_empty());
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FORMAT/DATAACCESS/SiriusFragmentAnnotation` against the upstream class
//! test (`SiriusFragmentAnnotation_test.cpp` at core `bc9cc12`) and its
//! project-space fixture `SiriusFragmentAnnotation_test/`, plus values read
//! off that fixture by hand and native error cases.

use openms::format::sirius_fragment_annotation::{
    extract_and_resolve_sirius_annotations, extract_annotations_from_sirius_file,
    extract_compound_ranking_and_filename, extract_compound_ranking_and_score,
    extract_concat_m_ids, extract_concat_native_ids, extract_feature_id,
};
use openms::kernel::MSSpectrum;
use openms::{Error, Result};
use std::path::{Path, PathBuf};

const NATIVE_IDS: &str = "sample=1 period=1 cycle=676 experiment=4|sample=1 period=1 cycle=677 experiment=5|sample=1 period=1 cycle=678 experiment=3";

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/SiriusFragmentAnnotation_test")
}

fn meta(spectrum: &MSSpectrum, key: &str) -> String {
    spectrum.metadata[key].to_string()
}

/// `TEST_REAL_SIMILAR`'s default: relative difference at most 1e-5.
fn similar(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 1e-5 * expected.abs().max(actual.abs())
}

/// Class test section 1: target, peak m/z from the `mz` column.
#[test]
fn target_annotations_match_class_test() -> Result<()> {
    let spectra = extract_annotations_from_sirius_file(workspace(), 1, false, false)?;
    let s = &spectra[0];
    assert_eq!(s.native_id, NATIVE_IDS);
    assert_eq!(s.ms_level, 2);
    assert!(!s.peaks.is_empty());
    assert!(similar(s.peaks[0].mz, 70.040098));
    assert_eq!(meta(s, "peak_mz"), "mz");
    assert_eq!(s.float_data_arrays[0].name, "exact_mass");
    assert!(similar(
        f64::from(s.float_data_arrays[0].data[0]),
        70.040098
    ));
    assert_eq!(s.string_data_arrays[0].data[0], "C2H3N3");
    assert_eq!(meta(s, "annotated_sumformula"), "C15H17ClN4");
    assert_eq!(meta(s, "annotated_adduct"), "[M+H]+");
    assert_eq!(s.metadata["decoy"].as_i64()?, 0);
    // Read off the fixture: first row of spectra/C15H17ClN4_[M+H]+.tsv is
    // `70.040098 24218.53 100.00 70.039974 C2H3N3 [M + H]+`, 22 rows in all.
    assert_eq!(s.peaks[0].mz, 70.040098);
    assert_eq!(s.peaks[0].intensity, 24218.53_f32);
    assert_eq!(
        s.float_data_arrays[0].data[0],
        "70.039974".parse::<f32>().unwrap()
    );
    assert_eq!(s.peaks.len(), 22);
    assert_eq!(s.string_data_arrays[1].data[0], "[M + H]+");
    assert_eq!(s.string_data_arrays[0].name, "explanation");
    assert_eq!(s.metadata["score"].as_f64()?, 0.9685426989556167);
    assert_eq!(meta(s, "feat_id"), "5053681185761206615");
    assert!(
        s.name.starts_with(
            "Myclobutanil_PestMix1_1ngSolventDDA20-50.wiff_sample=1 period=1 cycle=676"
        )
    );
    assert_eq!(s.name.matches('|').count(), 2);
    Ok(())
}

/// Class test section 2: target with exact masses as peak positions.
#[test]
fn exact_mass_annotations_match_class_test() -> Result<()> {
    let spectra = extract_annotations_from_sirius_file(workspace(), 1, false, true)?;
    let s = &spectra[0];
    assert_eq!(s.native_id, NATIVE_IDS);
    assert_eq!(s.ms_level, 2);
    assert!(!s.peaks.is_empty());
    // The class test compares with 70.040098 under TEST_REAL_SIMILAR; the
    // value read is the `exactmass` column, 70.039974, which is within 1e-5.
    assert!(similar(s.peaks[0].mz, 70.040098));
    assert_eq!(s.peaks[0].mz, 70.039974);
    assert_eq!(meta(s, "peak_mz"), "exact_mass");
    assert_eq!(s.float_data_arrays[0].name, "mz");
    assert!(similar(
        f64::from(s.float_data_arrays[0].data[0]),
        70.040098
    ));
    assert_eq!(
        s.float_data_arrays[0].data[0],
        "70.040098".parse::<f32>().unwrap()
    );
    assert_eq!(s.string_data_arrays[0].data[0], "C2H3N3");
    assert_eq!(meta(s, "annotated_sumformula"), "C15H17ClN4");
    assert_eq!(meta(s, "annotated_adduct"), "[M+H]+");
    assert_eq!(s.metadata["decoy"].as_i64()?, 0);
    Ok(())
}

/// Class test section 3: decoys.
#[test]
fn decoy_annotations_match_class_test() -> Result<()> {
    // use_exact_mass is forced off for decoys.
    for exact in [false, true] {
        let spectra = extract_annotations_from_sirius_file(workspace(), 1, true, exact)?;
        let s = &spectra[0];
        assert_eq!(s.native_id, NATIVE_IDS);
        assert_eq!(s.ms_level, 2);
        assert!(!s.peaks.is_empty());
        assert!(similar(s.peaks[0].mz, 53.013424));
        assert_eq!(meta(s, "peak_mz"), "mz");
        assert_eq!(s.string_data_arrays[0].data[0], "C2N2");
        assert_eq!(meta(s, "annotated_sumformula"), "C15H17ClN4");
        assert_eq!(meta(s, "annotated_adduct"), "[M+H]+");
        assert_eq!(s.metadata["decoy"].as_i64()?, 1);
        // Read off the fixture: `53.013424 100.00 C2N2 [M + H]+`; the
        // intensity is the rel.intensity column; no extra masses for decoys.
        assert_eq!(s.peaks[0].intensity, 100.0);
        assert!(s.float_data_arrays[0].data.is_empty());
        assert_eq!(s.peaks.len(), 22);
    }
    Ok(())
}

#[test]
fn spectrum_ms_and_candidate_helpers() -> Result<()> {
    assert_eq!(extract_concat_native_ids(workspace())?, NATIVE_IDS);
    assert_eq!(extract_concat_m_ids(workspace())?.split('|').count(), 3);
    assert_eq!(extract_feature_id(workspace())?, "5053681185761206615");
    let files = extract_compound_ranking_and_filename(workspace())?;
    assert_eq!(files.len(), 1);
    assert_eq!(files[&1], "C15H17ClN4_[M+H]+.tsv");
    let scores = extract_compound_ranking_and_score(workspace())?;
    assert_eq!(scores[&1], 0.9685426989556167);
    // Only one candidate exists, so a larger max_rank still yields one.
    assert_eq!(
        extract_annotations_from_sirius_file(workspace(), 5, false, false)?.len(),
        1
    );
    Ok(())
}

fn copy_dir(from: &Path, to: &Path) -> std::io::Result<()> {
    std::fs::create_dir_all(to)?;
    for entry in std::fs::read_dir(from)? {
        let entry = entry?;
        let target = to.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            copy_dir(&entry.path(), &target)?;
        } else {
            std::fs::copy(entry.path(), target)?;
        }
    }
    Ok(())
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("sirius-fa-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

#[test]
fn resolve_keeps_one_spectrum_per_native_id() -> Result<()> {
    let root = scratch("resolve");
    let (a, b) = (root.join("a"), root.join("b"));
    copy_dir(&workspace(), &a)?;
    copy_dir(&workspace(), &b)?;
    let resolved = extract_and_resolve_sirius_annotations(&[&a, &b], 0.5, false, true)?;
    assert_eq!(resolved.len(), 1);
    assert_eq!(resolved[0].target.native_id, NATIVE_IDS);
    assert_eq!(resolved[0].decoy.metadata["decoy"].as_i64()?, 1);
    // Without decoy generation the decoy only carries the identifiers.
    let resolved = extract_and_resolve_sirius_annotations(&[&a], 0.0, false, false)?;
    assert!(resolved[0].decoy.peaks.is_empty());
    assert_eq!(resolved[0].decoy.native_id, NATIVE_IDS);
    assert!(resolved[0].decoy.metadata.is_empty());
    // The fixture's score is 0.9685...: a higher threshold drops it.
    assert!(extract_and_resolve_sirius_annotations(&[&a], 0.99, false, true)?.is_empty());
    // A directory without spectra/ is skipped.
    let empty = root.join("empty");
    std::fs::create_dir_all(&empty)?;
    assert!(extract_and_resolve_sirius_annotations(&[&empty], 0.0, false, true)?.is_empty());
    std::fs::remove_dir_all(&root)?;
    Ok(())
}

#[test]
fn malformed_tables_are_errors() -> Result<()> {
    let root = scratch("bad");
    copy_dir(&workspace(), &root)?;
    // A fragment row too short for the exact-mass column (the source indexes
    // out of range).
    std::fs::write(
        root.join("spectra/C15H17ClN4_[M+H]+.tsv"),
        "mz\tintensity\n70.0\t1.0\n",
    )?;
    assert!(matches!(
        extract_annotations_from_sirius_file(&root, 1, false, false),
        Err(Error::InvalidValue(_))
    ));
    // A non-numeric rank (the source's ConversionError).
    std::fs::write(
        root.join("formula_candidates.tsv"),
        "formulaRank\tmolecularFormula\tadduct\texplainedIntensity\nfirst\tC\t[M+H]+\t0.5\n",
    )?;
    assert!(matches!(
        extract_compound_ranking_and_filename(&root),
        Err(Error::InvalidValue(_))
    ));
    // A gap in the ranks (the source's std::map::at throws).
    std::fs::write(
        root.join("formula_candidates.tsv"),
        "formulaRank\tmolecularFormula\tadduct\texplainedIntensity\n2\tC\t[M+H]+\t0.5\n",
    )?;
    assert!(matches!(
        extract_annotations_from_sirius_file(&root, 2, false, false),
        Err(Error::InvalidValue(_))
    ));
    std::fs::remove_dir_all(&root)?;
    Ok(())
}

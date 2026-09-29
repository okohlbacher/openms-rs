// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SpectrumAccessSqMass` over the retained upstream fixture
//! `SqliteMassFile_1.sqMass` (already registered by the sqMass handler group
//! under `tests/data/sqlite_s1_source_review/`).
//!
//! Expected values are transcribed from
//! `src/tests/class_tests/openms/source/SpectrumAccessSqMass_test.cpp` at the
//! core pin `bc9cc12` (evidence tier 3), cited per test by section. Tests
//! marked "derived" check independently derived behaviour (tier 4).

#![cfg(feature = "sqmass")]

use openms::Error;
use openms::analysis::openswath::spectrum_access::SpectrumAccess;
use openms::analysis::openswath::spectrum_access_sqmass::SpectrumAccessSqMass;
use openms::format::mzml_sqlite_handler::MzMLSqliteHandler;
use std::path::Path;
use std::sync::Arc;

fn handler() -> Arc<MzMLSqliteHandler> {
    Arc::new(MzMLSqliteHandler::new(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/sqlite_s1_source_review/SqliteMassFile_1.sqMass"),
        0,
    ))
}

/// Sections `SpectrumAccessSqMass(handler)`, `(handler, indices)`,
/// `getNrSpectra` and `lightClone`.
#[test]
fn construction_and_counts() {
    let all = SpectrumAccessSqMass::new(handler()).unwrap();
    assert_eq!(all.nr_spectra(), 2);
    assert_eq!(all.light_clone().unwrap().nr_spectra(), 2);
    let second = SpectrumAccessSqMass::with_indices(handler(), vec![1]).unwrap();
    assert_eq!(second.nr_spectra(), 1);
}

/// Section `SpectrumAccessSqMass(const SpectrumAccessSqMass& sp, indices)`.
#[test]
fn narrowing() {
    let sasm = SpectrumAccessSqMass::new(handler()).unwrap();
    assert_eq!(
        SpectrumAccessSqMass::narrowed(&sasm, &[0, 1])
            .unwrap()
            .nr_spectra(),
        2
    );
    let only_second = SpectrumAccessSqMass::narrowed(&sasm, &[1]).unwrap();
    assert_eq!(only_second.nr_spectra(), 1);
    assert!(matches!(
        SpectrumAccessSqMass::narrowed(&only_second, &[1]),
        Err(Error::InvalidValue(_))
    ));
    assert!(SpectrumAccessSqMass::narrowed(&only_second, &[50]).is_err());
    // Derived: an empty narrowing inherits the parent's view.
    let same = SpectrumAccessSqMass::narrowed(&only_second, &[]).unwrap();
    assert_eq!(same.subset(), &[1]);
}

/// Section `getAllSpectra`, all four blocks.
#[test]
fn all_spectra() {
    let check = |access: &SpectrumAccessSqMass, sizes: &[usize]| {
        assert_eq!(access.nr_spectra(), sizes.len());
        let (spectra, meta) = access.all_spectra().unwrap();
        assert_eq!(spectra.len(), sizes.len());
        assert_eq!(meta.len(), sizes.len());
        for (spectrum, &size) in spectra.iter().zip(sizes) {
            assert_eq!(spectrum.mz_array().data.len(), size);
            assert_eq!(spectrum.intensity_array().data.len(), size);
        }
    };
    check(
        &SpectrumAccessSqMass::new(handler()).unwrap(),
        &[19914, 19800],
    );
    check(
        &SpectrumAccessSqMass::with_indices(handler(), vec![0, 1]).unwrap(),
        &[19914, 19800],
    );
    let second = SpectrumAccessSqMass::with_indices(handler(), vec![1]).unwrap();
    check(&second, &[19800]);
    check(
        &SpectrumAccessSqMass::narrowed(&second, &[0]).unwrap(),
        &[19800],
    );
}

/// Derived: per-id reads agree with the batch read; the RT lookup maps SQL
/// ids back to visible positions; chromatograms are not served.
#[test]
fn per_id_reads_and_rt_lookup() {
    let all = SpectrumAccessSqMass::new(handler()).unwrap();
    let (spectra, meta) = all.all_spectra().unwrap();
    for id in 0..2 {
        assert_eq!(*all.spectrum_by_id(id).unwrap(), *spectra[id]);
        let m = all.spectrum_meta_by_id(id).unwrap();
        assert_eq!(
            (m.id.as_str(), m.rt, m.ms_level),
            (meta[id].id.as_str(), meta[id].rt, meta[id].ms_level)
        );
    }
    assert!(all.spectrum_by_id(2).is_err());
    let second = SpectrumAccessSqMass::with_indices(handler(), vec![1]).unwrap();
    assert_eq!(*second.spectrum_by_id(0).unwrap(), *spectra[1]);
    let rt1 = meta[1].rt;
    assert_eq!(all.spectra_by_rt(rt1, 0.0).unwrap(), vec![1]);
    assert_eq!(second.spectra_by_rt(rt1, 0.0).unwrap(), vec![0]);
    let wide = (meta[1].rt - meta[0].rt).abs() + 1.0;
    assert_eq!(all.spectra_by_rt(meta[0].rt, wide).unwrap(), vec![0, 1]);
    assert!(matches!(
        all.chromatogram_by_id(0),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        all.chromatogram_native_id(0),
        Err(Error::Unsupported(_))
    ));
    let memory = all.to_in_memory().unwrap();
    assert_eq!(memory.nr_spectra(), 2);
    assert_eq!(memory.nr_chromatograms(), 0);
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `NeighborSeq` against the literals of the upstream class test
//! `src/tests/class_tests/openms/source/NeighborSeq_test.cpp` at core `bc9cc12`
//! (evidence tier 3), plus native boundary checks (tier 4).

use openms::Result;
use openms::analysis::id::neighbor_seq::{NeighborSeq, NeighborStats};
use openms::chemistry::AASequence;
use openms::kernel::{MSSpectrum, Peak1D};

fn seq(text: &str) -> AASequence {
    AASequence::parse(text).expect("valid peptide")
}

fn spectrum(mzs: &[f64]) -> MSSpectrum {
    MSSpectrum {
        peaks: mzs.iter().map(|&mz| Peak1D::new(mz, 1.0)).collect(),
        ..MSSpectrum::default()
    }
}

/// `TEST_REAL_SIMILAR` default: relative 1e-5 (absolute 1e-5 near zero).
fn assert_similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} vs {expected}"
    );
}

// NeighborSeq_test.cpp:26-61
#[test]
fn generate_spectrum_matches_class_test() -> Result<()> {
    let ns = NeighborSeq::new(vec![seq("TEST")])?;
    let cases: [(&str, &[f64]); 3] = [
        (
            "PEPT",
            &[98.0600, 120.0655, 217.1182, 227.1026, 324.1553, 346.1608],
        ),
        ("AR", &[72.04439, 175.1189]),
        (
            "VGLPINQR",
            &[
                100.0756, 157.0971, 175.1189, 270.1812, 303.1775, 367.2339, 417.2204, 480.3180,
                530.3045, 594.3609, 627.3578, 722.4195, 740.4413, 797.4628,
            ],
        ),
    ];
    for (peptide, expected) in cases {
        let spec = ns.generate_spectrum(&seq(peptide))?;
        // The class test indexes these peaks; b and y at charge 1 including b1
        // give exactly 2 * (n - 1) peaks.
        assert_eq!(spec.peaks.len(), expected.len(), "{peptide}");
        for (peak, &mz) in spec.peaks.iter().zip(expected) {
            assert_similar(peak.mz, mz);
        }
    }
    Ok(())
}

fn four_spectra() -> [MSSpectrum; 4] {
    [
        spectrum(&[100.00, 200.00, 300.00]),
        spectrum(&[100.05, 200.05, 300.05]),
        spectrum(&[101.00, 201.00, 301.00]),
        spectrum(&[100.05, 201.00, 300.05, 301.00]),
    ]
}

// NeighborSeq_test.cpp:65-105
#[test]
fn is_neighbor_spectrum_matches_class_test() {
    let [s1, s2, s3, s4] = four_spectra();
    let f = 0.5;
    assert!(NeighborSeq::is_neighbor_spectrum(&s1, &s2, f, 1.0));
    assert!(!NeighborSeq::is_neighbor_spectrum(&s1, &s3, f, 1.0));
    assert!(NeighborSeq::is_neighbor_spectrum(&s1, &s4, f, 1.0));
    assert!(!NeighborSeq::is_neighbor_spectrum(&s2, &s3, f, 1.0));
    assert!(NeighborSeq::is_neighbor_spectrum(&s2, &s4, f, 1.0));
    assert!(NeighborSeq::is_neighbor_spectrum(&s3, &s4, f, 1.0));

    assert!(!NeighborSeq::is_neighbor_spectrum(&s1, &s2, f, 0.05));
    assert!(!NeighborSeq::is_neighbor_spectrum(&s1, &s3, f, 0.05));
    assert!(!NeighborSeq::is_neighbor_spectrum(&s1, &s4, f, 0.05));
    assert!(!NeighborSeq::is_neighbor_spectrum(&s2, &s3, f, 0.05));
    assert!(NeighborSeq::is_neighbor_spectrum(&s2, &s4, f, 0.05));
    assert!(NeighborSeq::is_neighbor_spectrum(&s3, &s4, f, 0.05));
}

// NeighborSeq_test.cpp:109-150
#[test]
fn compute_shared_ion_count_matches_class_test() {
    let [s1, s2, s3, s4] = four_spectra();
    let count = NeighborSeq::compute_shared_ion_count;
    assert_eq!(count(&s1, &s2, 2.0), 3);
    assert_eq!(count(&s1, &s3, 2.0), 3);
    assert_eq!(count(&s1, &s4, 2.0), 3);
    assert_eq!(count(&s2, &s3, 2.0), 3);
    assert_eq!(count(&s2, &s4, 2.0), 3);
    assert_eq!(count(&s3, &s4, 2.0), 3);
    for bin in [1.0, 0.1] {
        assert_eq!(count(&s1, &s2, bin), 3);
        assert_eq!(count(&s1, &s3, bin), 0);
        assert_eq!(count(&s1, &s4, bin), 2);
        assert_eq!(count(&s2, &s3, bin), 0);
        assert_eq!(count(&s2, &s4, bin), 2);
        assert_eq!(count(&s3, &s4, bin), 2);
    }
}

// NeighborSeq_test.cpp:152-199
#[test]
fn is_neighbor_peptide_matches_class_test() -> Result<()> {
    let velqsk = seq("VELQSK");
    let svqelk = seq("SVQELK");
    let tvdqlk = seq("TVDQLK");
    let vesqlk = seq("VESQLK");
    let mut ns = NeighborSeq::new(vec![
        seq("VELQSK"),
        seq("SVQELK"),
        seq("TVDQLK"),
        seq("VGEFK"),
    ])?;
    let pc_tolerance = 0.01;
    let mz_bin_size = 0.05;
    assert!((velqsk.mono_mass()? - vesqlk.mono_mass()?).abs() < pc_tolerance);
    let shared = |a: &AASequence, b: &AASequence| -> Result<usize> {
        Ok(NeighborSeq::compute_shared_ion_count(
            &ns.generate_spectrum(a)?,
            &ns.generate_spectrum(b)?,
            mz_bin_size,
        ))
    };
    assert_eq!(shared(&velqsk, &vesqlk)?, 6);
    assert_eq!(shared(&svqelk, &vesqlk)?, 4);
    assert_eq!(shared(&tvdqlk, &vesqlk)?, 6);

    let shared_ion_fraction = 6.0 * 2.0 / (((vesqlk.len() - 1) * 2 * 2) as f64);
    assert!(!ns.is_neighbor_peptide(
        &vesqlk,
        pc_tolerance,
        false,
        shared_ion_fraction + 0.1,
        mz_bin_size
    )?);
    assert!(ns.is_neighbor_peptide(
        &vesqlk,
        pc_tolerance,
        false,
        shared_ion_fraction - 0.1,
        mz_bin_size
    )?);
    assert!(ns.is_neighbor_peptide(&seq("GLDFK"), pc_tolerance, false, 0.25, mz_bin_size)?);

    let stats = ns.neighbor_stats();
    assert_eq!(stats.unfindable_peptides, 0);
    assert_eq!(stats.findable_no_neighbors, 1);
    assert_eq!(stats.findable_one_neighbor, 3);
    assert_eq!(stats.findable_multiple_neighbors, 0);

    assert!(ns.is_neighbor_peptide(&vesqlk, pc_tolerance, false, 0.25, mz_bin_size)?);
    let stats2 = ns.neighbor_stats();
    assert_eq!(stats2.unfindable_peptides, 0);
    assert_eq!(stats2.findable_no_neighbors, 0);
    assert_eq!(stats2.findable_one_neighbor, 2);
    assert_eq!(stats2.findable_multiple_neighbors, 2);
    Ok(())
}

// NeighborSeq_test.cpp:201-218 ([EXTRA] owns its moved-in peptides)
#[test]
fn owns_its_peptides() -> Result<()> {
    let mut source = vec![seq("VELQSK"), seq("TVDQLK")];
    let mut ns = NeighborSeq::new(std::mem::take(&mut source))?;
    source = vec![seq(&"A".repeat(100)), seq(&"W".repeat(100))];
    assert_eq!(source.len(), 2);
    assert!(ns.is_neighbor_peptide(&seq("VESQLK"), 0.01, false, 0.25, 0.05)?);
    assert_eq!(ns.neighbor_stats().total(), 2);
    Ok(())
}

// NeighborSeq_test.cpp:226-245 ([EXTRA] percentage formatters)
#[test]
fn neighbor_stats_formatters_match_class_test() {
    let s = NeighborStats::default();
    assert_eq!(s.total(), 0);
    assert_eq!(s.unfindable(), "0 (0%)");
    assert_eq!(s.no_nb(), "0 (0%)");
    assert_eq!(s.one_nb(), "0 (0%)");
    assert_eq!(s.multi_nb(), "0 (0%)");

    let s2 = NeighborStats {
        findable_no_neighbors: 1,
        findable_one_neighbor: 3,
        ..NeighborStats::default()
    };
    assert_eq!(s2.total(), 4);
    assert_eq!(s2.no_nb(), "1 (25%)");
    assert_eq!(s2.one_nb(), "3 (75%)");
}

// Tier 4: `X` peptides are unfindable and never matched; the source marks them -1.
#[test]
fn unknown_residue_peptides_are_unfindable() -> Result<()> {
    let mut ns = NeighborSeq::new(vec![seq("VELXSK"), seq("VELQSK")])?;
    assert_eq!(ns.relevant_peptides().len(), 2);
    assert!(ns.is_neighbor_peptide(&seq("VESQLK"), 0.01, false, 0.25, 0.05)?);
    let stats = ns.neighbor_stats();
    assert_eq!(stats.unfindable_peptides, 1);
    assert_eq!(stats.findable_one_neighbor, 1);
    assert_eq!(stats.unfindable(), "1 (50%)");
    Ok(())
}

// Tier 4: a ppm window equals the Da window it converts to; the window is
// inclusive at both ends like lower_bound/upper_bound. VELKSK differs from
// VELQSK by Q->K (0.036 Da) and shares b1-b3 and y1-y2 with it.
#[test]
fn ppm_window_and_boundaries() -> Result<()> {
    let candidate = seq("VELKSK");
    let relevant = seq("VELQSK");
    let delta = (candidate.mono_mass()? - relevant.mono_mass()?).abs();
    assert!(delta > 0.03);
    let ppm_needed = delta / candidate.mono_mass()? * 1e6;
    let mut ns = NeighborSeq::new(vec![relevant])?;
    assert!(!ns.is_neighbor_peptide(&candidate, ppm_needed * 0.99, true, 0.25, 0.05)?);
    assert!(ns.is_neighbor_peptide(&candidate, ppm_needed * 1.01, true, 0.25, 0.05)?);
    assert!(!ns.is_neighbor_peptide(&candidate, delta * 0.99, false, 0.25, 0.05)?);
    assert!(ns.is_neighbor_peptide(&candidate, delta * 1.01, false, 0.25, 0.05)?);
    // Identical mass: a zero tolerance still matches (closed window).
    assert!(ns.is_neighbor_peptide(&seq("VELQSK"), 0.0, false, 0.25, 0.05)?);
    assert_eq!(ns.neighbor_stats().findable_multiple_neighbors, 1);
    Ok(())
}

// Tier 4: invalid tolerances are refused instead of the source's debug assert,
// and a refused query leaves the counters unchanged.
#[test]
fn invalid_tolerance_is_refused() -> Result<()> {
    let mut ns = NeighborSeq::new(vec![seq("VELQSK")])?;
    for tolerance in [-0.01, f64::NAN, f64::INFINITY] {
        assert!(
            ns.is_neighbor_peptide(&seq("VESQLK"), tolerance, false, 0.25, 0.05)
                .is_err()
        );
    }
    assert_eq!(ns.neighbor_stats().findable_no_neighbors, 1);
    Ok(())
}

// Tier 4: empty spectra share nothing and are never neighbors (0/0 is NaN).
#[test]
fn empty_spectra_are_not_neighbors() {
    let empty = MSSpectrum::default();
    assert_eq!(
        NeighborSeq::compute_shared_ion_count(&empty, &empty, 0.05),
        0
    );
    assert!(!NeighborSeq::is_neighbor_spectrum(
        &empty, &empty, 0.0, 0.05
    ));
}

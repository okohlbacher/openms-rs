// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! ANALYSIS/OPENSWATH/DATAACCESS, `ISpectrumAccess`, `OpenSwathHelper` and
//! `SpectrumAddition`.
//!
//! Expected values are transcribed from the upstream class tests at the core
//! pin `bc9cc12` (evidence tier 3), cited per test as `<file>:<section>`, all
//! under `src/tests/class_tests/openms/source/`. Tests marked "derived" check
//! an independently derived invariant or a native guard (tier 4).

use openms::Error;
use openms::analysis::openswath::data_access_helper as dah;
use openms::analysis::openswath::data_structures::{BinaryDataArray, Chromatogram, Spectrum};
use openms::analysis::openswath::itransition::{
    FeatureAccess, MRMFeatureAccess, SignalToNoiseAccess, TransitionGroupAccess,
};
use openms::analysis::openswath::mrm_feature_access::{
    FeatureOpenMS, MRMFeatureOpenMS, SignalToNoiseOpenMS, TransitionGroupOpenMS,
};
use openms::analysis::openswath::openswath_helper as helper;
use openms::analysis::openswath::spectrum_access::{
    SpectrumAccess, SpectrumAccessPtr, filter_by_drift,
};
use openms::analysis::openswath::spectrum_access_openms::{
    QuadMZTransform, SimpleOpenMSSpectraFactory, SpectrumAccessOpenMS,
    SpectrumAccessOpenMSInMemory, SpectrumAccessQuadMZTransforming, SpectrumAccessTransforming,
};
use openms::analysis::openswath::spectrum_addition as addition;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::{
    LightCompound, LightProtein, LightTargetedExperiment, LightTransition,
};
use openms::kernel::features::Feature;
use openms::kernel::mrm::{MRMFeature, MRMTransitionGroup, SimpleTransition};
use openms::kernel::ranges::RangeBase;
use openms::kernel::{
    ChromatogramPeak, ConvexHull2D, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
    Point2D, Precursor,
};
use openms::metadata::{DataProcessing, MetaValue};
use std::collections::BTreeMap;
use std::sync::Arc;

fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

/// The two-spectrum, one-chromatogram experiment of SpectrumAccessOpenMS_test.cpp
/// and SpectrumAccessOpenMSInMemory_test.cpp (`mkExp`).
fn mk_exp() -> Arc<MSExperiment> {
    let mut exp = MSExperiment::new();
    exp.spectra.push(MSSpectrum {
        rt: 10.0,
        ms_level: 1,
        peaks: vec![Peak1D::new(100.0, 500.0), Peak1D::new(200.0, 800.0)],
        ..MSSpectrum::default()
    });
    exp.spectra.push(MSSpectrum {
        rt: 20.0,
        ms_level: 2,
        ..MSSpectrum::default()
    });
    exp.chromatograms.push(MSChromatogram {
        native_id: "chrom_42".into(),
        peaks: vec![ChromatogramPeak::new(5.0, 99.0)],
        ..MSChromatogram::default()
    });
    Arc::new(exp)
}

// ---------------------------------------------------------------------------
// SpectrumAccessOpenMS
// ---------------------------------------------------------------------------

/// SpectrumAccessOpenMS_test.cpp: every section.
#[test]
fn spectrum_access_openms() {
    let sa = SpectrumAccessOpenMS::new(mk_exp());
    assert_eq!(sa.nr_spectra(), 2);
    assert_eq!(sa.nr_chromatograms(), 1);
    let s = sa.spectrum_by_id(0).unwrap();
    assert_eq!(s.mz_array().data.len(), 2);
    real_similar(s.mz_array().data[0], 100.0);
    real_similar(s.mz_array().data[1], 200.0);
    real_similar(s.intensity_array().data[0], 500.0);
    real_similar(s.intensity_array().data[1], 800.0);
    let m = sa.spectrum_meta_by_id(0).unwrap();
    real_similar(m.rt, 10.0);
    assert_eq!(m.ms_level, 1);
    assert_eq!(sa.spectra_by_rt(10.0, 0.0).unwrap(), vec![0]);
    assert_eq!(sa.spectra_by_rt(15.0, 6.0).unwrap().len(), 2);
    let c = sa.chromatogram_by_id(0).unwrap();
    assert_eq!(c.time_array().data.len(), 1);
    real_similar(c.time_array().data[0], 5.0);
    real_similar(c.intensity_array().data[0], 99.0);
    assert_eq!(sa.chromatogram_native_id(0).unwrap(), "chrom_42");
    let clone = sa.light_clone().unwrap();
    assert_eq!(clone.nr_spectra(), 2);
    assert_eq!(clone.nr_chromatograms(), 1);
}

/// Derived from SpectrumAccessOpenMS.cpp: the RT window keeps the first
/// spectrum at or after `rt - delta` even when it lies beyond `rt + delta`,
/// includes the upper bound, and the view carries data arrays; out-of-range ids
/// and unsorted experiments are refused.
#[test]
fn spectrum_access_openms_rules() {
    let sa = SpectrumAccessOpenMS::new(mk_exp());
    assert_eq!(sa.spectra_by_rt(12.0, 0.0).unwrap(), vec![1]);
    assert_eq!(sa.spectra_by_rt(15.0, 5.0).unwrap(), vec![0, 1]);
    assert!(sa.spectra_by_rt(25.0, 1.0).unwrap().is_empty());
    assert!(matches!(
        sa.spectra_by_rt(1.0, -1.0),
        Err(Error::InvalidValue(_))
    ));
    assert!(sa.spectrum_by_id(2).is_err());
    assert!(sa.chromatogram_by_id(1).is_err());
    let meta = sa.spectrum_meta_by_id(1).unwrap();
    assert_eq!(
        (meta.index, meta.id.as_str()),
        (0, ""),
        "only RT and MS level are filled"
    );
    assert_eq!(sa.spectra_meta_info(1).unwrap().ms_level, 2);
    assert_eq!(sa.chromatogram_meta_info(0).unwrap().native_id, "chrom_42");

    let mut exp = (*mk_exp()).clone();
    exp.spectra[0]
        .float_data_arrays
        .push(DataArray::new("Ion Mobility", vec![1.5f32, 2.5]));
    exp.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("charge", vec![1, 2]));
    exp.spectra[0].string_data_arrays.push(DataArray::new(
        "names",
        vec!["a".to_string(), "b".to_string()],
    ));
    let sa = SpectrumAccessOpenMS::new(Arc::new(exp.clone()));
    let s = sa.spectrum_by_id(0).unwrap();
    assert_eq!(
        s.data_arrays().len(),
        4,
        "float and integer arrays, not strings"
    );
    assert_eq!(s.drift_time_array().unwrap().data, vec![1.5, 2.5]);
    assert_eq!(s.data_arrays()[3].description, "charge");
    exp.spectra.swap(0, 1);
    let unsorted = SpectrumAccessOpenMS::new(Arc::new(exp));
    assert!(matches!(
        unsorted.spectra_by_rt(10.0, 1.0),
        Err(Error::UnsortedData)
    ));
    assert!(unsorted.spectrum_by_id(0).is_ok());
}

// ---------------------------------------------------------------------------
// SpectrumAccessOpenMSInMemory
// ---------------------------------------------------------------------------

/// SpectrumAccessOpenMSInMemory_test.cpp: every section.
#[test]
fn spectrum_access_in_memory() {
    let origin = SpectrumAccessOpenMS::new(mk_exp());
    let sa = SpectrumAccessOpenMSInMemory::new(&origin).unwrap();
    assert_eq!(sa.nr_spectra(), 2);
    assert_eq!(sa.nr_chromatograms(), 1);
    let s = sa.spectrum_by_id(0).unwrap();
    assert_eq!(s.mz_array().data.len(), 2);
    real_similar(s.mz_array().data[0], 100.0);
    real_similar(s.mz_array().data[1], 200.0);
    real_similar(s.intensity_array().data[1], 800.0);
    real_similar(sa.spectrum_meta_by_id(0).unwrap().rt, 10.0);
    let c = sa.chromatogram_by_id(0).unwrap();
    assert_eq!(c.time_array().data.len(), 1);
    real_similar(c.time_array().data[0], 5.0);
    real_similar(c.intensity_array().data[0], 99.0);
    assert_eq!(sa.light_clone().unwrap().nr_spectra(), 2);
}

/// Derived from SpectrumAccessOpenMSInMemory.cpp: shared payloads, exclusive
/// upper bound, and the end of the list is checked before it is read.
#[test]
fn spectrum_access_in_memory_rules() {
    let origin = SpectrumAccessOpenMS::new(mk_exp());
    let sa = SpectrumAccessOpenMSInMemory::new(&origin).unwrap();
    let a = sa.spectrum_by_id(1).unwrap();
    let b = sa.spectrum_by_id(1).unwrap();
    assert!(Arc::ptr_eq(&a, &b), "the cached handle is shared");
    assert_eq!(
        sa.spectra_by_rt(15.0, 5.0).unwrap(),
        vec![0],
        "upper bound exclusive"
    );
    assert_eq!(
        sa.spectra_by_rt(20.0, 0.0).unwrap(),
        vec![1],
        "last spectrum, no overrun"
    );
    assert_eq!(sa.chromatogram_native_id(0).unwrap(), "chrom_42");
    let spectra_only = SpectrumAccessOpenMSInMemory::from_spectra(&origin).unwrap();
    assert_eq!(spectra_only.nr_chromatograms(), 0);
    assert!(
        SpectrumAccessOpenMSInMemory::from_parts(vec![a], Vec::new(), Vec::new(), Vec::new())
            .is_err()
    );
}

/// Derived from ISpectrumAccess.cpp: `getMultipleSpectra` picks the closest
/// spectrum and alternates outward, lower neighbour first, `n / 2` steps.
#[test]
fn multiple_spectra_order() {
    let mut exp = MSExperiment::new();
    for i in 0..6 {
        exp.spectra.push(MSSpectrum {
            rt: 10.0 * i as f64,
            peaks: vec![Peak1D::new(100.0 + i as f64, 1.0)],
            ..MSSpectrum::default()
        });
    }
    let sa = SpectrumAccessOpenMS::new(Arc::new(exp));
    let first_mz = |seq: Vec<Arc<Spectrum>>| -> Vec<f64> {
        seq.iter().map(|s| s.mz_array().data[0]).collect()
    };
    // RT 24: first at or after is index 3 (RT 30); index 2 (RT 20) is closer.
    assert_eq!(
        first_mz(sa.multiple_spectra(24.0, 3).unwrap()),
        vec![102.0, 101.0, 103.0]
    );
    // An even count fetches one more (4 / 2 = 2 steps).
    assert_eq!(
        first_mz(sa.multiple_spectra(24.0, 4).unwrap()),
        vec![102.0, 101.0, 103.0, 100.0, 104.0]
    );
    // Neighbours outside the view are skipped.
    assert_eq!(
        first_mz(sa.multiple_spectra(0.0, 3).unwrap()),
        vec![100.0, 101.0]
    );
    assert_eq!(first_mz(sa.multiple_spectra(0.0, 0).unwrap()), vec![100.0]);
    assert!(sa.multiple_spectra(1000.0, 3).unwrap().is_empty());
    // Equidistant: the later spectrum wins (the earlier must be strictly closer).
    assert_eq!(first_mz(sa.multiple_spectra(25.0, 1).unwrap()), vec![103.0]);
}

/// Derived from ISpectrumAccess.h: `filterByDrift` keeps peaks with the drift
/// time inside the inclusive window, in order, and refuses a spectrum without
/// an ion-mobility array.
#[test]
fn filter_by_drift_rules() {
    let mut s = Spectrum::from_arrays(vec![1.0, 2.0, 3.0, 4.0], vec![10.0, 20.0, 30.0, 40.0]);
    s.push_data_array(BinaryDataArray::with_description(vec![0.0], "other"));
    assert!(matches!(
        filter_by_drift(&s, 0.0, 1.0),
        Err(Error::MissingInformation(_))
    ));
    let mut s = Spectrum::from_arrays(vec![1.0, 2.0, 3.0, 4.0], vec![10.0, 20.0, 30.0, 40.0]);
    s.push_data_array(BinaryDataArray::with_description(
        vec![0.5, 0.9, 0.7, 1.1],
        "mean inverse reduced ion mobility array",
    ));
    let f = filter_by_drift(&s, 0.7, 0.9).unwrap();
    assert_eq!(f.mz_array().data, vec![2.0, 3.0]);
    assert_eq!(f.intensity_array().data, vec![20.0, 30.0]);
    let im = f.drift_time_array().unwrap();
    assert_eq!(im.data, vec![0.9, 0.7]);
    assert_eq!(im.description, "mean inverse reduced ion mobility array");
    let mut short = Spectrum::from_arrays(vec![1.0, 2.0], vec![1.0, 2.0]);
    short.set_drift_time_array(BinaryDataArray::new(vec![0.5]));
    assert!(matches!(
        filter_by_drift(&short, 0.0, 1.0),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// SpectrumAccessQuadMZTransforming / SimpleOpenMSSpectraFactory
// ---------------------------------------------------------------------------

/// SpectrumAccessQuadMZTransforming_test.cpp `getData`.
fn quad_data() -> Arc<MSExperiment> {
    let mut exp = MSExperiment::new();
    exp.spectra.push(MSSpectrum {
        peaks: vec![Peak1D::new(100.0, 50.0), Peak1D::new(500.0, 150.0)],
        ..MSSpectrum::default()
    });
    Arc::new(exp)
}

/// SpectrumAccessQuadMZTransforming_test.cpp: every section.
#[test]
fn quad_mz_transforming() {
    let empty =
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(MSExperiment::new()))
            .unwrap();
    let t = SpectrumAccessQuadMZTransforming::quad_mz(empty, 0.0, 0.0, 0.0, false);
    assert_eq!(t.nr_spectra(), 0);
    assert_eq!(t.light_clone().unwrap().nr_spectra(), 0);

    let expptr2 = SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(quad_data()).unwrap();
    let t = SpectrumAccessQuadMZTransforming::quad_mz(expptr2.clone(), 0.0, 0.0, 0.0, false);
    assert_eq!(t.nr_spectra(), 1);
    let s = t.spectrum_by_id(0).unwrap();
    assert_eq!(s.mz_array().data.len(), 2);
    real_similar(s.intensity_array().data[0], 50.0);
    real_similar(s.intensity_array().data[1], 150.0);
    real_similar(s.mz_array().data[0], 0.0);
    real_similar(s.mz_array().data[1], 0.0);

    let t = SpectrumAccessQuadMZTransforming::quad_mz(expptr2, 10.0, 5.0, 2.0, false);
    let check = |a: &dyn SpectrumAccess| {
        let s = a.spectrum_by_id(0).unwrap();
        real_similar(s.intensity_array().data[0], 50.0);
        real_similar(s.intensity_array().data[1], 150.0);
        real_similar(
            s.mz_array().data[0],
            10.0 + 100.0 * 5.0 + 100.0 * 100.0 * 2.0,
        );
        real_similar(
            s.mz_array().data[1],
            10.0 + 500.0 * 5.0 + 500.0 * 500.0 * 2.0,
        );
    };
    check(&t);
    let clone = t.light_clone().unwrap();
    assert_eq!(clone.nr_spectra(), t.nr_spectra());
    check(&t);
    check(clone.as_ref());
}

/// Derived: over a materialised view, repeated requests transform once (the
/// source transforms the shared cached spectrum in place, so its second
/// request would return a doubly transformed m/z).
#[test]
fn quad_mz_does_not_corrupt_a_shared_cache() {
    let origin = SpectrumAccessOpenMS::new(quad_data());
    let cache: SpectrumAccessPtr = Arc::new(SpectrumAccessOpenMSInMemory::new(&origin).unwrap());
    let t = SpectrumAccessTransforming::new(
        cache.clone(),
        QuadMZTransform {
            a: 1.0,
            b: 1.0,
            c: 0.0,
            ppm: false,
        },
    );
    for _ in 0..3 {
        assert_eq!(
            t.spectrum_by_id(0).unwrap().mz_array().data,
            vec![101.0, 501.0]
        );
    }
    assert_eq!(
        cache.spectrum_by_id(0).unwrap().mz_array().data,
        vec![100.0, 500.0]
    );
    let ppm = QuadMZTransform {
        a: 10.0,
        b: 0.0,
        c: 0.0,
        ppm: true,
    };
    real_similar(ppm.apply(500.0), 500.0 - 10.0 * 500.0 / 1e6);
}

/// Derived from SimpleOpenMSSpectraAccessFactory.cpp: a `cached_data`
/// data-processing entry selects the (unported) cached view.
#[test]
fn factory_detects_cached_experiments() {
    let mut exp = (*mk_exp()).clone();
    assert!(!SimpleOpenMSSpectraFactory::is_experiment_cached(&exp));
    let mut dp = DataProcessing::default();
    dp.metadata
        .insert("cached_data".into(), MetaValue::from("true"));
    exp.chromatograms[0].data_processing.push(Arc::new(dp));
    assert!(SimpleOpenMSSpectraFactory::is_experiment_cached(&exp));
    assert!(matches!(
        SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(exp)),
        Err(Error::Unsupported(_))
    ));
}

// ---------------------------------------------------------------------------
// DataAccessHelper
// ---------------------------------------------------------------------------

/// Derived from DataAccessHelper.cpp: the four conversions and the RT filter.
#[test]
fn data_access_helper_conversions() {
    let spectrum = MSSpectrum {
        rt: 7.0,
        native_id: "keep".into(),
        peaks: vec![Peak1D::new(100.0, 5.0), Peak1D::new(101.0, 6.0)],
        float_data_arrays: vec![DataArray::new("x", vec![1.0f32, 2.0])],
        ..MSSpectrum::default()
    };
    let light = dah::convert_to_spectrum(&spectrum);
    assert_eq!(
        light.data_arrays().len(),
        2,
        "data arrays are not converted"
    );
    assert_eq!(light.intensity_array().data, vec![5.0, 6.0]);
    let mut back = spectrum.clone();
    dah::convert_to_openms_spectrum(&Spectrum::from_arrays(vec![1.0], vec![2.5]), &mut back)
        .unwrap();
    assert_eq!(back.peaks, vec![Peak1D::new(1.0, 2.5)]);
    assert_eq!((back.rt, back.native_id.as_str()), (7.0, "keep"));
    assert!(back.float_data_arrays.is_empty());
    let bad = Spectrum::from_arrays(vec![1.0, 2.0], vec![1.0]);
    assert!(dah::convert_to_openms_spectrum(&bad, &mut back).is_err());
    assert_eq!(back.peaks.len(), 1, "unchanged on error");

    let chrom = MSChromatogram {
        peaks: vec![
            ChromatogramPeak::new(1.0, 10.0),
            ChromatogramPeak::new(2.0, 20.0),
            ChromatogramPeak::new(3.0, 30.0),
        ],
        ..MSChromatogram::default()
    };
    let light = dah::convert_to_chromatogram(&chrom);
    assert_eq!(light.time_array().data, vec![1.0, 2.0, 3.0]);
    let mut out = MSChromatogram::default();
    dah::convert_to_openms_chromatogram(&light, &mut out).unwrap();
    assert_eq!(out.peaks, chrom.peaks);
    let kept = dah::convert_to_openms_chromatogram_filter(&mut out, &light, 1.5, 3.0).unwrap();
    assert_eq!(kept, 2);
    assert_eq!(out.peaks[0].rt, 2.0);
    let bad = Chromatogram::from_arrays(vec![1.0], vec![]);
    assert!(dah::convert_to_openms_chromatogram(&bad, &mut out).is_err());
}

// ---------------------------------------------------------------------------
// MRMFeatureAccessOpenMS
// ---------------------------------------------------------------------------

/// MRMFeatureAccessOpenMS_test.cpp `mkFeat`.
fn mk_feat(rt: f64, intensity: f32, with_hull: bool) -> Feature {
    let mut f = Feature::new(rt, 0.0, intensity);
    if with_hull {
        let mut hull = ConvexHull2D::new();
        hull.set_hull_points(&[
            Point2D { rt: 95.0, mz: 10.0 },
            Point2D {
                rt: 105.0,
                mz: 20.0,
            },
        ])
        .unwrap();
        f.convex_hulls = vec![hull];
    }
    f
}

/// MRMFeatureAccessOpenMS_test.cpp: the FeatureOpenMS sections.
#[test]
fn feature_openms() {
    let f = mk_feat(100.0, 500.0, true);
    let fo = FeatureOpenMS::new(&f);
    real_similar(fo.rt(), 100.0);
    real_similar(f64::from(fo.intensity()), 500.0);
    let rt = fo.rt_values().unwrap();
    assert_eq!(rt.len(), 2);
    real_similar(rt[0], 95.0);
    real_similar(rt[1], 105.0);
    let intens = fo.intensity_values().unwrap();
    real_similar(intens[0], 10.0);
    real_similar(intens[1], 20.0);
    let no_hull = mk_feat(50.0, 5.0, false);
    let fo = FeatureOpenMS::new(&no_hull);
    assert!(matches!(fo.rt_values(), Err(Error::InvalidValue(_))));
    assert!(matches!(fo.intensity_values(), Err(Error::InvalidValue(_))));
}

/// MRMFeatureAccessOpenMS_test.cpp: the MRMFeatureOpenMS sections.
#[test]
fn mrm_feature_openms() {
    let mut mrm = MRMFeature::new();
    mrm.rt = 200.0;
    mrm.intensity = 1000.0;
    mrm.metadata
        .insert("score".into(), MetaValue::try_from(0.75).unwrap());
    mrm.add_feature(mk_feat(198.0, 300.0, true), "tr1").unwrap();
    mrm.add_feature(mk_feat(202.0, 700.0, true), "tr2").unwrap();
    let mo = MRMFeatureOpenMS::new(&mrm);
    real_similar(mo.rt(), 200.0);
    real_similar(f64::from(mo.intensity()), 1000.0);
    assert_eq!(mo.len(), 2);
    real_similar(mo.meta_value("score").unwrap(), 0.75);
    let ids = mo.native_ids();
    assert_eq!(ids.len(), 2);
    assert!(ids.contains(&"tr1".to_string()) && ids.contains(&"tr2".to_string()));
    let f1 = mo.feature("tr1").unwrap();
    real_similar(f1.rt(), 198.0);
    real_similar(f64::from(f1.intensity()), 300.0);
    real_similar(f64::from(mo.feature("tr2").unwrap().intensity()), 700.0);
    // Derived: unknown IDs and missing meta values.
    assert!(mo.feature("tr3").is_none());
    assert!(matches!(
        mo.meta_value("nope"),
        Err(Error::MissingInformation(_))
    ));
    assert!(mo.precursor_ids().is_empty());
}

/// Derived from MRMFeatureAccessOpenMS.cpp: explicit and aligned lanes.
#[test]
fn mrm_feature_openms_lanes() {
    let mut mrm = MRMFeature::new();
    mrm.add_feature(mk_feat(198.0, 300.0, true), "tr1").unwrap();
    mrm.add_feature(mk_feat(202.0, 700.0, true), "tr2").unwrap();
    mrm.add_precursor_feature(mk_feat(200.0, 50.0, true), "p0")
        .unwrap();
    let ids = |v: &[&str]| -> Vec<String> { v.iter().map(|s| s.to_string()).collect() };
    let explicit = MRMFeatureOpenMS::with_ids(&mrm, &ids(&["tr2", "tr1"]), &ids(&["p0"])).unwrap();
    assert_eq!(explicit.native_ids(), ["tr2", "tr1"]);
    assert_eq!(
        explicit.feature_intensity("tr1", 0).unwrap(),
        300.0,
        "hint not trusted"
    );
    assert_eq!(
        explicit.feature_intensities(&ids(&["tr1"])).unwrap(),
        vec![vec![10.0, 20.0]]
    );
    assert!(
        explicit
            .precursor_feature_intensities(&ids(&["missing"]))
            .is_err()
    );
    assert!(MRMFeatureOpenMS::with_ids(&mrm, &ids(&["nope"]), &[]).is_err());

    let aligned = MRMFeatureOpenMS::with_aligned_ids(
        &mrm,
        &ids(&["A", "B"]),
        &ids(&["P"]),
        &ids(&["tr1", "tr2"]),
        &ids(&["p0"]),
    )
    .unwrap();
    assert_eq!(aligned.native_ids(), ["A", "B"]);
    assert_eq!(aligned.precursor_ids(), ["P"]);
    assert_eq!(aligned.feature_intensity("B", 1).unwrap(), 700.0);
    assert_eq!(aligned.feature_intensity("B", 0).unwrap(), 700.0);
    assert_eq!(aligned.precursor_feature("P").unwrap().intensity(), 50.0);
    assert_eq!(
        aligned
            .feature_intensities(&ids(&["A", "B"]))
            .unwrap()
            .len(),
        2
    );
    assert!(
        MRMFeatureOpenMS::with_aligned_ids(&mrm, &ids(&["A"]), &[], &ids(&["tr1"]), &[]).is_err()
    );
    assert!(
        MRMFeatureOpenMS::with_aligned_ids(
            &mrm,
            &ids(&["A", "B"]),
            &ids(&["P"]),
            &ids(&["x", "y"]),
            &[]
        )
        .is_err()
    );
}

/// MRMFeatureAccessOpenMS_test.cpp: the TransitionGroupOpenMS sections.
#[test]
fn transition_group_openms() {
    let mut tg: MRMTransitionGroup<MSChromatogram, SimpleTransition> = MRMTransitionGroup::new();
    tg.add_transition(SimpleTransition::new("tr1", 0.8), "tr1")
        .unwrap();
    tg.add_transition(SimpleTransition::new("tr2", 0.3), "tr2")
        .unwrap();
    for id in ["tr1", "tr2"] {
        tg.add_chromatogram(
            MSChromatogram {
                native_id: id.into(),
                ..MSChromatogram::default()
            },
            id,
        )
        .unwrap();
    }
    let tgo = TransitionGroupOpenMS::new(&tg);
    assert_eq!(tgo.len(), 2);
    assert_eq!(tgo.native_ids(), ["tr1", "tr2"]);
    let libs = tgo.library_intensities().unwrap();
    real_similar(libs[0], 0.8);
    real_similar(libs[1], 0.3);
}

/// MRMFeatureAccessOpenMS_test.cpp: the SignalToNoiseOpenMS section.
#[test]
fn signal_to_noise_openms() {
    let mut chrom = MSChromatogram::default();
    for i in 0..30 {
        chrom.peaks.push(ChromatogramPeak::new(i as f64, 10.0));
    }
    chrom.peaks[15].intensity = 1000.0;
    let sn = SignalToNoiseOpenMS::from_chromatogram(&chrom, 30.0, 30, false).unwrap();
    let apex = sn.value_at_rt(15.0);
    let base = sn.value_at_rt(2.0);
    assert!(apex > base);
    real_similar(apex / base, 102.352941176471);
    let empty = SignalToNoiseOpenMS::from_chromatogram(&MSChromatogram::default(), 30.0, 30, false)
        .unwrap();
    real_similar(empty.value_at_rt(5.0), -1.0);
}

// ---------------------------------------------------------------------------
// OpenSwathHelper
// ---------------------------------------------------------------------------

/// OpenSwathHelper_test.cpp: computePrecursorId and computeTransitionGroupId.
#[test]
fn precursor_ids() {
    assert_eq!(
        helper::compute_precursor_id("tr_gr2", 0),
        "tr_gr2_Precursor_i0"
    );
    assert_eq!(
        helper::compute_precursor_id("tr_gr2__test", 0),
        "tr_gr2__test_Precursor_i0"
    );
    assert_eq!(
        helper::compute_transition_group_id("tr_gr2_Precursor_i0"),
        "tr_gr2"
    );
    assert_eq!(
        helper::compute_transition_group_id("tr_gr2__test_Precursor_i0"),
        "tr_gr2__test"
    );
    // Derived: three parts, fewer parts, empty.
    assert_eq!(helper::compute_transition_group_id("g_Precursor_i1"), "g");
    assert_eq!(helper::compute_transition_group_id("g_i1"), "");
    assert_eq!(helper::compute_transition_group_id(""), "");
}

/// OpenSwathHelper_test.cpp: selectSwathTransitions (LightTargetedExperiment).
#[test]
fn select_swath_transitions() {
    let mut exp1 = LightTargetedExperiment::new();
    for mz in [100.0, 200.0, 300.0] {
        exp1.transitions.push(LightTransition {
            precursor_mz: mz,
            ..LightTransition::default()
        });
    }
    let mut exp2 = LightTargetedExperiment::new();
    helper::select_swath_transitions(&exp1, &mut exp2, 1.0, 199.9, 500.0);
    assert_eq!(exp2.transitions.len(), 2);
}

/// Derived from OpenSwathHelper.cpp: compounds and proteins follow the selected
/// transitions; the upper-edge distance applies.
#[test]
fn select_swath_transitions_follows_references() {
    let mut exp = LightTargetedExperiment::new();
    for (name, pep, mz) in [
        ("t1", "p1", 450.0),
        ("t2", "p2", 520.0),
        ("t3", "p3", 499.5),
    ] {
        exp.transitions.push(LightTransition {
            transition_name: name.into(),
            peptide_ref: pep.into(),
            precursor_mz: mz,
            ..LightTransition::default()
        });
    }
    for (id, prot) in [("p1", "A"), ("p2", "B"), ("p3", "A")] {
        exp.compounds.push(LightCompound {
            id: id.into(),
            protein_refs: vec![prot.into()],
            ..LightCompound::default()
        });
    }
    for id in ["A", "B"] {
        exp.proteins.push(LightProtein {
            id: id.into(),
            ..LightProtein::default()
        });
    }
    let mut out = LightTargetedExperiment::new();
    helper::select_swath_transitions(&exp, &mut out, 1.0, 400.0, 500.0);
    assert_eq!(
        out.transitions.len(),
        1,
        "t3 is within 1 Th of the upper edge"
    );
    assert_eq!(out.compounds.len(), 1);
    assert_eq!(out.proteins.len(), 1);
    assert_eq!(out.proteins[0].id, "A");
}

/// OpenSwathHelper_test.cpp: checkSwathMap.
#[test]
fn check_swath_map() {
    let mut map = MSExperiment::new();
    map.spectra.push(MSSpectrum {
        precursors: vec![Precursor {
            mz: 250.0,
            isolation_window_lower_offset: 50.0,
            isolation_window_upper_offset: 50.0,
            ..Precursor::default()
        }],
        ..MSSpectrum::default()
    });
    let (lower, upper, center) = helper::check_swath_map(&map).unwrap();
    real_similar(lower, 200.0);
    real_similar(upper, 300.0);
    real_similar(center, 250.0);
    // Derived: the sanity checks.
    assert!(helper::check_swath_map(&MSExperiment::new()).is_err());
    let mut other = map.clone();
    other.spectra.push(MSSpectrum {
        ms_level: 2,
        ..map.spectra[0].clone()
    });
    assert!(helper::check_swath_map(&other).is_err());
    let mut shifted = map.clone();
    let mut s = map.spectra[0].clone();
    s.precursors[0].mz = 250.2;
    shifted.spectra.push(s);
    assert!(helper::check_swath_map(&shifted).is_err());
    let mut exp = LightTargetedExperiment::new();
    exp.transitions.push(LightTransition {
        precursor_mz: 260.0,
        ..LightTransition::default()
    });
    let mut out = LightTargetedExperiment::new();
    assert!(helper::check_swath_map_and_select_transitions(&map, &exp, &mut out, 1.0).unwrap());
    assert!(
        !helper::check_swath_map_and_select_transitions(
            &MSExperiment::new(),
            &exp,
            &mut LightTargetedExperiment::new(),
            1.0
        )
        .unwrap()
    );
}

/// OpenSwathHelper_test.cpp: estimateRTRange.
#[test]
fn estimate_rt_range() {
    let mut exp = LightTargetedExperiment::new();
    for rt in [-100.0, 900.0, 300.0] {
        exp.compounds.push(LightCompound {
            rt,
            ..LightCompound::default()
        });
    }
    let (min, max) = helper::estimate_rt_range(&exp).unwrap();
    real_similar(min, -100.0);
    real_similar(max, 900.0);
    assert!(helper::estimate_rt_range(&LightTargetedExperiment::new()).is_err());
}

/// Derived from OpenSwathHelper.cpp: the PASEF window choice.
#[test]
fn select_swath_transitions_pasef() {
    let mut exp = LightTargetedExperiment::new();
    for (mz, im) in [(450.0, 0.95), (450.0, 1.3), (700.0, 1.0)] {
        exp.transitions.push(LightTransition {
            precursor_mz: mz,
            precursor_im: im,
            ..LightTransition::default()
        });
    }
    let maps = [
        SwathMap::with_im(400.0, 500.0, 450.0, 0.8, 1.2, false),
        SwathMap::with_im(400.0, 500.0, 450.0, 0.9, 1.1, false),
        SwathMap::with_im(400.0, 500.0, 450.0, 1.1, 1.5, false),
    ];
    let map = helper::select_swath_transitions_pasef(&exp, 1.0, &maps);
    // IM 0.95: windows 0 (centre 1.0) and 1 (centre 1.0) tie; the first stays.
    assert_eq!(map, vec![Some(0), Some(2), None]);
}

/// Derived from OpenSwathHelper.cpp: the best feature per group, keyed by the
/// first transition's compound reference, with the quality cutoff.
#[test]
fn simple_find_best_feature() {
    let mut groups: BTreeMap<String, MRMTransitionGroup<MSChromatogram, LightTransition>> =
        BTreeMap::new();
    let mut g = MRMTransitionGroup::new();
    g.add_transition(
        LightTransition {
            transition_name: "t".into(),
            peptide_ref: "PEP".into(),
            ..LightTransition::default()
        },
        "t",
    )
    .unwrap();
    for (rt, q) in [(10.0, 0.2f32), (20.0, 0.9), (30.0, 0.9)] {
        let mut f = MRMFeature::new();
        f.rt = rt;
        f.feature.base.quality = q;
        g.add_feature(f).unwrap();
    }
    groups.insert("g1".into(), g);
    groups.insert("empty".into(), MRMTransitionGroup::new());
    let best = helper::simple_find_best_feature(&groups, false, 0.0).unwrap();
    assert_eq!(best.len(), 1);
    assert_eq!(best["PEP"], 20.0);
    assert!(
        helper::simple_find_best_feature(&groups, true, 0.95)
            .unwrap()
            .is_empty()
    );
}

// ---------------------------------------------------------------------------
// SpectrumAddition
// ---------------------------------------------------------------------------

const INT_SORTED: [f64; 11] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
const MZ_SORTED: [f64; 11] = [
    100.0, 101.5, 101.9, 102.0, 102.1, 102.11, 102.2, 102.25, 102.3, 102.4, 102.45,
];
const INT_UNSORTED: [f64; 11] = [11.0, 4.0, 3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 10.0];
const MZ_UNSORTED: [f64; 11] = [
    102.45, 102.0, 101.9, 102.1, 102.11, 102.2, 102.25, 102.3, 100.0, 101.5, 102.4,
];

/// SpectrumAddition_test.cpp: sortSpectrumByMZ, without and with IM.
#[test]
fn sort_spectrum_by_mz() {
    let mut spec = Spectrum::from_arrays(MZ_UNSORTED.to_vec(), INT_UNSORTED.to_vec());
    addition::sort_spectrum_by_mz(&mut spec).unwrap();
    for i in 0..11 {
        real_similar(MZ_SORTED[i], spec.mz_array().data[i]);
        real_similar(INT_SORTED[i], spec.intensity_array().data[i]);
    }
    let mut spec = Spectrum::from_arrays(MZ_UNSORTED.to_vec(), INT_UNSORTED.to_vec());
    spec.set_drift_time_array(BinaryDataArray::new(INT_UNSORTED.to_vec()));
    addition::sort_spectrum_by_mz(&mut spec).unwrap();
    let im = &spec.drift_time_array().unwrap().data;
    assert_eq!(im.len(), 11);
    for i in 0..11 {
        real_similar(MZ_SORTED[i], spec.mz_array().data[i]);
        real_similar(INT_SORTED[i], spec.intensity_array().data[i]);
        real_similar(INT_SORTED[i], im[i]);
    }
    // Derived: a misaligned array is refused, an empty one skipped.
    let mut bad = Spectrum::from_arrays(vec![2.0, 1.0], vec![1.0]);
    assert!(addition::sort_spectrum_by_mz(&mut bad).is_err());
    let mut skip = Spectrum::from_arrays(vec![2.0, 1.0], vec![]);
    addition::sort_spectrum_by_mz(&mut skip).unwrap();
    assert_eq!(skip.mz_array().data, vec![1.0, 2.0]);
}

const INT1: [f64; 11] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0];
const INT2: [f64; 11] = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 9.0, 7.0, 5.0, 3.0, 1.0];
const MZ1: [f64; 11] = [
    100.0, 101.5, 101.9, 102.0, 102.1, 102.11, 102.2, 102.25, 102.3, 102.4, 102.45,
];
const MZ2: [f64; 11] = [
    100.0, 101.6, 101.95, 102.0, 102.05, 102.1, 102.12, 102.15, 102.2, 102.25, 102.30,
];

/// SpectrumAddition_test.cpp: addUpSpectra over OpenSWATH spectra.
#[test]
fn add_up_spectra_light() {
    let empty = addition::add_up_spectra(&[], 0.1, false).unwrap();
    assert_eq!(empty.mz_array().data.len(), 0);
    let two_empty = [Arc::new(Spectrum::new()), Arc::new(Spectrum::new())];
    assert_eq!(
        addition::add_up_spectra(&two_empty, 0.1, false)
            .unwrap()
            .len(),
        0
    );
    let all = [
        Arc::new(Spectrum::from_arrays(MZ1.to_vec(), INT1.to_vec())),
        Arc::new(Spectrum::from_arrays(MZ2.to_vec(), INT2.to_vec())),
    ];
    let result = addition::add_up_spectra(&all, 0.1, false).unwrap();
    assert_eq!(result.mz_array().data.len(), 25);
    let filtered = addition::add_up_spectra(&all, 0.1, true).unwrap();
    assert_eq!(filtered.mz_array().data.len(), 9);
    real_similar(filtered.mz_array().data[0], 100.0);
    real_similar(filtered.intensity_array().data[0], 2.0);
    real_similar(filtered.mz_array().data[3], 101.9);
    real_similar(filtered.intensity_array().data[3], 3.0 + 5.0 / 2.0);
    // Derived: total intensity is conserved; a single spectrum is returned as is.
    let total: f64 = result.intensity_array().data.iter().sum();
    real_similar(total, INT1.iter().sum::<f64>() + INT2.iter().sum::<f64>());
    let single = addition::add_up_spectra(&all[..1], 0.1, true).unwrap();
    assert!(Arc::ptr_eq(&single, &all[0]));
    assert!(addition::add_up_spectra(&all, 0.0, false).is_err());
    assert!(addition::add_up_spectra(&all, 1e-12, false).is_err());
}

/// SpectrumAddition_test.cpp: addUpSpectra over MSSpectrum.
#[test]
fn add_up_spectra_ms() {
    let s = |mz: &[f64], int: &[f64]| MSSpectrum {
        peaks: mz
            .iter()
            .zip(int)
            .map(|(&m, &i)| Peak1D::new(m, i as f32))
            .collect(),
        ..MSSpectrum::default()
    };
    assert!(
        addition::add_up_ms_spectra(&[], 0.1, false)
            .unwrap()
            .peaks
            .is_empty()
    );
    let two_empty = [MSSpectrum::default(), MSSpectrum::default()];
    assert_eq!(
        addition::add_up_ms_spectra(&two_empty, 0.1, false)
            .unwrap()
            .peaks
            .len(),
        0
    );
    let all = [s(&MZ1, &INT1), s(&MZ2, &INT2)];
    let result = addition::add_up_ms_spectra(&all, 0.1, false).unwrap();
    assert_eq!(result.peaks.len(), 25);
    let filtered = addition::add_up_ms_spectra(&all, 0.1, true).unwrap();
    assert_eq!(filtered.peaks.len(), 9);
    real_similar(filtered.peaks[0].mz, 100.0);
    real_similar(f64::from(filtered.peaks[0].intensity), 2.0);
    real_similar(filtered.peaks[3].mz, 101.9);
    real_similar(f64::from(filtered.peaks[3].intensity), 3.0 + 5.0 / 2.0);
    let fine = addition::add_up_ms_spectra(&all, 0.01, true).unwrap();
    let auto = addition::add_up_ms_spectra(&all, -1.0, true).unwrap();
    assert_eq!(fine.peaks.len(), 16);
    real_similar(fine.peaks[0].mz, 100.0);
    real_similar(f64::from(fine.peaks[0].intensity), 2.0);
    real_similar(fine.peaks[3].mz, 101.9);
    real_similar(f64::from(fine.peaks[3].intensity), 3.0);
    assert_eq!(auto.peaks.len(), 28);
}

/// Derived from SpectrumAddition.cpp: the ion-mobility variant filters first;
/// concatenation keeps extra arrays and sorts.
#[test]
fn add_up_in_mobility_range_and_concatenate() {
    let mk = |mz: Vec<f64>, int: Vec<f64>, im: Vec<f64>| {
        let mut s = Spectrum::from_arrays(mz, int);
        s.set_drift_time_array(BinaryDataArray::new(im));
        Arc::new(s)
    };
    let all = [
        mk(vec![100.0, 101.0], vec![1.0, 2.0], vec![0.5, 1.5]),
        mk(vec![100.0, 101.0], vec![3.0, 4.0], vec![0.6, 0.7]),
    ];
    let range = RangeBase::from_min_max(0.4, 1.0).unwrap();
    let summed = addition::add_up_spectra_in_mobility_range(&all, &range, 1.0, true).unwrap();
    assert_eq!(summed.mz_array().data, vec![100.0, 101.0]);
    assert_eq!(summed.intensity_array().data, vec![4.0, 4.0]);
    let unfiltered =
        addition::add_up_spectra_in_mobility_range(&all, &RangeBase::new(), 1.0, true).unwrap();
    assert_eq!(unfiltered.intensity_array().data, vec![4.0, 6.0]);
    let plain = Arc::new(Spectrum::from_arrays(vec![1.0], vec![1.0]));
    let no_im = vec![Arc::clone(&plain), plain];
    assert!(addition::add_up_spectra_in_mobility_range(&no_im, &range, 1.0, true).is_err());

    let cat = addition::concatenate_spectra(&all).unwrap();
    assert_eq!(cat.mz_array().data, vec![100.0, 100.0, 101.0, 101.0]);
    assert_eq!(cat.intensity_array().data, vec![1.0, 3.0, 2.0, 4.0]);
    assert_eq!(
        cat.drift_time_array().unwrap().data,
        vec![0.5, 0.6, 1.5, 0.7]
    );
    let wider = [Arc::new(Spectrum::new()), all[0].clone()];
    assert!(addition::concatenate_spectra(&wider).is_err());
    assert!(addition::concatenate_spectra(&[]).unwrap().is_empty());
}

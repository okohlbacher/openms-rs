// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (evidence tier 1) for the OpenSWATH scoring port:
//! DIAHelper, DIAPrescoring, DIAScoring, MRMScoring, IonMobilityScoring and
//! OpenSwathScoring against the C++ Release build.
//!
//! The fixture `tests/data/openswath_scoring/oracle_all.stdout` is the output
//! of the driver `../oracle/O2-scoring/driver.cpp`, built with
//! x86_64-conda-linux-gnu-g++ 14.4.0 `-O3 -DNDEBUG` against the Linux x86-64
//! Release install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
//! ibminode06 with byte-identical output. Every record is `key<TAB>%.17g`; this
//! test recomputes each key with the port and compares. The tolerance classes
//! are declared in `docs/OPENSWATH_SCORING_SUPPORT.md`.

use openms::analysis::openswath::data_structures::{BinaryDataArray, Spectrum, SpectrumPtr};
use openms::analysis::openswath::dia_helper as dh;
use openms::analysis::openswath::dia_prescoring::DiaPrescore;
use openms::analysis::openswath::dia_scoring::DiaScoring;
use openms::analysis::openswath::ion_mobility_scoring::{DriftExtraction, IonMobilityScoring};
use openms::analysis::openswath::itransition::SignalToNoisePtr;
use openms::analysis::openswath::mock_objects::{MockFeature, MockMRMFeature, MockSignalToNoise};
use openms::analysis::openswath::mrm_scoring::MRMScoring;
use openms::analysis::openswath::openswath_scores::OpenSwathScoresUsage;
use openms::analysis::openswath::openswath_scoring::OpenSwathScoring;
use openms::analysis::openswath::spectrum_access::SpectrumAccessPtr;
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::swath_map::SwathMap;
use openms::analysis::openswath::transition_experiment::{LightCompound, LightTransition};
use openms::chemistry::{
    AASequence, EmpiricalFormula, PeptideFragmentType, TheoreticalSpectrumGenerator,
};
use openms::kernel::mrm::{OpenSwathIndScores, OpenSwathScores};
use openms::kernel::ranges::RangeBase;
use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
use openms::param::ParamValue;
use std::collections::BTreeMap;
use std::sync::Arc;

const ORACLE: &str = include_str!("data/openswath_scoring/oracle_all.stdout");
const ORACLE_ERR: &str = include_str!("data/openswath_scoring/oracle_all.stderr");

/// Records `key -> value`, in the order the port produces them.
#[derive(Default)]
struct Rec(BTreeMap<String, f64>);

impl Rec {
    fn p(&mut self, k: impl Into<String>, v: f64) {
        let k = k.into();
        assert!(self.0.insert(k.clone(), v).is_none(), "duplicate key {k}");
    }
    fn pv(&mut self, k: &str, v: &[f64]) {
        for (i, x) in v.iter().enumerate() {
            self.p(format!("{k}.{i}"), *x);
        }
        self.p(format!("{k}.n"), v.len() as f64);
    }
    fn pairs(&mut self, k: &str, v: &[(f64, f64)]) {
        for (j, (mz, int)) in v.iter().enumerate() {
            self.p(format!("{k}.{j}.mz"), *mz);
            self.p(format!("{k}.{j}.int"), *int);
        }
    }
}

fn oracle() -> BTreeMap<String, String> {
    ORACLE
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

/// Keys compared with relative tolerance 1e-12 instead of exactly (see
/// `docs/OPENSWATH_SCORING_SUPPORT.md`, "Tolerance").
///
/// Only the fragment m/z values of the chemistry port's
/// `TheoreticalSpectrumGenerator` remain (`by2.*`, `theor1.*`): it adds the
/// protons after the residues where the source starts from
/// `PROTON_MASS_U * charge` and adds residues
/// (`TheoreticalSpectrumGenerator.cpp:724-794`), which moves the last bits;
/// outside this port's scope and reported to the chemistry owner.
///
/// The Eigen reductions of `Scoring.cpp` and `StatsHelpers.cpp` (every `shape`
/// score, `prescore*`, `with_iso*`, `lib.dot`, `lib.man`, `*dotprod_score_dia`,
/// `*manhatt_score_dia`) were a second class here while the port added left to
/// right; the port now adds in Eigen's own order and they are compared exactly.
fn reordered(key: &str) -> bool {
    ["by2.", "theor1."].iter().any(|p| key.starts_with(p))
}

fn compare(port: &Rec, prefix: &str) {
    let oracle = oracle();
    let mut checked = 0;
    let mut failures = Vec::new();
    for (key, text) in oracle.iter().filter(|(k, _)| k.starts_with(prefix)) {
        let Ok(expected) = text.parse::<f64>() else {
            continue;
        };
        let Some(&got) = port.0.get(key) else {
            failures.push(format!("{key}: missing in port (oracle {text})"));
            continue;
        };
        let ok = if expected.is_nan() {
            got.is_nan()
        } else if reordered(key) {
            (got - expected).abs() <= 1e-12 * expected.abs().max(1e-300)
        } else {
            got == expected
        };
        if !ok {
            failures.push(format!("{key}: port {got:.17e} oracle {expected:.17e}"));
        }
        checked += 1;
    }
    for key in port.0.keys().filter(|k| k.starts_with(prefix)) {
        if !oracle.contains_key(key) {
            failures.push(format!("{key}: not in oracle"));
        }
    }
    assert!(checked > 0, "no oracle records for {prefix}");
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn mk(mz: &[f64], int: &[f64], im: Option<&[f64]>) -> SpectrumPtr {
    let mut s = Spectrum::from_arrays(mz.to_vec(), int.to_vec());
    if let Some(im) = im {
        s.push_data_array(BinaryDataArray::with_description(
            im.to_vec(),
            "Ion Mobility",
        ));
    }
    Arc::new(s)
}

fn rng(c: f64, w: f64) -> RangeBase {
    let mut r = RangeBase::singular(c).unwrap();
    r.min_span_if_singular(w).unwrap();
    r
}

const SMZ: [f64; 48] = [
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03, 499.97, 499.98, 499.99, 500.0, 500.01,
    500.02, 500.03, 500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03, 501.99, 502.0, 502.01,
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03, 600.97, 600.98, 600.99, 601.0, 601.01,
    601.02, 601.03, 601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03, 602.99, 603.0, 603.01,
];
const SINT: [f64; 48] = [
    10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 15., 7., 3., 1., 3., 9., 15., 9., 3., 1.,
    3., 9., 3., 10., 20., 50., 100., 50., 20., 10., 3., 7., 15., 30., 15., 7., 3., 1., 3., 9., 15.,
    9., 3., 1., 3., 9., 3.,
];

fn tr(n: &str, pmz: f64, ch: i8, lib: f64, prec: f64) -> LightTransition {
    LightTransition {
        transition_name: n.into(),
        product_mz: pmz,
        fragment_charge: ch,
        library_intensity: lib,
        precursor_mz: prec,
        ..LightTransition::default()
    }
}

fn mf(iv: &[f64], apex: f32, rt: f64) -> Arc<MockFeature> {
    Arc::new(MockFeature {
        intensity_vec: iv.to_vec(),
        rt_vec: vec![0.0; iv.len()],
        intensity: apex,
        rt,
    })
}

fn dia(window: f64) -> DiaScoring {
    let mut p = DiaScoring::defaults().unwrap();
    p.set_value("dia_extraction_window", ParamValue::Float(window), "", &[])
        .unwrap();
    p.set_value(
        "dia_byseries_intensity_min",
        ParamValue::Float(30.0),
        "",
        &[],
    )
    .unwrap();
    p.set_value("dia_byseries_ppm_diff", ParamValue::Float(50.0), "", &[])
        .unwrap();
    DiaScoring::from_param(&p).unwrap()
}

/// driver.cpp `helper_case`.
#[test]
fn oracle_dia_helper() {
    let mut r = Rec::default();
    let cases = [
        (100., 1, 4, 1.00048),
        (30., 1, 4, 1.00048),
        (110., 1, 4, 1.00048),
        (500., 2, 4, 1.00048),
        (500., -3, 6, 1.0033548378),
        (1500.5, 1, 0, 1.00048),
    ];
    for (i, (mz, ch, n, mann)) in cases.iter().enumerate() {
        let v = dh::averagine_isotope_distribution(*mz, *ch, *n, *mann).unwrap();
        r.pairs(&format!("averagine.{i}"), &v);
        r.p(format!("averagine.{i}.n"), v.len() as f64);
    }
    let mut v = Vec::new();
    dh::add_single_peak_isotopes_to_spec(500.2, 7.5, &mut v, 5, 2).unwrap();
    r.pairs("single", &v);
    r.p("single.n", v.len() as f64);
    let mut out = Vec::new();
    dh::add_isotopes_to_spec(&[(100., 100.), (200., 300.)], &mut out, 3, 1).unwrap();
    r.pairs("iso2spec", &out);
    r.p("iso2spec.n", out.len() as f64);
    let mut spec = vec![(501., 1.)];
    dh::add_preisotope_weights(&[500., 600.5], &mut spec, 2, -0.5, 1.000482, -2).unwrap();
    r.pairs("preiso", &spec);
    r.p("preiso.n", spec.len() as f64);
    let a = dh::create_mz_range_ppm(500., 20., true).unwrap();
    r.p("range_ppm.min", a.min().unwrap());
    r.p("range_ppm.max", a.max().unwrap());
    let b = dh::create_mz_range_ppm(500., 0.05, false).unwrap();
    r.p("range_th.min", b.min().unwrap());
    r.p("range_th.max", b.max().unwrap());
    let g = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..TheoreticalSpectrumGenerator::default()
    };
    let seq = AASequence::parse("SYVAWDR").unwrap();
    let (bs, ys) = dh::by_series(&seq, &g, 2).unwrap();
    r.pv("by2.b", &bs);
    r.pv("by2.y", &ys);
    r.pv("theor1", &dh::theoretical_masses(&seq, &g, 1).unwrap());
    for prefix in [
        "averagine",
        "single",
        "iso2spec",
        "preiso",
        "range_",
        "by2",
        "theor1",
    ] {
        compare(&r, prefix);
    }
}

/// driver.cpp `diascoring_case`.
#[test]
fn oracle_dia_scoring() {
    let mut r = Rec::default();
    let seq = vec![mk(&SMZ, &SINT, None)];
    let d = dia(0.05);
    let none = RangeBase::new();
    let mut features = BTreeMap::new();
    features.insert("group1".to_string(), mf(&[], 0.3, 0.0));
    features.insert("group2".to_string(), mf(&[], 0.7, 0.0));
    let f = MockMRMFeature {
        features,
        intensity: 1.0,
        ..MockMRMFeature::default()
    };
    let ts = [tr("group1", 500., 1, 5., 0.), tr("group2", 600., 1, 5., 0.)];
    let s = d.dia_isotope_scores(&ts, &seq, &f, &none).unwrap();
    r.p("iso.corr", s.isotope_corr);
    r.p("iso.overlap", s.isotope_overlap);
    let s = d
        .dia_isotope_scores(&[tr("group2", 600., -2, 5., 0.)], &seq, &f, &none)
        .unwrap();
    r.p("iso_neg2.corr", s.isotope_corr);
    r.p("iso_neg2.overlap", s.isotope_overlap);
    for (k, mz, ch) in [("500_1", 500., 1), ("500_2", 500., 2), ("499_1", 499., 1)] {
        let s = d
            .dia_ms1_isotope_scores_averagine(mz, &seq, ch, &none)
            .unwrap();
        r.p(format!("ms1avg.{k}.corr"), s.isotope_corr);
        r.p(format!("ms1avg.{k}.overlap"), s.isotope_overlap);
    }
    let glucose = EmpiricalFormula::parse("C6H12O6").unwrap();
    let s = d
        .dia_ms1_isotope_scores(499., &seq, &none, &glucose.clone().with_charge(1))
        .unwrap();
    r.p("ms1formula.glucose1.corr", s.isotope_corr);
    r.p("ms1formula.glucose1.overlap", s.isotope_overlap);
    let s = d
        .dia_ms1_isotope_scores(499., &seq, &none, &glucose)
        .unwrap();
    r.p("ms1formula.glucose0.corr", s.isotope_corr);
    r.p("ms1formula.glucose0.overlap", s.isotope_overlap);
    let pep = AASequence::parse("PEPTIDEK")
        .unwrap()
        .formula_for(PeptideFragmentType::Full, 2)
        .unwrap();
    // The source's EmpiricalFormula::toString omits the charge; getFormula(Full,
    // charge) sets it as metadata without adding hydrogens, as the port does.
    assert_eq!(pep.charge(), 2);
    assert_eq!(
        pep.clone().with_charge(0).to_string(),
        oracle()["pep_formula"],
        "peptide formula differs from the oracle"
    );
    let s = d.dia_ms1_isotope_scores(500., &seq, &none, &pep).unwrap();
    r.p("ms1formula.pep2.corr", s.isotope_corr);
    r.p("ms1formula.pep2.overlap", s.isotope_overlap);
    let w = d.score_with_isotopes(&seq, &ts, &none).unwrap();
    r.p("with_iso.dot", w.dotprod);
    r.p("with_iso.man", w.manhattan);
    let mut smz = SMZ.to_vec();
    let half = smz.len() as f64 / 2.0;
    for (i, v) in smz.iter_mut().enumerate() {
        let ppm = if (i as f64) < half { 15.0 } else { 10.0 };
        *v += *v / 1_000_000.0 * ppm;
    }
    let sh = vec![mk(&smz, &SINT, None)];
    let dl = dia(0.5);
    let m = dl.dia_massdiff_score(&ts, &sh, &[0.7, 0.3], &none).unwrap();
    r.p("massdiff.ppm", m.ppm_score);
    r.p("massdiff.ppmw", m.ppm_score_weighted);
    r.pv("massdiff.diff", &m.diff_ppm);
    let x = dl.dia_ms1_massdiff_score(500., &sh, &none).unwrap();
    r.p("ms1diff.500", x.ppm_score);
    r.p("ms1diff.500.found", f64::from(u8::from(x.signal_found)));
    let x = dl.dia_ms1_massdiff_score(100., &sh, &none).unwrap();
    r.p("ms1diff.100", x.ppm_score);
    r.p("ms1diff.100.found", f64::from(u8::from(x.signal_found)));
    let mut dppm = dia(20.0);
    dppm.dia_extraction_ppm = true;
    let x = dppm.dia_ms1_massdiff_score(600.004, &sh, &none).unwrap();
    r.p("ms1diff_ppm.600", x.ppm_score);
    r.p("ms1diff_ppm.600.found", f64::from(u8::from(x.signal_found)));
    let bmz = [
        350.17164,
        421.20875,
        421.20875 + 79.9657,
        547.26291,
        646.33133,
        809.39466 + 79.9657,
    ];
    let bs = vec![mk(&bmz, &[100.0; 6], None)];
    let by = d
        .dia_by_ion_score(&bs, &AASequence::parse("SYVAWDR").unwrap(), 1, &none)
        .unwrap();
    r.p("by.b", by.bseries_score);
    r.p("by.y", by.yseries_score);
    let p = DiaPrescore::with_window(0.05, 4, 4)
        .score(&seq, &ts, &none)
        .unwrap();
    r.p("prescore.dot", p.dotprod);
    r.p("prescore.man", p.manhattan);
    let p = DiaPrescore::with_window(0.1, 0, 4)
        .score(&seq, &ts, &none)
        .unwrap();
    r.p("prescore_all.dot", p.dotprod);
    r.p("prescore_all.man", p.manhattan);
    for prefix in [
        "iso.",
        "iso_neg2",
        "ms1avg",
        "ms1formula",
        "with_iso",
        "massdiff",
        "ms1diff",
        "by.",
        "prescore",
    ] {
        compare(&r, prefix);
    }
    // The only stderr line is the source's once-per-process deprecation
    // warning of CoarseIsotopePatternGenerator::run for a charged formula.
    assert!(ORACLE_ERR.starts_with("Warning: CoarseIsotopePatternGenerator::run()"));
}

fn mrm_feature() -> MockMRMFeature {
    let i1 = [
        5.97543668746948,
        4.2749171257019,
        3.3301842212677,
        4.08597040176392,
        5.50307035446167,
        5.24326848983765,
        8.40812492370605,
        2.83419919013977,
        6.94378805160522,
        7.69957494735718,
        4.08597040176392,
    ];
    let i2 = [
        15.8951349258423,
        41.5446395874023,
        76.0746307373047,
        109.069435119629,
        111.90364074707,
        169.79216003418,
        121.043930053711,
        63.0136985778809,
        44.6150207519531,
        21.4926776885986,
        7.93575811386108,
    ];
    let m1 = [
        0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0,
    ];
    let m2 = [
        10.0, 115.0, 180.0, 280.0, 330.0, 340.0, 390.0, 320.0, 300.0, 250.0, 100.0,
    ];
    let mut features = BTreeMap::new();
    features.insert("group1".to_string(), mf(&i1, 782.38073, 1200.));
    features.insert("group2".to_string(), mf(&i2, 58.384506, 1200.));
    let mut pre = BTreeMap::new();
    pre.insert("ms1trace1".to_string(), mf(&m1, 0.0, 1200.));
    pre.insert("ms1trace2".to_string(), mf(&m2, 0.0, 1200.));
    pre.insert("ms1trace3".to_string(), mf(&[0.0; 11], 0.0, 1200.));
    MockMRMFeature {
        features,
        precursor_features: pre,
        intensity: 1.0,
        rt: 1200.0,
        metavalue: 0.0,
    }
}

fn ids(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| (*s).to_string()).collect()
}

/// driver.cpp `mrm_case`.
#[test]
fn oracle_mrm_scoring() {
    let mut r = Rec::default();
    let f = mrm_feature();
    let nat = ids(&["group1", "group2"]);
    let pre = ids(&["ms1trace1", "ms1trace2", "ms1trace3"]);
    let mut s = MRMScoring::new();
    s.initialize_xcorr_matrix(&f, &nat).unwrap();
    r.p("xcorr.coelution", s.calc_xcorr_coelution_score());
    r.p("xcorr.shape", s.calc_xcorr_shape_score());
    r.p(
        "xcorr.coelution_w",
        s.calc_xcorr_coelution_weighted_score(&[0.3, 0.7]).unwrap(),
    );
    r.p(
        "xcorr.shape_w",
        s.calc_xcorr_shape_weighted_score(&[0.3, 0.7]).unwrap(),
    );
    s.initialize_xcorr_contrast_matrix(&f, &ids(&["group1"]), &nat)
        .unwrap();
    r.p("contrast.shape_sum", s.calc_xcorr_contrast_shape_score());
    r.pv(
        "contrast.sep_coel",
        &s.calc_separate_xcorr_contrast_coelution_score(),
    );
    r.pv(
        "contrast.sep_shape",
        &s.calc_separate_xcorr_contrast_shape_score(),
    );
    s.initialize_xcorr_precursor_matrix(&f, &pre).unwrap();
    r.p("pre.coelution", s.calc_xcorr_precursor_coelution_score());
    r.p("pre.shape", s.calc_xcorr_precursor_shape_score());
    s.initialize_xcorr_precursor_contrast_matrix(&f, &pre, &nat)
        .unwrap();
    r.p(
        "precon.coelution",
        s.calc_xcorr_precursor_contrast_coelution_score(),
    );
    r.p(
        "precon.shape",
        s.calc_xcorr_precursor_contrast_shape_score(),
    );
    s.initialize_xcorr_precursor_combined_matrix(&f, &pre, &nat)
        .unwrap();
    r.p(
        "precomb.coelution",
        s.calc_xcorr_precursor_combined_coelution_score(),
    );
    r.p(
        "precomb.shape",
        s.calc_xcorr_precursor_combined_shape_score(),
    );
    s.initialize_mi_matrix(&f, &nat).unwrap();
    r.p("mi.score", s.calc_mi_score());
    r.p("mi.score_w", s.calc_mi_weighted_score(&[0.3, 0.7]).unwrap());
    s.initialize_mi_precursor_matrix(&f, &pre).unwrap();
    r.p("mipre.score", s.calc_mi_precursor_score());
    s.initialize_mi_precursor_contrast_matrix(&f, &pre, &nat)
        .unwrap();
    r.p("miprecon.score", s.calc_mi_precursor_contrast_score());
    s.initialize_mi_precursor_combined_matrix(&f, &pre, &nat)
        .unwrap();
    r.p("miprecomb.score", s.calc_mi_precursor_combined_score());
    s.initialize_mi_contrast_matrix(&f, &nat, &ids(&["group2", "group1"]))
        .unwrap();
    r.pv("micon.sep", &s.calc_separate_mi_contrast_score());
    let ts = [
        tr("group1", 500., 1, 1., 0.),
        tr("group2", 600., 1, 10000., 0.),
    ];
    let lib = MRMScoring::calc_library_score(&f, &ts).unwrap();
    r.p("lib.corr", lib.correlation);
    r.p("lib.norm_man", lib.norm_manhattan);
    r.p("lib.man", lib.manhattan);
    r.p("lib.dot", lib.dotprod);
    r.p("lib.sangle", lib.spectral_angle);
    r.p("lib.rmsd", lib.rmsd);
    let sn: Vec<SignalToNoisePtr> = vec![
        Arc::new(MockSignalToNoise { sn_value: 0.5 }),
        Arc::new(MockSignalToNoise { sn_value: 37.25 }),
    ];
    r.p("sn.mean", MRMScoring::calc_sn_score(&f, &sn));
    r.pv("sn.sep", &MRMScoring::calc_separate_sn_score(&f, &sn));
    for prefix in [
        "xcorr.",
        "contrast.",
        "pre.",
        "precon.",
        "precomb.",
        "mi.",
        "mipre.",
        "miprecon.",
        "miprecomb.",
        "micon.",
        "lib.",
        "sn.",
    ] {
        compare(&r, prefix);
    }
}

fn imspec() -> SpectrumPtr {
    let mut mz = vec![500.2; 8];
    mz.extend([500.3; 8]);
    mz.extend([600.2, 600.3, 600.4, 600.5, 600.6, 600.7, 600.8, 600.9]);
    let mut int = vec![10., 20., 30., 40., 40., 30., 20., 10.];
    int.extend([10., 20., 30., 40., 40., 30., 20., 10.]);
    int.extend([20.; 8]);
    let im = [
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 0.1, 0.2,
        0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    ];
    mk(&mz, &int, Some(&im))
}

fn ms1spec() -> SpectrumPtr {
    mk(
        &[700.2, 700.3, 700.4, 700.5, 700.6, 700.7, 700.8, 700.9],
        &[10., 20., 30., 40., 40., 30., 20., 10.],
        Some(&[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
    )
}

fn im_scores(r: &mut Rec, k: &str, s: &OpenSwathScores) {
    for (name, v) in [
        ("im_drift", s.im_drift),
        ("im_drift_weighted", s.im_drift_weighted),
        ("im_delta_score", s.im_delta_score),
        ("im_delta", s.im_delta),
        ("im_log_intensity", s.im_log_intensity),
        ("im_xcorr_coelution_score", s.im_xcorr_coelution_score),
        ("im_xcorr_shape_score", s.im_xcorr_shape_score),
        ("im_ms1_drift", s.im_ms1_drift),
        ("im_ms1_delta_score", s.im_ms1_delta_score),
        ("im_ms1_delta", s.im_ms1_delta),
        ("im_ms1_contrast_coelution", s.im_ms1_contrast_coelution),
        ("im_ms1_contrast_shape", s.im_ms1_contrast_shape),
        (
            "im_ms1_sum_contrast_coelution",
            s.im_ms1_sum_contrast_coelution,
        ),
        ("im_ms1_sum_contrast_shape", s.im_ms1_sum_contrast_shape),
        ("im_ind_contrast_coelution", s.im_ind_contrast_coelution),
        ("im_ind_contrast_shape", s.im_ind_contrast_shape),
        (
            "im_ind_sum_contrast_coelution",
            s.im_ind_sum_contrast_coelution,
        ),
        ("im_ind_sum_contrast_shape", s.im_ind_sum_contrast_shape),
    ] {
        r.p(format!("{k}.{name}"), v);
    }
}

fn ext(window: f64, ppm: bool, extra: f64) -> DriftExtraction {
    DriftExtraction {
        dia_extraction_window: window,
        dia_extraction_ppm: ppm,
        drift_extra: extra,
    }
}

/// driver.cpp `im_case`.
#[test]
fn oracle_ion_mobility_scoring() {
    let mut r = Rec::default();
    let sp = vec![imspec()];
    let ms1 = vec![ms1spec()];
    let ts = [
        tr("group1", 500.2, 1, 1., 700.2),
        tr("group2", 600.5, 1, 1., 700.2),
    ];
    let mut s = OpenSwathScores::default();
    let e = ext(0.3, false, 0.25);
    IonMobilityScoring::drift_scoring(&sp, &ts, &mut s, 1.0, rng(1.0, 1.0), &e, false, None)
        .unwrap();
    IonMobilityScoring::drift_scoring_ms1(&ms1, &ts, &mut s, 1.0, rng(1.0, 1.0), &e, None).unwrap();
    IonMobilityScoring::drift_scoring_ms1_contrast(&sp, &ms1, &ts, &mut s, rng(1.0, 1.0), &e, None)
        .unwrap();
    im_scores(&mut r, "im_a", &s);
    let mut s = OpenSwathScores::default();
    IonMobilityScoring::drift_scoring(
        &sp,
        &ts,
        &mut s,
        0.9,
        rng(0.8, 0.6),
        &ext(20.0, true, 0.1),
        false,
        None,
    )
    .unwrap();
    IonMobilityScoring::drift_scoring_ms1_contrast(
        &sp,
        &ms1,
        &ts,
        &mut s,
        rng(0.8, 0.6),
        &ext(400.0, true, 0.1),
        None,
    )
    .unwrap();
    im_scores(&mut r, "im_b", &s);
    let mut s = OpenSwathScores::default();
    IonMobilityScoring::drift_id_scoring(
        &sp,
        &[tr("id1", 500.3, 1, 1., 700.2)],
        &ts,
        &mut s,
        1.0,
        rng(1.0, 1.0),
        &e,
        false,
        None,
    )
    .unwrap();
    im_scores(&mut r, "im_id", &s);
    let m =
        IonMobilityScoring::compute_ion_mobilogram(&sp, &rng(500.25, 0.3), &rng(1.0, 1.5), 1.5e-5)
            .unwrap();
    r.p("mobilogram.im", m.im);
    r.p("mobilogram.int", m.intensity);
    r.p("mobilogram.n", m.mobilogram.peaks.len() as f64);
    for (i, p) in m.mobilogram.peaks.iter().enumerate() {
        r.p(format!("mobilogram.{i}.mob"), p.mobility);
        r.p(format!("mobilogram.{i}.int"), f64::from(p.intensity));
    }
    for prefix in ["im_a.", "im_b.", "im_id.", "mobilogram."] {
        compare(&r, prefix);
    }
}

fn msspec(mz: &[f64], int: &[f64], rt: f64, level: u32) -> MSSpectrum {
    MSSpectrum {
        peaks: mz
            .iter()
            .zip(int)
            .map(|(&m, &i)| Peak1D::new(m, i as f32))
            .collect(),
        rt,
        ms_level: level,
        ..MSSpectrum::default()
    }
}

fn access(spectra: Vec<MSSpectrum>) -> SpectrumAccessPtr {
    SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(Arc::new(MSExperiment {
        spectra,
        ..MSExperiment::default()
    }))
    .unwrap()
}

fn scores_all(r: &mut Rec, k: &str, s: &OpenSwathScores) {
    for (name, v) in [
        ("massdev_score", s.massdev_score),
        ("weighted_massdev_score", s.weighted_massdev_score),
        ("dotprod_score_dia", s.dotprod_score_dia),
        ("manhatt_score_dia", s.manhatt_score_dia),
        ("isotope_correlation", s.isotope_correlation),
        ("isotope_overlap", s.isotope_overlap),
        ("bseries_score", s.bseries_score),
        ("yseries_score", s.yseries_score),
        ("ms1_ppm_score", s.ms1_ppm_score),
        ("ms1_isotope_correlation", s.ms1_isotope_correlation),
        ("ms1_isotope_overlap", s.ms1_isotope_overlap),
        ("xcorr_coelution_score", s.xcorr_coelution_score),
        ("weighted_coelution_score", s.weighted_coelution_score),
        ("xcorr_shape_score", s.xcorr_shape_score),
        ("weighted_xcorr_shape", s.weighted_xcorr_shape),
        ("ms1_xcorr_coelution_score", s.ms1_xcorr_coelution_score),
        ("ms1_xcorr_shape_score", s.ms1_xcorr_shape_score),
        (
            "ms1_xcorr_coelution_contrast_score",
            s.ms1_xcorr_coelution_contrast_score,
        ),
        (
            "ms1_xcorr_shape_contrast_score",
            s.ms1_xcorr_shape_contrast_score,
        ),
        (
            "ms1_xcorr_coelution_combined_score",
            s.ms1_xcorr_coelution_combined_score,
        ),
        (
            "ms1_xcorr_shape_combined_score",
            s.ms1_xcorr_shape_combined_score,
        ),
        ("nr_peaks", s.nr_peaks),
        ("sn_ratio", s.sn_ratio),
        ("log_sn_score", s.log_sn_score),
        ("mi_score", s.mi_score),
        ("weighted_mi_score", s.weighted_mi_score),
        ("ms1_mi_score", s.ms1_mi_score),
        ("ms1_mi_contrast_score", s.ms1_mi_contrast_score),
        ("ms1_mi_combined_score", s.ms1_mi_combined_score),
        ("library_corr", s.library_corr),
        ("library_norm_manhattan", s.library_norm_manhattan),
        ("norm_rt_score", s.norm_rt_score),
        ("raw_rt_score", s.raw_rt_score),
    ] {
        r.p(format!("{k}.{name}"), v);
    }
    im_scores(r, k, s);
}

/// driver.cpp `osw_case`.
#[test]
fn oracle_openswath_scoring() {
    let mut r = Rec::default();
    let half: Vec<f64> = SINT.iter().map(|v| v / 2.0).collect();
    let ms2 = access(vec![
        msspec(&SMZ, &half, 19.0, 2),
        msspec(&SMZ, &SINT, 20.0, 2),
        msspec(&SMZ, &half, 21.0, 2),
    ]);
    let ms1 = access(vec![msspec(
        &[499.99, 500.0, 500.01, 500.49, 500.5, 500.51, 501.0, 501.01],
        &[50., 100., 50., 30., 60., 30., 20., 5.],
        20.0,
        1,
    )]);
    let maps = [SwathMap {
        sptr: Some(ms2),
        ..SwathMap::default()
    }];
    let mut features = BTreeMap::new();
    features.insert(
        "group1".to_string(),
        mf(&[1., 3., 6., 3., 1., 0.5], 6.0, 20.),
    );
    features.insert(
        "group2".to_string(),
        mf(&[2., 5., 9., 4., 1., 0.2], 9.0, 20.),
    );
    let mut pre = BTreeMap::new();
    pre.insert("p0".to_string(), mf(&[0., 4., 8., 5., 2., 1.], 8.0, 20.));
    pre.insert("p1".to_string(), mf(&[1., 2., 5., 3., 1., 1.], 5.0, 20.));
    let f = MockMRMFeature {
        features,
        precursor_features: pre,
        intensity: 15.0,
        rt: 20.0,
        metavalue: 0.0,
    };
    let ts = [
        tr("group1", 500., 1, 3., 500.),
        tr("group2", 600., 1, 7., 500.),
    ];
    let pep = LightCompound {
        sequence: "PEPTIDEK".into(),
        charge: 2,
        rt: 25.0,
        ..LightCompound::default()
    };
    let met = LightCompound {
        sequence: "C6H12O6".into(),
        compound_name: "glucose".into(),
        charge: 1,
        ..LightCompound::default()
    };
    let d = dia(0.05);
    let su = OpenSwathScoresUsage::default();
    let nat = ids(&["group1", "group2"]);
    let prec = ids(&["p0", "p1"]);
    for method in ["simple", "resample"] {
        let mut sc = OpenSwathScoring::new();
        sc.initialize(
            100.0, 3, 0.005, 0.2, 0.25, &su, method, "fixed", false, false,
        )
        .unwrap();
        let mut s = OpenSwathScores::default();
        let ppm = sc
            .calculate_dia_scores(
                &f,
                &ts,
                &[0.3, 0.7],
                &maps,
                Some(&ms1),
                &d,
                &pep,
                &mut s,
                -1.0,
                &RangeBase::new(),
                None,
                -1,
            )
            .unwrap();
        sc.calculate_chromatographic_scores(&f, &nat, &prec, &[0.3, 0.7], &[], &mut s)
            .unwrap();
        sc.calculate_library_scores(&f, &ts, &pep, 20.0, &mut s)
            .unwrap();
        let k = format!("osw_{method}");
        scores_all(&mut r, &k, &s);
        r.pv(&format!("{k}.ppm"), &ppm);
        let mut s2 = OpenSwathScores::default();
        sc.calculate_precursor_dia_scores(
            Some(&ms1),
            &d,
            500.0,
            20.0,
            &met,
            RangeBase::new(),
            &mut s2,
        )
        .unwrap();
        r.p(format!("{k}.met.ms1_ppm_score"), s2.ms1_ppm_score);
        r.p(
            format!("{k}.met.ms1_isotope_correlation"),
            s2.ms1_isotope_correlation,
        );
        r.p(
            format!("{k}.met.ms1_isotope_overlap"),
            s2.ms1_isotope_overlap,
        );
        let fetched = sc
            .fetch_spectrum_swath_maps(&maps, 20.0, 3, &RangeBase::new())
            .unwrap();
        r.p(format!("{k}.fetch.n"), fetched.len() as f64);
        for (i, sp) in fetched.iter().enumerate() {
            r.p(
                format!("{k}.fetch.{i}.len"),
                sp.mz_array().data.len() as f64,
            );
            r.p(
                format!("{k}.fetch.{i}.sum"),
                sp.intensity_array().data.iter().fold(0.0, |a, v| a + v),
            );
        }
    }
    let mut sc = OpenSwathScoring::new();
    sc.initialize(
        100.0, 1, 0.005, 0.2, 0.25, &su, "simple", "fixed", false, false,
    )
    .unwrap();
    let mut ind = OpenSwathIndScores::default();
    sc.calculate_chromatographic_id_scores(&f, &ids(&["group1"]), &nat, &[], &mut ind)
        .unwrap();
    r.pv("ind.coel", &ind.ind_xcorr_coelution_score);
    r.pv("ind.shape", &ind.ind_xcorr_shape_score);
    r.pv("ind.mi", &ind.ind_mi_score);
    for prefix in ["osw_simple.", "osw_resample.", "ind."] {
        compare(&r, prefix);
    }
}

// ---------------------------------------------------------------------------
// Ion-mobilogram peak picking (group O integration)
// ---------------------------------------------------------------------------
//
// `../oracle/integ-O-impick/driver.cpp` ran `IonMobilityScoring::driftScoring`
// and `driftIdScoring` with `apply_im_peak_picking` true and false on five
// spectra (`inputs.tsv`) against the Linux x86-64 Release install on
// ibminode06, twice, byte-identical. The records are the expected values; the
// `.nopick` records of the same run are the control that the picking alone
// changes the result.

const IMPICK_INPUTS: &str = include_str!("data/openswath_scoring/impick_inputs.tsv");
const IMPICK_ORACLE: &str = include_str!("data/openswath_scoring/impick_oracle.stdout");

fn impick_spectra() -> BTreeMap<String, SpectrumPtr> {
    let mut out = BTreeMap::new();
    let mut lines = IMPICK_INPUTS.lines();
    while let Some(header) = lines.next() {
        let mut cols = header.split('\t');
        assert_eq!(cols.next(), Some("S"));
        let name = cols.next().unwrap().to_string();
        let n: usize = cols.next().unwrap().parse().unwrap();
        let (mut mz, mut im, mut int) = (Vec::new(), Vec::new(), Vec::new());
        for _ in 0..n {
            let row: Vec<f64> = lines
                .next()
                .unwrap()
                .split('\t')
                .map(|v| v.parse().unwrap())
                .collect();
            mz.push(row[0]);
            im.push(row[1]);
            int.push(row[2]);
        }
        out.insert(name, mk(&mz, &int, Some(&im)));
    }
    out
}

fn impick_scores(r: &mut Rec, k: &str, s: &OpenSwathScores) {
    for (name, v) in [
        ("im_drift", s.im_drift),
        ("im_drift_weighted", s.im_drift_weighted),
        ("im_delta_score", s.im_delta_score),
        ("im_delta", s.im_delta),
        ("im_log_intensity", s.im_log_intensity),
        ("im_drift_left", s.im_drift_left),
        ("im_drift_right", s.im_drift_right),
        ("im_xcorr_coelution_score", s.im_xcorr_coelution_score),
        ("im_xcorr_shape_score", s.im_xcorr_shape_score),
        ("im_ind_contrast_coelution", s.im_ind_contrast_coelution),
        ("im_ind_contrast_shape", s.im_ind_contrast_shape),
        (
            "im_ind_sum_contrast_coelution",
            s.im_ind_sum_contrast_coelution,
        ),
        ("im_ind_sum_contrast_shape", s.im_ind_sum_contrast_shape),
    ] {
        r.p(format!("{k}.{name}"), v);
    }
}

/// The driver's `fresh()`: the left and right drift start at the sentinel -2.
fn impick_fresh() -> OpenSwathScores {
    OpenSwathScores {
        im_drift_left: -2.0,
        im_drift_right: -2.0,
        ..OpenSwathScores::default()
    }
}

/// `driftScoring` and `driftIdScoring` with and without `apply_im_peak_picking`,
/// bit for bit against the Release build.
#[test]
fn oracle_im_peak_picking() {
    let spectra = impick_spectra();
    let det = [
        tr("det1", 500.2, 1, 1., 700.2),
        tr("det2", 600.5, 1, 1., 700.2),
    ];
    let id = [tr("id1", 700.1, 1, 1., 700.2)];
    let mut r = Rec::default();
    // (spectrum, drift target, range centre, range width, window, drift extra)
    let cases = [
        ("multi", 0.8, 1.0, 1.0, 0.3, 0.25),
        ("equal", 0.7, 1.0, 1.0, 0.3, 0.25),
        ("single", 0.8, 0.8, 0.5, 0.3, 0.25),
        ("flat", 1.0, 1.0, 1.0, 0.3, 0.25),
    ];
    for (name, target, centre, width, window, extra) in cases {
        let seq = vec![spectra[name].clone()];
        for apply in [false, true] {
            let tag = format!("{name}{}", if apply { ".pick" } else { ".nopick" });
            let mut s = impick_fresh();
            IonMobilityScoring::drift_scoring(
                &seq,
                &det,
                &mut s,
                target,
                rng(centre, width),
                &ext(window, false, extra),
                apply,
                None,
            )
            .unwrap();
            impick_scores(&mut r, &format!("{tag}.drift"), &s);
            let mut s = impick_fresh();
            IonMobilityScoring::drift_id_scoring(
                &seq,
                &id,
                &det,
                &mut s,
                target,
                rng(centre, width),
                &ext(window, false, extra),
                apply,
                None,
            )
            .unwrap();
            impick_scores(&mut r, &format!("{tag}.id"), &s);
        }
    }
    {
        let seq = vec![spectra["imspec"].clone()];
        let ts = [
            tr("group1", 500.2, 1, 1., 700.2),
            tr("group2", 600.5, 1, 1., 700.2),
        ];
        let one = [tr("id1", 500.3, 1, 1., 700.2)];
        for apply in [false, true] {
            let tag = format!("imspec{}", if apply { ".pick" } else { ".nopick" });
            let mut s = impick_fresh();
            IonMobilityScoring::drift_scoring(
                &seq,
                &ts,
                &mut s,
                1.0,
                rng(1.0, 1.0),
                &ext(0.3, false, 0.25),
                apply,
                None,
            )
            .unwrap();
            impick_scores(&mut r, &format!("{tag}.drift"), &s);
            let mut s = impick_fresh();
            IonMobilityScoring::drift_id_scoring(
                &seq,
                &one,
                &ts,
                &mut s,
                1.0,
                rng(1.0, 1.0),
                &ext(0.3, false, 0.25),
                apply,
                None,
            )
            .unwrap();
            impick_scores(&mut r, &format!("{tag}.id"), &s);
        }
    }

    let expected: BTreeMap<&str, &str> = IMPICK_ORACLE
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .collect();
    let mut failures = Vec::new();
    for (key, text) in &expected {
        let want = text.parse::<f64>().unwrap();
        match r.0.get(*key) {
            None => failures.push(format!("{key}: missing in port (oracle {text})")),
            Some(&got) => {
                let same = if want.is_nan() {
                    got.is_nan()
                } else {
                    got.to_bits() == want.to_bits() || got == want
                };
                if !same {
                    failures.push(format!("{key}: port {got:.17e} oracle {want:.17e}"));
                }
            }
        }
    }
    for key in r.0.keys() {
        if !expected.contains_key(key.as_str()) {
            failures.push(format!("{key}: not in oracle"));
        }
    }
    assert!(
        expected.len() > 100,
        "oracle holds only {} records",
        expected.len()
    );
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
    // The picking must have an effect on the multi-peak spectrum: the control
    // run leaves the borders at the sentinel, the picked run sets them.
    assert_eq!(r.0["multi.nopick.drift.im_drift_left"], -2.0);
    assert!(r.0["multi.pick.drift.im_drift_left"] > 0.0);
}

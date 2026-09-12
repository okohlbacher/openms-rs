// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `RangeUtils_test.cpp` (33 sections, tier 3 transcribed literals)
//! plus tier-4 coverage for the four predicates the class test never
//! exercises and for the native retain helpers and limits.

use openms::Error;
use openms::kernel::range_utils::{
    COLLISION_ENERGY_ACCESSION, COLLISION_ENERGY_KEY, HasActivationMethod, HasMetaValue,
    HasPrecursorCharge, HasScanMode, HasScanPolarity, InIntensityRange, InMzRange,
    InPrecursorMZRange, IsEmptySpectrum, IsInCollisionEnergyRange, IsInIsolationWindow,
    IsInIsolationWindowSizeRange, IsZoomSpectrum, PeakPredicate, RangeFilterLimits,
    SpectrumPredicate, collision_energy,
};
use openms::kernel::{DataArray, MSExperiment, MSSpectrum, Peak1D, Precursor};
use openms::metadata::{
    ActivationMethod, CVTerm, MetaInfo, MetaValue, MetaValueData, Polarity, ScanMode,
};

fn spectrum_at(rt: f64) -> MSSpectrum {
    let mut spectrum = MSSpectrum::new();
    spectrum.rt = rt;
    spectrum
}

fn spectrum_level(ms_level: u32) -> MSSpectrum {
    let mut spectrum = MSSpectrum::new();
    spectrum.ms_level = ms_level;
    spectrum
}

fn precursor_with_methods(methods: &[ActivationMethod]) -> Precursor {
    Precursor {
        activation_methods: methods.iter().copied().collect(),
        ..Precursor::default()
    }
}

fn float_value(value: f64) -> MetaValue {
    MetaValue::new(MetaValueData::Float(value)).unwrap()
}

// --- InRTRange: mapped to `spectra_in_rt_range` and a closure predicate ---

#[test]
fn in_rt_range_constructor() {
    // Section `InRTRange(double min, double max, bool reverse = false)`: (5, 10, false).
    let experiment = MSExperiment::new();
    assert!(
        experiment
            .spectra_in_rt_range(5.0, 10.0, 0)
            .unwrap()
            .is_empty()
    );
    let predicate = |s: &MSSpectrum| 5.0 <= s.rt && s.rt <= 10.0;
    assert!(!predicate.keep(&MSSpectrum::new()));
    assert!(matches!(
        experiment.spectra_in_rt_range(10.0, 5.0, 0),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn in_rt_range_destructor() {
    // Section `[EXTRA]~InRTRange()`: a predicate value is dropped without effect.
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(spectrum_at(7.5));
    let predicate = |s: &MSSpectrum| 5.0 <= s.rt && s.rt <= 10.0;
    assert_eq!(experiment.retain_spectra(&predicate).unwrap(), 0);
    assert_eq!(experiment.spectra.len(), 1);
}

#[test]
fn in_rt_range_operator() {
    // Section `bool operator()(const SpectrumType& s) const` for InRTRange.
    let r = |s: &MSSpectrum| 5.0 <= s.rt && s.rt <= 10.0;
    let r2 = |s: &MSSpectrum| !r(s);
    for (rt, inside) in [
        (4.9, false),
        (5.0, true),
        (7.5, true),
        (10.0, true),
        (10.1, false),
    ] {
        let s = spectrum_at(rt);
        assert_eq!(r.keep(&s), inside, "rt {rt}");
        assert_eq!(r2.keep(&s), !inside, "rt {rt}");
    }
    // The same closed boundaries through the existing sorted query.
    let mut experiment = MSExperiment::new();
    for rt in [4.9, 5.0, 7.5, 10.0, 10.1] {
        experiment.spectra.push(spectrum_at(rt));
    }
    let inside: Vec<f64> = experiment
        .spectra_in_rt_range(5.0, 10.0, 0)
        .unwrap()
        .iter()
        .map(|s| s.rt)
        .collect();
    assert_eq!(inside, [5.0, 7.5, 10.0]);
    // Source removal `erase_if(spectra, InRTRange(5, 10))` keeps the outside.
    assert_eq!(experiment.retain_spectra(&r2).unwrap(), 3);
    let left: Vec<f64> = experiment.spectra.iter().map(|s| s.rt).collect();
    assert_eq!(left, [4.9, 10.1]);
}

// --- InMSLevelRange: mapped to `ms_levels`/`contains_scan_of_level` and a closure ---

#[test]
fn in_ms_level_range_constructor() {
    // Section `MSLevelRange(const IntList& levels, bool reverse = false)`: empty list.
    let levels: Vec<u32> = Vec::new();
    let predicate = move |s: &MSSpectrum| levels.contains(&s.ms_level);
    assert!(!predicate.keep(&MSSpectrum::new()));
    let experiment = MSExperiment::new();
    assert!(!experiment.contains_scan_of_level(1).unwrap());
    assert!(experiment.ms_levels().is_empty());
}

#[test]
fn in_ms_level_range_destructor() {
    // Section `[EXTRA]~InMSLevelRange()`.
    let levels: Vec<u32> = Vec::new();
    let predicate = move |s: &MSSpectrum| levels.contains(&s.ms_level);
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(spectrum_level(2));
    assert_eq!(experiment.retain_spectra(&predicate).unwrap(), 1);
    assert!(experiment.spectra.is_empty());
}

#[test]
fn in_ms_level_range_operator() {
    // Section `bool operator()(const SpectrumType& s) const` for InMSLevelRange.
    let levels = [2_u32, 3, 4];
    let r = |s: &MSSpectrum| levels.contains(&s.ms_level);
    let r2 = |s: &MSSpectrum| !r(s);
    for (level, inside) in [(1, false), (2, true), (3, true), (4, true), (5, false)] {
        let s = spectrum_level(level);
        assert_eq!(r.keep(&s), inside, "level {level}");
        assert_eq!(r2.keep(&s), !inside, "level {level}");
    }
    let mut experiment = MSExperiment::new();
    for level in levels {
        experiment.spectra.push(spectrum_level(level));
    }
    assert_eq!(experiment.ms_levels(), vec![2, 3, 4]);
    for (level, inside) in [(1, false), (2, true), (3, true), (4, true), (5, false)] {
        assert_eq!(experiment.contains_scan_of_level(level).unwrap(), inside);
    }
}

// --- HasScanMode ---

#[test]
fn has_scan_mode_constructor() {
    // Section `HasScanMode(Int mode, bool reverse = false)`: mode 1 is MASSSPECTRUM.
    let predicate = HasScanMode::new(ScanMode::MassSpectrum, false);
    assert_eq!(predicate, HasScanMode::new(ScanMode::MassSpectrum, false));
}

#[test]
fn has_scan_mode_destructor() {
    // Section `[EXTRA]~HasScanMode()`.
    let _dropped_at_scope_end = HasScanMode::new(ScanMode::MassSpectrum, false);
}

#[test]
fn has_scan_mode_operator() {
    let r = HasScanMode::new(ScanMode::SelectedIonMonitoring, false);
    let r2 = HasScanMode::new(ScanMode::MassSpectrum, true);
    let mut s = MSSpectrum::new();
    s.instrument_settings.scan_mode = ScanMode::SelectedIonMonitoring;
    assert!(r.keep(&s));
    assert!(r2.keep(&s));
    s.instrument_settings.scan_mode = ScanMode::MassSpectrum;
    assert!(!r.keep(&s));
    assert!(!r2.keep(&s));
}

// --- InMzRange ---

#[test]
fn in_mz_range_constructor() {
    // Section `InMzRange(double min, double max, bool reverse = false)`: (5.0, 10.0, false).
    let predicate = InMzRange::new(5.0, 10.0, false).unwrap();
    assert_eq!(predicate, InMzRange::new(5.0, 10.0, false).unwrap());
}

#[test]
fn in_mz_range_destructor() {
    // Section `[EXTRA]~InMzRange()`.
    let _dropped_at_scope_end = InMzRange::new(5.0, 10.0, false).unwrap();
}

#[test]
fn in_mz_range_operator() {
    let r = InMzRange::new(5.0, 10.0, false).unwrap();
    let r2 = InMzRange::new(5.0, 10.0, true).unwrap();
    for (mz, inside) in [
        (4.9, false),
        (5.0, true),
        (7.5, true),
        (10.0, true),
        (10.1, false),
    ] {
        let p = Peak1D::new(mz, 0.0);
        assert_eq!(r.keep(&p), inside, "mz {mz}");
        assert_eq!(r2.keep(&p), !inside, "mz {mz}");
        assert_eq!(r.contains(mz), inside);
    }
}

// --- InIntensityRange ---

#[test]
fn in_intensity_range_constructor() {
    // Section `IntensityRange(double min, double max, bool reverse = false)`.
    let predicate = InIntensityRange::new(5.0, 10.0, false).unwrap();
    assert_eq!(predicate, InIntensityRange::new(5.0, 10.0, false).unwrap());
}

#[test]
fn in_intensity_range_destructor() {
    // Section `[EXTRA]~InIntensityRange()`.
    let _dropped_at_scope_end = InIntensityRange::new(5.0, 10.0, false).unwrap();
}

#[test]
fn in_intensity_range_operator() {
    let r = InIntensityRange::new(5.0, 10.0, false).unwrap();
    let r2 = InIntensityRange::new(5.0, 10.0, true).unwrap();
    // Source sets f32 intensities and widens them to double for the comparison.
    for (intensity, inside) in [
        (4.9_f32, false),
        (5.0, true),
        (7.5, true),
        (10.0, true),
        (10.1, false),
    ] {
        let p = Peak1D::new(0.0, intensity);
        assert_eq!(r.keep(&p), inside, "intensity {intensity}");
        assert_eq!(r2.keep(&p), !inside, "intensity {intensity}");
        assert_eq!(r.contains(intensity), inside);
    }
}

// --- IsEmptySpectrum ---

#[test]
fn is_empty_spectrum_constructor() {
    // Section `IsEmptySpectrum(bool reverse = false)`: default argument.
    assert_eq!(IsEmptySpectrum::new(false), IsEmptySpectrum::new(false));
}

#[test]
fn is_empty_spectrum_destructor() {
    // Section `[EXTRA]~IsEmptySpectrum()`.
    let _dropped_at_scope_end = IsEmptySpectrum::new(false);
}

#[test]
fn is_empty_spectrum_operator() {
    let s = IsEmptySpectrum::new(false);
    let s2 = IsEmptySpectrum::new(true);
    let mut spec = MSSpectrum::new();
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
    // Source `spec.resize(5)`: five default peaks.
    spec.peaks = vec![Peak1D::default(); 5];
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
}

// --- IsZoomSpectrum ---

#[test]
fn is_zoom_spectrum_constructor() {
    // Section `IsZoomSpectrum(bool reverse = false)`.
    assert_eq!(IsZoomSpectrum::new(false), IsZoomSpectrum::new(false));
}

#[test]
fn is_zoom_spectrum_destructor() {
    // Section `[EXTRA]~IsZoomSpectrum()`.
    let _dropped_at_scope_end = IsZoomSpectrum::new(false);
}

#[test]
fn is_zoom_spectrum_operator() {
    let s = IsZoomSpectrum::new(false);
    let s2 = IsZoomSpectrum::new(true);
    let mut spec = MSSpectrum::new();
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
    spec.instrument_settings.zoom_scan = true;
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

// --- HasActivationMethod ---

#[test]
fn has_activation_method_constructor() {
    // Section `HasActivationMethod(const StringList& methods, bool reverse = false)`:
    // the source builds it from `ListUtils::create<std::string>("")`, one empty
    // name that never matches. The port rejects an unknown name instead.
    assert!(matches!(
        HasActivationMethod::from_names([""], false),
        Err(Error::InvalidValue(_))
    ));
    let predicate = HasActivationMethod::new([], false);
    assert!(!predicate.keep(&MSSpectrum::new()));
}

#[test]
fn has_activation_method_destructor() {
    // Section `[EXTRA]~HasActivationMethod()`.
    let _dropped_at_scope_end = HasActivationMethod::new([], false);
}

#[test]
fn has_activation_method_operator() {
    // NamesOfActivationMethod[1] and [2]: "Post-source decay", "Plasma desorption".
    let names = [ActivationMethod::Psd.name(), ActivationMethod::Pd.name()];
    assert_eq!(names, ["Post-source decay", "Plasma desorption"]);
    let s = HasActivationMethod::from_names(names, false).unwrap();
    let s2 = HasActivationMethod::from_names(names, true).unwrap();
    assert_eq!(
        s,
        HasActivationMethod::new([ActivationMethod::Psd, ActivationMethod::Pd], false)
    );

    let mut spec = MSSpectrum::new();
    // PSD occurs; BIRD is just a dummy.
    spec.precursors = vec![precursor_with_methods(&[
        ActivationMethod::Psd,
        ActivationMethod::Bird,
    ])];
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));

    // Does not occur as activation method.
    spec.precursors[0] = precursor_with_methods(&[ActivationMethod::Bird]);
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    // Multiple precursors: adding another dummy.
    spec.precursors
        .push(precursor_with_methods(&[ActivationMethod::Lcid]));
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    // Adding a matching precursor.
    spec.precursors
        .push(precursor_with_methods(&[ActivationMethod::Pd]));
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

// --- InPrecursorMZRange ---

#[test]
fn in_precursor_mz_range_constructor() {
    // Section `InPrecursorMZRange(const double& mz_left, const double& mz_right, bool reverse)`.
    let predicate = InPrecursorMZRange::new(100.0, 200.0, false).unwrap();
    assert_eq!(
        predicate,
        InPrecursorMZRange::new(100.0, 200.0, false).unwrap()
    );
}

#[test]
fn in_precursor_mz_range_destructor() {
    // Section `[EXTRA]~InPrecursorMZRange()`.
    let _dropped_at_scope_end = InPrecursorMZRange::new(100.0, 200.0, false).unwrap();
}

#[test]
fn in_precursor_mz_range_operator() {
    let s = InPrecursorMZRange::new(100.0, 200.0, false).unwrap();
    let s2 = InPrecursorMZRange::new(100.0, 200.0, true).unwrap();
    let mut spec = MSSpectrum::new();
    spec.precursors = vec![Precursor::new(150.0, 0)];
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));

    // Outside of allowed window.
    spec.precursors[0] = Precursor::new(444.0, 0);
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    // Multiple precursors: the second is within limits, but all must be.
    spec.precursors.push(Precursor::new(150.0, 0));
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
}

// --- IsInIsolationWindow ---

#[test]
fn is_in_isolation_window_constructor() {
    // Section `IsInIsolationWindow(...)`: `ListUtils::create<double>("100.0, 200.0")`.
    let predicate = IsInIsolationWindow::new([100.0, 200.0], false).unwrap();
    assert_eq!(predicate.mz_values(), [100.0, 200.0]);
}

#[test]
fn is_in_isolation_window_destructor() {
    // Section `[EXTRA]~IsInIsolationWindow()`.
    let _dropped_at_scope_end = IsInIsolationWindow::new([100.0, 200.0], false).unwrap();
}

#[test]
fn is_in_isolation_window_operator() {
    // Unsorted on purpose, as in the source.
    let values = [300.0, 100.0, 200.0, 400.0];
    let s = IsInIsolationWindow::new(values, false).unwrap();
    let s2 = IsInIsolationWindow::new(values, true).unwrap();
    assert_eq!(s.mz_values(), [100.0, 200.0, 300.0, 400.0]);

    let mut spec = MSSpectrum::new();
    spec.ms_level = 2;
    let mut p = Precursor::new(200.3, 0);
    p.isolation_window_lower_offset = 0.5;
    p.isolation_window_upper_offset = 0.5;
    spec.precursors = vec![p.clone()];
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));

    // Outside of allowed window.
    p.mz = 201.1;
    spec.precursors[0] = p.clone();
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    // Multiple precursors: the second is within limits, so it's a hit (any PC).
    p.mz = 299.9;
    spec.precursors.push(p);
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

// --- HasScanPolarity ---

#[test]
fn has_scan_polarity_constructor() {
    // Section `HasScanPolarity(IonSource::Polarity polarity, bool reverse = false)`: POLNULL.
    let predicate = HasScanPolarity::new(Polarity::Unknown, false);
    assert_eq!(predicate, HasScanPolarity::new(Polarity::Unknown, false));
}

#[test]
fn has_scan_polarity_destructor() {
    // Section `[EXTRA]~HasScanPolarity()`.
    let _dropped_at_scope_end = HasScanPolarity::new(Polarity::Unknown, false);
}

#[test]
fn has_scan_polarity_operator() {
    let s = HasScanPolarity::new(Polarity::Positive, false);
    let s2 = HasScanPolarity::new(Polarity::Positive, true);
    let mut spec = MSSpectrum::new();
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
    spec.instrument_settings.polarity = Polarity::Positive;
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

// --- Tier 4: predicates without a class-test section ---

#[test]
fn has_precursor_charge_any_precursor_matches() {
    let s = HasPrecursorCharge::new([2, 3], false);
    let s2 = HasPrecursorCharge::new([2, 3], true);
    let mut spec = MSSpectrum::new();
    // No precursors: nothing matches.
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
    spec.precursors = vec![Precursor::new(500.0, 2)];
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
    spec.precursors[0].charge = 1;
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
    // Any listed charge across several precursors matches.
    spec.precursors.push(Precursor::new(600.0, 3));
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
    // Unknown charge zero can be listed explicitly.
    let zero = HasPrecursorCharge::new([0], false);
    assert!(zero.keep(&MSSpectrum {
        precursors: vec![Precursor::default()],
        ..MSSpectrum::new()
    }));
}

#[test]
fn is_in_collision_energy_range_reads_source_key_and_term() {
    let s = IsInCollisionEnergyRange::new(20.0, 40.0, false).unwrap();
    let s2 = IsInCollisionEnergyRange::new(20.0, 40.0, true).unwrap();

    // MS1 spectra: false regardless of reverse (source code, not source note).
    let mut ms1 = MSSpectrum::new();
    ms1.precursors = vec![Precursor::default()];
    ms1.precursors[0]
        .cv_terms
        .metadata
        .insert(COLLISION_ENERGY_KEY.into(), float_value(30.0));
    assert!(!s.keep(&ms1));
    assert!(!s2.keep(&ms1));

    // MS2 without any collision energy: false regardless of reverse.
    let mut spec = spectrum_level(2);
    spec.precursors = vec![Precursor::new(500.0, 2)];
    assert!(!s.keep(&spec));
    assert!(!s2.keep(&spec));
    assert_eq!(collision_energy(&spec.precursors[0]), None);

    // Source metadata key, closed boundaries.
    for (energy, inside) in [
        (19.9, false),
        (20.0, true),
        (30.0, true),
        (40.0, true),
        (40.1, false),
    ] {
        spec.precursors[0]
            .cv_terms
            .metadata
            .insert(COLLISION_ENERGY_KEY.into(), float_value(energy));
        assert_eq!(collision_energy(&spec.precursors[0]), Some(energy));
        assert_eq!(s.keep(&spec), inside, "energy {energy}");
        assert_eq!(s2.keep(&spec), !inside, "energy {energy}");
    }

    // Integer-typed values convert as the source DataValue does.
    spec.precursors[0]
        .cv_terms
        .metadata
        .insert(COLLISION_ENERGY_KEY.into(), MetaValue::from(25_i64));
    assert_eq!(collision_energy(&spec.precursors[0]), Some(25.0));
    assert!(s.keep(&spec));

    // A string value counts as absent here; the source conversion misbehaves.
    spec.precursors[0]
        .cv_terms
        .metadata
        .insert(COLLISION_ENERGY_KEY.into(), MetaValue::from("35"));
    assert_eq!(collision_energy(&spec.precursors[0]), None);
    assert!(!s.keep(&spec));
    assert!(!s2.keep(&spec));

    // Retained MS:1000045 term as the fallback.
    let mut precursor = Precursor::new(500.0, 2);
    let mut term = CVTerm::new(COLLISION_ENERGY_ACCESSION, "collision energy", "MS");
    term.value = float_value(50.0);
    precursor.cv_terms.add(term).unwrap();
    spec.precursors = vec![precursor];
    assert_eq!(collision_energy(&spec.precursors[0]), Some(50.0));
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    // Any precursor in range makes the spectrum match.
    let mut second = Precursor::new(600.0, 2);
    second
        .cv_terms
        .metadata
        .insert(COLLISION_ENERGY_KEY.into(), float_value(35.0));
    spec.precursors.push(second);
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

#[test]
fn is_in_isolation_window_size_range_sums_offsets() {
    let s = IsInIsolationWindowSizeRange::new(0.5, 2.0, false).unwrap();
    let s2 = IsInIsolationWindowSizeRange::new(0.5, 2.0, true).unwrap();

    // MS1: false regardless of reverse.
    let mut ms1 = MSSpectrum::new();
    let mut p = Precursor::new(500.0, 2);
    p.isolation_window_lower_offset = 0.5;
    p.isolation_window_upper_offset = 0.5;
    ms1.precursors = vec![p.clone()];
    assert!(!s.keep(&ms1));
    assert!(!s2.keep(&ms1));

    // MS2 without precursors: nothing is in range, so the answer is `reverse`.
    let mut spec = spectrum_level(2);
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));

    for (lower, upper, inside) in [
        (0.2, 0.2, false),
        (0.25, 0.25, true),
        (0.5, 0.5, true),
        (1.0, 1.0, true),
        (1.0, 1.5, false),
    ] {
        p.isolation_window_lower_offset = lower;
        p.isolation_window_upper_offset = upper;
        spec.precursors = vec![p.clone()];
        assert_eq!(s.keep(&spec), inside, "width {}", lower + upper);
        assert_eq!(s2.keep(&spec), !inside, "width {}", lower + upper);
    }

    // Any precursor in range matches.
    let mut narrow = Precursor::new(600.0, 2);
    narrow.isolation_window_lower_offset = 0.5;
    narrow.isolation_window_upper_offset = 0.5;
    spec.precursors.push(narrow);
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
}

#[test]
fn has_meta_value_checks_key_presence() {
    let s = HasMetaValue::new("label", false);
    let s2 = HasMetaValue::new("label", true);
    let mut spec = MSSpectrum::new();
    assert!(!s.keep(&spec));
    assert!(s2.keep(&spec));
    spec.metadata.insert("label".into(), MetaValue::from("x"));
    assert!(s.keep(&spec));
    assert!(!s2.keep(&spec));
    // An empty value still exists as a key.
    spec.metadata.insert("label".into(), MetaValue::default());
    assert!(s.keep(&spec));
    // Any metadata owner through `evaluate`.
    let mut meta = MetaInfo::new();
    assert!(!s.evaluate(&meta));
    meta.insert("label".into(), MetaValue::from(1_i64));
    assert!(s.evaluate(&meta));
}

// --- Tier 4: native construction checks and retain helpers ---

#[test]
fn inverted_and_nonfinite_ranges_are_rejected() {
    assert!(matches!(
        InMzRange::new(10.0, 5.0, false),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        InIntensityRange::new(10.0, 5.0, true),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        InPrecursorMZRange::new(200.0, 100.0, false),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        IsInCollisionEnergyRange::new(40.0, 20.0, false),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        IsInIsolationWindowSizeRange::new(2.0, 0.5, false),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        InMzRange::new(f64::NAN, 5.0, false),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        InIntensityRange::new(0.0, f64::INFINITY, false),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        IsInIsolationWindow::new([100.0, f64::NAN], false),
        Err(Error::InvalidValue(_))
    ));
    // A degenerate closed range is valid and matches exactly its point.
    let point = InMzRange::new(5.0, 5.0, false).unwrap();
    assert!(point.contains(5.0));
    assert!(!point.contains(5.000001));
}

#[test]
fn retain_spectra_counts_removals_and_leaves_chromatograms() {
    let mut experiment = MSExperiment::new();
    for level in [1, 2, 2, 1, 3] {
        experiment.spectra.push(spectrum_level(level));
    }
    experiment
        .chromatograms
        .push(openms::kernel::MSChromatogram::new());
    let removed = experiment
        .retain_spectra(&IsEmptySpectrum::new(true))
        .unwrap();
    assert_eq!(removed, 5);
    assert_eq!(experiment.chromatograms.len(), 1);

    let mut experiment = MSExperiment::new();
    for level in [1, 2, 2, 1, 3] {
        experiment.spectra.push(spectrum_level(level));
    }
    let removed = experiment
        .retain_spectra(&|s: &MSSpectrum| s.ms_level != 2)
        .unwrap();
    assert_eq!(removed, 2);
    assert_eq!(experiment.ms_levels(), vec![1, 3]);
}

#[test]
fn retain_spectra_limit_is_checked_before_mutation() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(spectrum_level(1));
    experiment.spectra.push(spectrum_level(2));
    let limits = RangeFilterLimits { max_items: 1 };
    let result = experiment.retain_spectra_with_limits(&IsEmptySpectrum::new(true), limits);
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(experiment.spectra.len(), 2);
    assert_eq!(
        experiment
            .retain_spectra_with_limits(
                &IsEmptySpectrum::new(false),
                RangeFilterLimits { max_items: 2 }
            )
            .unwrap(),
        0
    );
}

#[test]
fn retain_peaks_where_moves_aligned_arrays() {
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(100.0, 1.0),
        Peak1D::new(200.0, 20.0),
        Peak1D::new(300.0, 300.0),
    ]);
    spectrum
        .float_data_arrays
        .push(DataArray::new("width", vec![0.1, 0.2, 0.3]));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("charge", vec![1, 2, 3]));
    spectrum.string_data_arrays.push(DataArray::new(
        "label",
        vec!["a".into(), "b".into(), "c".into()],
    ));
    // Source `erase_if(spectrum, InIntensityRange(0.0, 50.0))` removes the low peaks.
    let removed = spectrum
        .retain_peaks_where(&InIntensityRange::new(0.0, 50.0, true).unwrap())
        .unwrap();
    assert_eq!(removed, 2);
    assert_eq!(spectrum.peaks, [Peak1D::new(300.0, 300.0)]);
    assert_eq!(spectrum.float_data_arrays[0].data, [0.3]);
    assert_eq!(spectrum.integer_data_arrays[0].data, [3]);
    assert_eq!(spectrum.string_data_arrays[0].data, ["c".to_string()]);

    // Closures work on peaks as well.
    let mut spectrum =
        MSSpectrum::from_peaks(vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)]);
    assert_eq!(
        spectrum
            .retain_peaks_where(&|p: &Peak1D| p.mz > 150.0)
            .unwrap(),
        1
    );
    assert_eq!(spectrum.peaks[0].mz, 200.0);
}

#[test]
fn retain_peaks_where_errors_leave_the_spectrum_unchanged() {
    let mut spectrum =
        MSSpectrum::from_peaks(vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)]);
    let limits = RangeFilterLimits { max_items: 1 };
    let result = spectrum
        .retain_peaks_where_with_limits(&InMzRange::new(150.0, 250.0, false).unwrap(), limits);
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(spectrum.peaks.len(), 2);

    // A misaligned nonempty data array is rejected by the shared retain path.
    spectrum
        .float_data_arrays
        .push(DataArray::new("width", vec![0.1]));
    let result = spectrum.retain_peaks_where(&InMzRange::new(150.0, 250.0, false).unwrap());
    assert!(matches!(result, Err(Error::InvalidValue(_))));
    assert_eq!(spectrum.peaks.len(), 2);
    assert_eq!(spectrum.float_data_arrays[0].data, [0.1]);
}

#[test]
fn ms_level_zero_is_not_treated_as_survey_scan() {
    // The source tests `getMSLevel() == 1` only; level 0 falls through.
    let mut spec = spectrum_level(0);
    let mut p = Precursor::new(500.0, 2);
    p.isolation_window_lower_offset = 1.0;
    p.isolation_window_upper_offset = 1.0;
    spec.precursors = vec![p];
    assert!(
        IsInIsolationWindowSizeRange::new(1.0, 3.0, false)
            .unwrap()
            .keep(&spec)
    );
    assert!(
        IsInIsolationWindow::new([500.5], false)
            .unwrap()
            .keep(&spec)
    );
}

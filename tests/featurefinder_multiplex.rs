// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of the multiplex data-structure class tests at
//! source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`): `MultiplexDeltaMasses_test.cpp`,
//! `MultiplexDeltaMassesGenerator_test.cpp`,
//! `MultiplexIsotopicPeakPattern_test.cpp`,
//! `MultiplexSatelliteCentroided_test.cpp`,
//! `MultiplexSatelliteProfile_test.cpp`, `MultiplexFilteredPeak_test.cpp` and
//! `MultiplexFilteredMSExperiment_test.cpp`.
//!
//! Evidence tier 3 (transcribed class-test literals) unless a test says it is
//! tier 4, an invariant derived from the source's documented formula. No
//! expected value was derived from this port's output.

use std::collections::BTreeMap;

use openms::analysis::featurefinder::multiplex::delta_masses_generator::{
    KNOWN_LABELS, LabelMassShifts, Labelling, label_master_list,
};
use openms::analysis::featurefinder::multiplex::{
    DeltaMass, LabelSet, MultiplexDeltaMasses, MultiplexDeltaMassesGenerator,
    MultiplexFilteredMSExperiment, MultiplexFilteredPeak, MultiplexIsotopicPeakPattern,
    MultiplexSatelliteCentroided, MultiplexSatelliteProfile,
};
use openms::constants::C13C12_MASSDIFF_U;

/// `TEST_REAL_SIMILAR` at the ClassTest default relative tolerance 1e-5.
fn similar(actual: f64, expected: f64) {
    let tolerance = 1e-5 * expected.abs().max(actual.abs());
    assert!(
        (actual - expected).abs() <= tolerance || actual == expected,
        "{actual} is not similar to {expected}"
    );
}

fn arg6_pair() -> MultiplexDeltaMasses {
    let mut pattern = MultiplexDeltaMasses::new();
    pattern
        .delta_masses
        .push(DeltaMass::with_label(0.0, "no_label"));
    pattern
        .delta_masses
        .push(DeltaMass::with_label(6.031817, "Arg6"));
    pattern
}

// ------------------------------------------------------------ MultiplexDeltaMasses

#[test]
fn delta_masses_class_test() {
    // MultiplexDeltaMasses()
    assert_eq!(MultiplexDeltaMasses::new().delta_masses().len(), 0);
    // std::vector<DeltaMass>& getDeltaMasses()
    let pattern = arg6_pair();
    similar(pattern.delta_masses()[0].delta_mass, 0.0);
    similar(pattern.delta_masses()[1].delta_mass, 6.031817);
}

/// Tier 4: `labelSetToString` joins the multiset in ascending order with
/// single spaces; the multiset keeps repeated labels.
#[test]
fn label_set_is_an_ordered_multiset() {
    let set: LabelSet = ["Lys8", "Arg10", "Lys8"].into_iter().collect();
    assert_eq!(set.len(), 3);
    assert_eq!(set.count("Lys8"), 2);
    assert_eq!(
        MultiplexDeltaMasses::label_set_to_string(&set),
        "Arg10 Lys8 Lys8"
    );
    assert_eq!(
        MultiplexDeltaMasses::label_set_to_string(&LabelSet::new()),
        ""
    );
}

/// Tier 4: the source `operator<` puts longer patterns first, then compares the
/// shifts relative to the first shift.
#[test]
fn delta_masses_source_order() {
    let pair = arg6_pair();
    let mut triple = arg6_pair();
    triple.delta_masses.push(DeltaMass::with_label(12.0, "x"));
    assert!(triple.source_less(&pair));
    assert!(!pair.source_less(&triple));
    let mut larger = arg6_pair();
    larger.delta_masses[1].delta_mass = 7.0;
    assert!(pair.source_less(&larger));
    assert!(!larger.source_less(&pair));
    // Only relative shifts count, and labels play no part.
    let mut shifted = MultiplexDeltaMasses::new();
    shifted.delta_masses.push(DeltaMass::with_label(1.0, "a"));
    shifted
        .delta_masses
        .push(DeltaMass::with_label(7.031817, "b"));
    assert!(!shifted.source_less(&pair));
    assert!(!pair.source_less(&shifted));
}

// ------------------------------------------------------------ MultiplexIsotopicPeakPattern

#[test]
fn isotopic_peak_pattern_class_test() {
    let pattern = MultiplexIsotopicPeakPattern::new(2, 4, arg6_pair(), 3).unwrap();
    assert_eq!(pattern.charge(), 2);
    assert_eq!(pattern.peaks_per_peptide(), 4);
    assert_eq!(pattern.mass_shifts().delta_masses()[0].delta_mass, 0.0);
    assert_eq!(pattern.mass_shifts().delta_masses()[1].delta_mass, 6.031817);
    assert_eq!(pattern.mass_shift_index(), 3);
    assert_eq!(pattern.mass_shift_count(), 2);
    assert_eq!(pattern.mass_shift_at(0), Some(0.0));
    assert_eq!(pattern.mass_shift_at(1), Some(6.031817));
}

/// Tier 4. The class test's `getMZShiftAt` and `getMZShiftCount` sections are
/// commented out upstream (their literals describe an older pattern with a
/// zeroth peak), so the expectations here come from the constructor's formula
/// `(delta_mass_i + j * C13C12_MASSDIFF_U) / charge`.
#[test]
fn isotopic_peak_pattern_mz_shifts_follow_the_formula() {
    let pattern = MultiplexIsotopicPeakPattern::new(2, 4, arg6_pair(), 3).unwrap();
    assert_eq!(pattern.mz_shift_count(), 8);
    for i in 0..2 {
        let delta = [0.0, 6.031817][i];
        for j in 0..4 {
            assert_eq!(
                pattern.mz_shift_at(i * 4 + j),
                Some((delta + j as f64 * C13C12_MASSDIFF_U) / 2.0)
            );
        }
    }
    assert_eq!(pattern.mz_shift_at(8), None);
    assert_eq!(pattern.mass_shift_at(2), None);
    // A zero charge is refused where the source divides by zero.
    assert!(MultiplexIsotopicPeakPattern::new(0, 4, arg6_pair(), 0).is_err());
    // A negative peak count yields no shifts, as the source loop does.
    assert_eq!(
        MultiplexIsotopicPeakPattern::new(2, -1, arg6_pair(), 0)
            .unwrap()
            .mz_shift_count(),
        0
    );
}

// ------------------------------------------------------------ satellites

#[test]
fn satellite_centroided_class_test() {
    let satellite = MultiplexSatelliteCentroided::new(4, 7);
    assert_eq!(satellite.mz_idx, 7);
    assert_eq!(satellite.rt_idx, 4);
}

#[test]
fn satellite_profile_class_test() {
    let satellite = MultiplexSatelliteProfile::new(2565.3, 618.4, 1000.0);
    assert_eq!(satellite.mz, 618.4);
    similar(satellite.mz, 618.4);
    similar(f64::from(satellite.rt), 2565.3);
    similar(f64::from(satellite.intensity), 1000.0);
}

// ------------------------------------------------------------ MultiplexFilteredPeak

#[test]
fn filtered_peak_class_test() {
    assert_eq!(
        MultiplexFilteredPeak::new(654.32, 2345.67, 24, 42).mz_idx(),
        24
    );

    let mut peak = MultiplexFilteredPeak::new(654.32, 2345.67, 24, 42);
    peak.add_satellite(25, 43, 3);
    peak.add_satellite_from(MultiplexSatelliteCentroided::new(26, 44), 3);
    peak.add_satellite_profile(2347.67, 656.32, 1010.0, 4);
    peak.add_satellite_profile_from(MultiplexSatelliteProfile::new(2346.67, 655.32, 1000.0), 4);

    similar(peak.mz(), 654.32);
    similar(f64::from(peak.rt()), 2345.67);
    assert_eq!(peak.mz_idx(), 24);
    assert_eq!(peak.rt_idx(), 42);

    // addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx)
    let n = peak.satellite_count();
    peak.add_satellite(25, 43, 3);
    assert_eq!(peak.satellite_count(), n + 1);
    // addSatellite(const MultiplexSatelliteCentroided&, size_t)
    let n = peak.satellite_count();
    peak.add_satellite_from(MultiplexSatelliteCentroided::new(27, 45), 3);
    assert_eq!(peak.satellite_count(), n + 1);
    // addSatelliteProfile(double rt, double mz, double intensity, size_t)
    let n = peak.satellite_profile_count();
    peak.add_satellite_profile(2348.67, 657.32, 1020.0, 5);
    assert_eq!(peak.satellite_profile_count(), n + 1);
    // addSatelliteProfile(const MultiplexSatelliteProfile&, size_t)
    let n = peak.satellite_profile_count();
    peak.add_satellite_profile_from(MultiplexSatelliteProfile::new(2349.67, 658.32, 1030.0), 6);
    assert_eq!(peak.satellite_profile_count(), n + 1);

    // getSatellites, getSatellitesProfile, size, sizeProfile
    assert_eq!(peak.satellite_entries().count(), 4);
    assert_eq!(peak.satellite_profile_entries().count(), 4);
    assert_eq!(peak.satellite_count(), 4);
    assert_eq!(peak.satellite_profile_count(), 4);
}

/// Tier 4: `checkSatellite` searches every mass trace; the multimap iterates
/// keys in ascending order and equal keys in insertion order.
#[test]
fn filtered_peak_satellites_behave_as_a_multimap() {
    let mut peak = MultiplexFilteredPeak::new(1.0, 1.0, 0, 0);
    peak.add_satellite(5, 1, 7);
    peak.add_satellite(3, 2, 2);
    peak.add_satellite(4, 3, 7);
    assert!(peak.check_satellite(4, 3));
    assert!(!peak.check_satellite(3, 1));
    let order: Vec<_> = peak
        .satellite_entries()
        .map(|(k, s)| (k, s.rt_idx, s.mz_idx))
        .collect();
    assert_eq!(order, vec![(2, 3, 2), (7, 5, 1), (7, 4, 3)]);
}

// ------------------------------------------------------------ MultiplexFilteredMSExperiment

#[test]
fn filtered_ms_experiment_class_test() {
    assert_eq!(MultiplexFilteredMSExperiment::new().len(), 0);
    let mut exp = MultiplexFilteredMSExperiment::new();
    exp.add_peak(MultiplexFilteredPeak::new(654.32, 2345.67, 24, 110));
    // addPeak
    let n = exp.len();
    exp.add_peak(MultiplexFilteredPeak::new(655.32, 2346.67, 25, 111));
    assert_eq!(exp.len(), n + 1);
    // getPeak, getMZ(i), getMZ(), getRT(i), getRT(), size
    similar(exp.peak(0).unwrap().mz(), 654.32);
    similar(exp.mz_at(0).unwrap(), 654.32);
    similar(exp.mz()[0], 654.32);
    similar(exp.rt_at(0).unwrap(), 2345.67);
    similar(exp.rt()[0], 2345.67);
    assert_eq!(exp.len(), 2);
    assert!(exp.peak(2).is_none());
}

// ------------------------------------------------------------ MultiplexDeltaMassesGenerator

fn triple_silac_masses() -> BTreeMap<String, f64> {
    [
        ("Arg6", 6.0201290268),
        ("Arg10", 10.008268600),
        ("Lys4", 4.0251069836),
        ("Lys8", 8.0141988132),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect()
}

#[test]
fn delta_masses_generator_class_test() {
    let labels = "[][Lys4,Arg6][Lys8,Arg10]";
    let mut list = MultiplexDeltaMassesGenerator::new(labels, 1, triple_silac_masses()).unwrap();
    assert_eq!(list.labelling(), Labelling::Silac);

    // getDeltaMassesList
    let masses = list.delta_masses_list();
    assert_eq!(masses.len(), 5);
    similar(masses[2].delta_masses()[1].delta_mass, 8.0502139672);
    similar(masses[4].delta_masses()[2].delta_mass, 20.0165372);

    // getSamplesLabelsList
    let samples = list.samples_labels_list();
    assert_eq!(samples.len(), 3);
    assert_eq!(samples[1][0], "Lys4");
    assert_eq!(samples[2][1], "Arg10");

    // generateKnockoutDeltaMasses
    list.generate_knockout_delta_masses().unwrap();
    let knockout = list.delta_masses_list();
    assert_eq!(knockout.len(), 21);
    similar(knockout[6].delta_masses()[1].delta_mass, 8.0141988132);
    similar(knockout[19].delta_masses()[1].delta_mass, 20.0165372);
    assert_eq!(knockout[20].delta_masses().len(), 1);

    // getLabelShort, getLabelLong
    assert_eq!(list.label_short("Label:13C(6)15N(2)"), Some("Lys8"));
    assert_eq!(list.label_short("Dimethyl:2H(6)13C(2)"), Some("Dimethyl8"));
    assert_eq!(list.label_long("Lys8"), Some("Label:13C(6)15N(2)"));
    assert_eq!(list.label_long("Dimethyl8"), Some("Dimethyl:2H(6)13C(2)"));

    // extractLabelSet
    let sequence =
        openms::chemistry::AASequence::parse("LAPITSDPTEAAAVGAVEASFK(Label:13C(6)15N(2))").unwrap();
    let label_set = list.extract_label_set(&sequence);
    assert_eq!(label_set.len(), 1);
    assert_eq!(label_set.iter().next(), Some("Lys8"));
}

/// Tier 4, from the header's documentation of `printDeltaMassesList`: a
/// standard SILAC experiment allowing one missed cleavage yields the five
/// patterns 8.0142 (Lys8), 10.0083 (Arg10), 16.0284 (Lys8,Lys8),
/// 18.0225 (Arg10,Lys8) and 20.0165 (Arg10,Arg10).
#[test]
fn delta_masses_generator_documented_silac_list() {
    let masses = LabelMassShifts::default().masses;
    let generator = MultiplexDeltaMassesGenerator::new("[][Lys8,Arg10]", 1, masses).unwrap();
    let text = generator.delta_masses_text();
    assert_eq!(
        text,
        "\nmass shift 1:    0 (no_label)    8.0142 (Lys8)    \n\
         mass shift 2:    0 (no_label)    10.0083 (Arg10)    \n\
         mass shift 3:    0 (no_label)    16.0284 (Lys8,Lys8)    \n\
         mass shift 4:    0 (no_label)    18.0225 (Arg10,Lys8)    \n\
         mass shift 5:    0 (no_label)    20.0165 (Arg10,Arg10)    \n\n"
    );
    assert_eq!(
        generator.samples_labels_text(),
        "\nsample 1:    no_label    \nsample 2:    Lys8    Arg10    \n"
    );
    let mut out = Vec::new();
    generator.print_delta_masses_list(&mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), text);
}

/// Tier 4: the documented labelling kinds and the source's parsing rules.
#[test]
fn delta_masses_generator_labelling_kinds() {
    let masses = LabelMassShifts::default().masses;
    // Dimethyl: each sample's first label against the first sample's.
    let g =
        MultiplexDeltaMassesGenerator::new("[Dimethyl0][Dimethyl8]", 1, masses.clone()).unwrap();
    assert_eq!(g.labelling(), Labelling::Dimethyl);
    assert_eq!(g.delta_masses_list().len(), 2);
    let first = g.delta_masses_list()[0].delta_masses();
    assert_eq!(first[0].delta_mass, 0.0);
    assert_eq!(first[1].delta_mass, 36.075670 - 28.031300);
    assert_eq!(first[1].label_set.as_slice(), ["Dimethyl8"]);
    let second = g.delta_masses_list()[1].delta_masses();
    assert_eq!(second[1].delta_mass, 2.0 * (36.075670 - 28.031300));
    assert_eq!(second[1].label_set.count("Dimethyl8"), 2);

    // No labelling: one zero shift.
    for labels in ["", "[]", "()", "{}"] {
        let g = MultiplexDeltaMassesGenerator::new(labels, 1, masses.clone()).unwrap();
        assert_eq!(g.labelling(), Labelling::None);
        assert_eq!(g.delta_masses_list().len(), 1);
        assert_eq!(g.samples_labels_list(), [vec!["no_label".to_owned()]]);
    }

    // Numeric labels, including a numeric prefix as std::stod reads it.
    let g = MultiplexDeltaMassesGenerator::new("[0][3.5abc]", 0, BTreeMap::new()).unwrap();
    assert_eq!(g.labelling(), Labelling::Numeric);
    assert_eq!(g.delta_masses_list()[0].delta_masses()[1].delta_mass, 3.5);
    assert_eq!(g.label_long("3.5abc"), Some("3.5abc"));

    // Leu with a single sample reads past the samples in the source.
    assert!(MultiplexDeltaMassesGenerator::new("[Leu3]", 1, masses.clone()).is_err());
    let g = MultiplexDeltaMassesGenerator::new("[][Leu3]", 1, masses.clone()).unwrap();
    assert_eq!(g.labelling(), Labelling::Leu);
    assert_eq!(
        g.delta_masses_list()[0].delta_masses()[1].delta_mass,
        3.018830
    );

    // Mixed labelling kinds and unknown labels are refused.
    assert!(MultiplexDeltaMassesGenerator::new("[Dimethyl0][Lys8]", 1, masses.clone()).is_err());
    assert!(MultiplexDeltaMassesGenerator::new("[][Lys9]", 1, masses.clone()).is_err());
    // The known-label test is a substring test, as in the source.
    assert!(KNOWN_LABELS.contains("Arg"));
    assert!(MultiplexDeltaMassesGenerator::new("[][Arg]", 1, masses.clone()).is_ok());
    // Negative counts wrap in the source; they are refused here.
    assert!(MultiplexDeltaMassesGenerator::new("[][Lys8]", -1, masses).is_err());
}

/// Tier 4: `Label:13C(6)` is the long name of Arg6 and Lys6, and the source map
/// keeps the first insertion.
#[test]
fn delta_masses_generator_shared_long_name_maps_to_arg6() {
    let g = MultiplexDeltaMassesGenerator::new("[][Lys6]", 0, LabelMassShifts::default().masses)
        .unwrap();
    assert_eq!(g.label_short("Label:13C(6)"), Some("Arg6"));
    assert_eq!(g.label_short("unknown"), None);
}

#[test]
fn delta_masses_generator_knockout_errors() {
    let masses = LabelMassShifts::default().masses;
    let mut singlet = MultiplexDeltaMassesGenerator::new("", 0, masses.clone()).unwrap();
    assert!(singlet.generate_knockout_delta_masses().is_err());
    assert_eq!(singlet.delta_masses_list().len(), 1);
    let mut pair = MultiplexDeltaMassesGenerator::new("[][Lys8]", 0, masses).unwrap();
    pair.generate_knockout_delta_masses().unwrap();
    let list = pair.delta_masses_list();
    assert_eq!(list.len(), 2);
    assert_eq!(
        list[1].delta_masses()[0].label_set.as_slice(),
        ["any_label_set"]
    );
}

/// The `labels` parameter section is one float parameter per master-list
/// label, restricted to `0:` (source default constructor).
#[test]
fn label_mass_shift_parameters() {
    let defaults = LabelMassShifts::defaults().unwrap();
    for label in label_master_list() {
        assert_eq!(
            defaults.value(&label.short_name).unwrap().to_f64().unwrap(),
            label.delta_mass
        );
    }
    let parsed = LabelMassShifts::from_param(&defaults).unwrap();
    assert_eq!(parsed, LabelMassShifts::default());
    assert_eq!(parsed.masses.len(), 14);
}

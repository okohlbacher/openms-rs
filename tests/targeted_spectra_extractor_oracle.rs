// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (tier 1) for `TargetedSpectraExtractor` on the class
//! test's data.
//!
//! `tests/data/targeted_spectra_extractor/oracle_expected.tsv` is the output of
//! `../oracle/O4b-tse/driver.cpp` linked against the Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
//! ibminode06 (byte-identical runs). The driver loads the target list with the
//! source's `TransitionTSVFile`; this test builds it from the same CSV and
//! first checks that the two agree on every field the extractor reads.
//!
//! Names, counts, annotation coordinates and selection membership compare
//! exactly. Values computed through smoothing, picking and noise estimation
//! (picked m/z, intensity and FWHM, the scores and their parts, library match
//! scores) compare within `1e-9` relative (absolute below `1e-9`), a bound
//! fixed before the first run; the number of bit-identical values is printed.

#![cfg(feature = "mzml")]

use openms::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::analysis::openswath::targeted_experiment_helper::{Peptide, RTUnit, RetentionTime};
use openms::analysis::openswath::targeted_spectra_extractor::{
    BinnedSpectrumComparator, Comparator, TargetedSpectraExtractor,
};
use openms::kernel::features::FeatureMap;
use openms::kernel::{MSExperiment, MSSpectrum};
use openms::metadata::MetaValue;
use openms::param::ParamValue;
use std::collections::BTreeMap;

const EXPECTED: &str = include_str!("data/targeted_spectra_extractor/oracle_expected.tsv");

fn experiment() -> MSExperiment {
    openms::format::mzml::read(
        include_bytes!(
            "data/targeted_spectra_extractor/TargetedSpectraExtractor_13C1_spectra0to100.mzML"
        )
        .as_slice(),
    )
    .unwrap()
}

/// The target list, built as in `tests/targeted_spectra_extractor.rs`.
fn targeted_experiment() -> TargetedExperiment {
    let text =
        include_str!("data/targeted_spectra_extractor/TargetedSpectraExtractor_13CFlux_TraML.csv");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split(',').collect();
    let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (group, name, rt, precursor, product) = (
        column("transition_group_id"),
        column("transition_name"),
        column("RetentionTime"),
        column("PrecursorMz"),
        column("ProductMz"),
    );
    let mut exp = TargetedExperiment::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let fields: Vec<&str> = line.split(',').collect();
        let mut transition = ReactionMonitoringTransition::new();
        transition.set_name(fields[name]);
        transition.set_peptide_ref(fields[group]);
        transition.set_precursor_mz(fields[precursor].parse().unwrap());
        transition.set_product_mz(fields[product].parse().unwrap());
        exp.add_transition(transition);
        if !exp.has_peptide(fields[group]) {
            let mut peptide = Peptide::default();
            peptide.id = fields[group].to_owned();
            let mut retention_time = RetentionTime::default();
            retention_time.set_rt(fields[rt].parse().unwrap());
            retention_time.retention_time_unit = RTUnit::Minute;
            peptide.rts.push(retention_time);
            exp.add_peptide(peptide);
        }
    }
    exp
}

fn expected() -> BTreeMap<String, Vec<(String, String)>> {
    let mut map: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    for line in EXPECTED.lines().filter(|l| !l.is_empty()) {
        let mut f = line.splitn(3, '\t');
        let case = f.next().unwrap().to_owned();
        let key = f.next().unwrap().to_owned();
        let value = f.next().unwrap_or("").to_owned();
        map.entry(case).or_default().push((key, value));
    }
    map
}

fn bits(v: f64) -> String {
    format!("{:016x}", v.to_bits())
}

fn from_bits(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

/// Compares records: exact unless the key is in a floating field list, where
/// comma- or bar-separated hex values compare within 1e-9 relative.
struct Comparison {
    compared: usize,
    identical_floats: usize,
    floats: usize,
    worst: f64,
}

impl Comparison {
    fn new() -> Self {
        Self {
            compared: 0,
            identical_floats: 0,
            floats: 0,
            worst: 0.0,
        }
    }

    fn close(&mut self, case: &str, key: &str, got: f64, want: f64) {
        self.floats += 1;
        if got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()) {
            self.identical_floats += 1;
            return;
        }
        let diff = (got - want).abs();
        let relative = if want.abs() > 1e-9 {
            diff / want.abs()
        } else {
            diff
        };
        self.worst = self.worst.max(relative);
        assert!(
            relative <= 1e-9,
            "{case} {key}: {got:e} vs {want:e} (relative {relative:e})"
        );
    }

    /// `got` and `want` are `|`/`,`-separated fields; `float_fields` marks the
    /// positions compared as floats, the others exactly.
    fn fields(&mut self, case: &str, key: &str, got: &str, want: &str, float_fields: &[usize]) {
        let split = |s: &str| s.split(['|', ',']).map(str::to_owned).collect::<Vec<_>>();
        let (g, w) = (split(got), split(want));
        assert_eq!(g.len(), w.len(), "{case} {key}: {got} vs {want}");
        for (i, (a, b)) in g.iter().zip(&w).enumerate() {
            if float_fields.contains(&i) && a != "-" && b != "-" {
                self.close(case, key, from_bits(a), from_bits(b));
            } else {
                assert_eq!(a, b, "{case} {key} field {i}");
            }
        }
        self.compared += 1;
    }

    fn exact(&mut self, case: &str, key: &str, got: &str, want: &str) {
        assert_eq!(got, want, "{case} {key}");
        self.compared += 1;
    }
}

fn set(values: &[(&str, ParamValue)]) -> TargetedSpectraExtractor {
    let mut params = TargetedSpectraExtractor::new().unwrap().to_param().unwrap();
    for (key, value) in values {
        let description = params.description(key).unwrap_or("").to_owned();
        params
            .set_value(key, value.clone(), &description, &[])
            .unwrap();
    }
    TargetedSpectraExtractor::from_param(&params).unwrap()
}

fn class_params() -> Vec<(&'static str, ParamValue)> {
    vec![
        ("GaussFilter:gaussian_width", ParamValue::Float(0.25)),
        ("peak_height_min", ParamValue::Float(15000.0)),
        ("peak_height_max", ParamValue::Float(110000.0)),
        ("fwhm_threshold", ParamValue::Float(0.23)),
    ]
}

fn picked_records(prefix: &str, picked: &MSSpectrum) -> Vec<(String, String)> {
    let mut out = vec![
        (format!("{prefix}.count"), picked.len().to_string()),
        (
            format!("{prefix}.arrays"),
            picked.float_data_arrays.len().to_string(),
        ),
    ];
    for (j, peak) in picked.peaks.iter().enumerate() {
        let fwhm = picked
            .float_data_arrays
            .first()
            .map_or("-".to_owned(), |a| bits(f64::from(a.data[j])));
        out.push((
            format!("{prefix}.{j}"),
            format!(
                "{},{},{}",
                bits(peak.mz),
                bits(f64::from(peak.intensity)),
                fwhm
            ),
        ));
    }
    out
}

#[test]
fn every_oracle_record_matches() {
    let expected = expected();
    let experiment = experiment();
    let targets = targeted_experiment();
    let mut cmp = Comparison::new();
    let check = |cmp: &mut Comparison,
                 case: &str,
                 got: Vec<(String, String)>,
                 floats: &dyn Fn(&str) -> Vec<usize>| {
        let want = &expected[case];
        assert_eq!(
            got.iter().map(|r| &r.0).collect::<Vec<_>>(),
            want.iter().map(|r| &r.0).collect::<Vec<_>>(),
            "{case}: record keys"
        );
        for ((key, value), (_, want_value)) in got.iter().zip(want) {
            let positions = floats(key);
            if positions.is_empty() {
                cmp.exact(case, key, value, want_value);
            } else {
                cmp.fields(case, key, value, want_value, &positions);
            }
        }
    };

    // The target list agrees with TransitionTSVFile's.
    let mut got = vec![(
        "transitions".to_owned(),
        targets.transitions().len().to_string(),
    )];
    for (i, t) in targets.transitions().iter().enumerate() {
        got.push((
            format!("transition.{i}"),
            format!(
                "{}|{}|{}",
                t.name(),
                t.peptide_ref(),
                bits(t.precursor_mz())
            ),
        ));
    }
    got.push(("peptides".to_owned(), targets.peptides().len().to_string()));
    for p in targets.peptides() {
        got.push((
            format!("peptide.{}", p.id),
            format!(
                "{}|{}",
                bits(p.retention_time().unwrap()),
                p.retention_time_unit().unwrap() as i8
            ),
        ));
    }
    check(&mut cmp, "targets", got, &|_| vec![]);

    let tse = set(&class_params());
    let (annotated, features) = tse.annotate_spectra(&experiment.spectra, &targets).unwrap();
    let mut got = vec![("count".to_owned(), annotated.len().to_string())];
    for (i, (s, f)) in annotated.iter().zip(&features.features).enumerate() {
        got.push((
            i.to_string(),
            format!("{}|{}|{}|{}", s.name, s.len(), bits(f.rt), bits(f.mz)),
        ));
    }
    check(&mut cmp, "annotate", got, &|_| vec![]);

    let pick_float = |key: &str| {
        if key.ends_with(".count") || key.ends_with(".arrays") {
            vec![]
        } else {
            vec![0, 1, 2]
        }
    };
    for (case, values) in [
        ("pick_class", class_params()),
        (
            "pick_nofilter",
            vec![
                ("GaussFilter:gaussian_width", ParamValue::Float(0.25)),
                ("peak_height_max", ParamValue::Float(200000.0)),
            ],
        ),
        (
            "pick_sgolay",
            vec![
                ("use_gauss", ParamValue::String("false".into())),
                ("PeakPickerHiRes:signal_to_noise", ParamValue::Float(0.5)),
            ],
        ),
    ] {
        let tse = set(&values);
        let annotated = tse
            .annotate_spectra_only(&experiment.spectra, &targets)
            .unwrap();
        let mut got = Vec::new();
        for (i, spectrum) in annotated.iter().enumerate() {
            got.extend(picked_records(
                &i.to_string(),
                &tse.pick_spectrum(spectrum).unwrap(),
            ));
        }
        check(&mut cmp, case, got, &pick_float);
    }

    // Scoring and selection.
    let tse = set(&class_params());
    let (annotated, mut features) = tse.annotate_spectra(&experiment.spectra, &targets).unwrap();
    let mut kept = (Vec::new(), Vec::new(), FeatureMap::new());
    for (i, spectrum) in annotated.iter().enumerate() {
        let picked = tse.pick_spectrum(spectrum).unwrap();
        if !picked.is_empty() {
            kept.0.push(spectrum.clone());
            kept.1.push(picked);
            kept.2.features.push(features.features[i].clone());
        }
    }
    features = kept.2;
    let scored = tse.score_spectra(&kept.0, &kept.1, &mut features).unwrap();
    let mut got = vec![("count".to_owned(), scored.len().to_string())];
    for (i, (s, f)) in scored.iter().zip(&features.features).enumerate() {
        got.push((format!("{i}.name"), s.name.clone()));
        got.push((format!("{i}.arrays"), s.float_data_arrays.len().to_string()));
        for a in &s.float_data_arrays[1..5] {
            got.push((format!("{i}.{}", a.name), bits(f64::from(a.data[0]))));
        }
        let p = format!("feature.{i}");
        got.push((
            format!("{p}.transition_name"),
            f.metadata["transition_name"].as_str().unwrap().to_owned(),
        ));
        got.push((format!("{p}.intensity"), bits(f64::from(f.intensity))));
        for key in ["log10_total_tic", "inverse_avgFWHM", "avgFWHM", "avgSNR"] {
            got.push((
                format!("{p}.{key}"),
                bits(f.metadata[key].as_f64().unwrap()),
            ));
        }
        got.push((
            format!("{p}.subordinates"),
            f.subordinates.len().to_string(),
        ));
        for (j, sub) in f.subordinates.iter().enumerate() {
            got.push((
                format!("{p}.sub.{j}"),
                format!(
                    "{},{},{}",
                    bits(sub.mz),
                    bits(f64::from(sub.intensity)),
                    bits(sub.metadata["FWHM"].as_f64().unwrap())
                ),
            ));
        }
    }
    for threshold in [0.7, 12.0, 15.0] {
        let tse = set(&[
            &class_params()[..],
            &[("min_select_score", ParamValue::Float(threshold))],
        ]
        .concat());
        let (selected, _) = tse.select_spectra(&scored, &features).unwrap();
        // std::to_string(double) spells the threshold with six decimals.
        let key = format!("select.{threshold:.6}");
        got.push((format!("{key}.count"), selected.len().to_string()));
        for (i, s) in selected.iter().enumerate() {
            got.push((
                format!("{key}.{i}"),
                format!(
                    "{}|{}",
                    s.name,
                    bits(f64::from(s.float_data_arrays[1].data[0]))
                ),
            ));
        }
    }
    let score_float = |key: &str| {
        if key.ends_with(".name")
            || key.ends_with(".arrays")
            || key.ends_with(".count")
            || key.ends_with(".transition_name")
            || key.ends_with(".subordinates")
        {
            vec![]
        } else if key.starts_with("select.") {
            vec![1]
        } else {
            vec![0, 1, 2]
        }
    };
    check(&mut cmp, "score", got, &score_float);

    let tse = set(&[
        &class_params()[..],
        &[("min_select_score", ParamValue::Float(12.0))],
    ]
    .concat());
    let (extracted, extracted_features) = tse.extract_spectra(&experiment, &targets).unwrap();
    let mut got = vec![("count".to_owned(), extracted.len().to_string())];
    for (i, (s, f)) in extracted
        .iter()
        .zip(&extracted_features.features)
        .enumerate()
    {
        got.push((
            i.to_string(),
            format!(
                "{}|{}|{}",
                s.name,
                bits(f64::from(s.float_data_arrays[1].data[0])),
                bits(f64::from(f.intensity))
            ),
        ));
    }
    check(&mut cmp, "extract", got, &|key| {
        if key == "count" { vec![] } else { vec![1, 2] }
    });

    // Library matching with the experiment's own spectra as the library.
    let tse = set(&[
        ("top_matches_to_report", ParamValue::Integer(3)),
        ("min_match_score", ParamValue::Float(0.5)),
    ]);
    let mut comparator = BinnedSpectrumComparator::new().unwrap();
    let options = BTreeMap::from([
        ("bin_size".to_owned(), MetaValue::try_from(1.0).unwrap()),
        ("peak_spread".to_owned(), MetaValue::from(0_i64)),
        ("bin_offset".to_owned(), MetaValue::try_from(0.4).unwrap()),
    ]);
    comparator.init(&experiment.spectra, &options).unwrap();
    let mut got = Vec::new();
    for i in (0..experiment.spectra.len()).step_by(7) {
        let matches = tse
            .match_spectrum(&experiment.spectra[i], &comparator)
            .unwrap();
        got.push((format!("{i}.count"), matches.len().to_string()));
        for (j, m) in matches.iter().enumerate() {
            got.push((
                format!("{i}.{j}"),
                format!("{}|{}", m.spectrum.native_id, bits(m.score)),
            ));
        }
    }
    check(&mut cmp, "match", got, &|key| {
        if key.ends_with(".count") {
            vec![]
        } else {
            vec![1]
        }
    });

    assert_eq!(expected.len(), 8, "every oracle case is replayed");
    assert!(
        cmp.compared > 1000,
        "only {} records compared",
        cmp.compared
    );
    eprintln!(
        "{} records; {} of {} floating values bit-identical; worst relative difference {:e}",
        cmp.compared, cmp.identical_floats, cmp.floats, cmp.worst
    );
}

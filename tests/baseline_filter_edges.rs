// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Executed C++ evidence for `MorphologicalFilter` at the ends of a spectrum.
//!
//! Every expected number in this file was produced by the OpenMS4 Release build
//! (core `bc9cc12`, `x86_64` Linux) running the pinned, header-only
//! `MorphologicalFilter` through the oracle driver
//! `../oracle/baseline-filter-edges/drivers/morph_edges.cpp`, one forked child
//! per case so the source's process-wide `static` buffers start empty. The
//! fixtures, their derivation and the run's hashes are recorded in
//! `tests/data/baseline_filter_edges_provenance.json`; the findings are
//! explained in `docs/MORPHOLOGICAL_FILTER_SUPPORT.md`.

use openms::kernel::SpectrumType;
use openms::processing::SpectrumFilter;
use openms::processing::baseline::{
    MorphologicalFilter, MorphologicalMethod as Method, StructuringElement as Element,
};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};
use std::collections::BTreeMap;

const SHAPES: &str = include_str!("data/baseline_filter_edges_shapes.tsv");
const EXPECTED: &str = include_str!("data/baseline_filter_edges_expected.tsv");
const SWEEP: &str = include_str!("data/baseline_filter_edges_sweep.tsv");

/// One executed case: the run identity and the C++ result.
struct Case<'a> {
    api: &'a str,
    group: &'a str,
    shape: &'a str,
    method: Method,
    unit: &'a str,
    length: f64,
    output: Vec<f32>,
    spectrum_type: &'a str,
}

fn source_name(method: Method) -> &'static str {
    match method {
        Method::Identity => "identity",
        Method::Erosion => "erosion",
        Method::Dilation => "dilation",
        Method::Opening => "opening",
        Method::Closing => "closing",
        Method::Gradient => "gradient",
        Method::TopHat => "tophat",
        Method::BottomHat => "bothat",
        Method::ErosionSimple => "erosion_simple",
        Method::DilationSimple => "dilation_simple",
    }
}

fn method(name: &str) -> Method {
    match name {
        "identity" => Method::Identity,
        "erosion" => Method::Erosion,
        "dilation" => Method::Dilation,
        "opening" => Method::Opening,
        "closing" => Method::Closing,
        "gradient" => Method::Gradient,
        "tophat" => Method::TopHat,
        "bothat" => Method::BottomHat,
        "erosion_simple" => Method::ErosionSimple,
        "dilation_simple" => Method::DilationSimple,
        other => panic!("unknown method {other}"),
    }
}

fn floats<T: std::str::FromStr>(text: &str) -> Vec<T>
where
    T::Err: std::fmt::Debug,
{
    if text == "-" {
        return Vec::new();
    }
    text.split(',').map(|v| v.parse().unwrap()).collect()
}

/// `shape -> (type, m/z, intensities)`, as the driver read them.
fn shapes() -> BTreeMap<&'static str, (SpectrumType, Vec<f64>, Vec<f32>)> {
    SHAPES
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let kind = match f[1] {
                "centroid" => SpectrumType::Centroid,
                "profile" => SpectrumType::Profile,
                other => panic!("unknown type {other}"),
            };
            (f[0], (kind, floats(f[2]), floats(f[3])))
        })
        .collect()
}

fn cases() -> Vec<Case<'static>> {
    EXPECTED
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let id: Vec<_> = f[0].split('|').collect();
            Case {
                api: id[0],
                group: id[1],
                shape: id[2],
                method: method(id[3]),
                unit: id[4],
                length: id[5].parse().unwrap(),
                output: floats(f[1]),
                spectrum_type: f[2],
            }
        })
        .collect()
}

impl Case<'_> {
    /// The element the source built from `struc_elem_unit` and
    /// `struc_elem_length`: `Thomson` keeps the width, `DataPoints` truncates
    /// the double to a count, as the source's `(UInt)(double)` cast does.
    fn element(&self) -> Element {
        match self.unit {
            "Thomson" => Element::Thomson(self.length),
            "DataPoints" => Element::DataPoints(self.length as usize),
            other => panic!("unknown unit {other}"),
        }
    }

    /// This port refuses two inputs the source filters; see the support
    /// document. Both are recorded here so the C++ result stays visible.
    fn refusal(&self) -> Option<Refusal> {
        if self.shape.starts_with("unsorted") {
            Some(Refusal::Unsorted)
        } else if matches!(self.element(), Element::DataPoints(0))
            || matches!(self.element(), Element::Thomson(w) if w <= 0.0)
        {
            Some(Refusal::ZeroElement)
        } else {
            None
        }
    }
}

enum Refusal {
    /// The source documents ascending m/z as a precondition and does not check
    /// it; this port returns [`Error::UnsortedData`].
    Unsorted,
    /// The source turns an element of zero samples into one sample; this port
    /// refuses the parameter.
    ZeroElement,
}

fn spectrum(kind: SpectrumType, mz: &[f64], intensity: &[f32]) -> MSSpectrum {
    let mut s = MSSpectrum::from_peaks(
        mz.iter()
            .zip(intensity)
            .map(|(&m, &i)| Peak1D::new(m, i))
            .collect(),
    );
    s.spectrum_type = kind;
    s
}

fn values(s: &MSSpectrum) -> Vec<f32> {
    s.peaks.iter().map(|p| p.intensity).collect()
}

/// Bitwise equality, so a sign or a last-sample difference cannot pass.
fn assert_same(got: &[f32], want: &[f32], at: &str) {
    assert_eq!(got.len(), want.len(), "{at}: length");
    for (index, (g, w)) in got.iter().zip(want).enumerate() {
        assert_eq!(g.to_bits(), w.to_bits(), "{at}: sample {index}: {g} != {w}");
    }
}

/// The element length source `filter(MSSpectrum&)` computes, replicated for the
/// unsorted shapes, whose positions this port refuses to convert.
fn source_window(unit: &str, length: f64, mz: &[f64]) -> usize {
    let points = match unit {
        "DataPoints" => length as usize,
        _ => (length * (mz.len() - 1) as f64 / (mz[mz.len() - 1] - mz[0])).ceil() as usize,
    };
    if points % 2 == 0 { points + 1 } else { points }
}

#[test]
fn executed_spectrum_cases_match_the_release_build() {
    let shapes = shapes();
    let mut compared = 0;
    for case in cases().iter().filter(|c| c.api == "spectrum") {
        let (kind, mz, intensity) = &shapes[case.shape];
        let at = format!(
            "{} {:?} {} {}",
            case.shape, case.method, case.unit, case.length
        );
        let mut produced = spectrum(*kind, mz, intensity);
        let before = produced.clone();
        let filter = MorphologicalFilter {
            method: case.method,
            structuring_element: case.element(),
        };
        match case.refusal() {
            None => {
                filter.filter_spectrum(&mut produced).expect(&at);
                assert_same(&values(&produced), &case.output, &at);
                assert_eq!(case.spectrum_type, "profile");
                assert_eq!(produced.spectrum_type, SpectrumType::Profile, "{at}");
                assert_eq!(
                    produced.peaks.iter().map(|p| p.mz).collect::<Vec<_>>(),
                    *mz,
                    "{at}: m/z"
                );
                compared += 1;
            }
            Some(Refusal::Unsorted) => {
                assert!(
                    matches!(
                        filter.filter_spectrum(&mut produced),
                        Err(Error::UnsortedData)
                    ),
                    "{at}: expected the unsorted refusal"
                );
                assert_eq!(produced, before, "{at}: refused call changed the spectrum");
                // The same operation over the same sample order, without the
                // precondition check, reproduces the C++ result exactly.
                if intensity.len() > 1 {
                    let points = source_window(case.unit, case.length, mz);
                    let ranged = MorphologicalFilter {
                        method: case.method,
                        structuring_element: Element::DataPoints(points),
                    };
                    assert_same(&ranged.filter_range(intensity).unwrap(), &case.output, &at);
                    compared += 1;
                }
            }
            Some(Refusal::ZeroElement) => {
                assert!(filter.filter_spectrum(&mut produced).is_err(), "{at}");
                assert_eq!(produced, before, "{at}: refused call changed the spectrum");
                // The source's own substitute for the refused parameter.
                let single = MorphologicalFilter {
                    method: case.method,
                    structuring_element: Element::DataPoints(1),
                };
                let mut substitute = spectrum(*kind, mz, intensity);
                single.filter_spectrum(&mut substitute).expect(&at);
                assert_same(&values(&substitute), &case.output, &at);
                compared += 1;
            }
        }
    }
    assert_eq!(compared, 1240, "every recorded spectrum case");
}

#[test]
fn executed_range_cases_match_the_release_build() {
    let shapes = shapes();
    let all = cases();
    let odd: BTreeMap<_, _> = all
        .iter()
        .filter(|c| c.api == "range")
        .map(|c| {
            (
                (c.shape, source_name(c.method), c.length as usize),
                &c.output,
            )
        })
        .collect();
    let mut compared = 0;
    for case in all.iter().filter(|c| c.api == "range") {
        let (_, _, intensity) = &shapes[case.shape];
        let points = case.length as usize;
        let at = format!("{} {:?} DataPoints {points}", case.shape, case.method);
        let filter = MorphologicalFilter {
            method: case.method,
            structuring_element: Element::DataPoints(points),
        };
        // Source `filterRange` uses an even count as given, producing the
        // asymmetric windows of van Herk's blocks; this port rounds up to odd,
        // so an even count must equal the executed result for count + 1.
        let expected = if points % 2 == 0 {
            odd[&(case.shape, source_name(case.method), points + 1)]
        } else {
            &case.output
        };
        assert_same(&filter.filter_range(intensity).unwrap(), expected, &at);
        compared += 1;
    }
    assert_eq!(compared, 350, "every recorded range case");
}

#[test]
fn executed_experiment_groups_reproduce_the_source_buffer_history() {
    let shapes = shapes();
    let all = cases();
    let mut groups: Vec<Vec<&Case>> = Vec::new();
    for case in all.iter().filter(|c| c.api == "experiment") {
        match groups.last() {
            Some(last) if last[0].group == case.group => groups.last_mut().unwrap().push(case),
            _ => groups.push(vec![case]),
        }
    }
    assert_eq!(groups.len(), 12, "one group per method plus two orderings");
    for group in groups {
        let first = group[0];
        let filter = MorphologicalFilter {
            method: first.method,
            structuring_element: first.element(),
        };
        let mut experiment = MSExperiment {
            spectra: group
                .iter()
                .map(|c| {
                    let (kind, mz, intensity) = &shapes[c.shape];
                    spectrum(*kind, mz, intensity)
                })
                .collect(),
            ..Default::default()
        };
        filter.filter_experiment(&mut experiment).unwrap();
        for (case, produced) in group.iter().zip(&experiment.spectra) {
            assert_same(
                &values(produced),
                &case.output,
                &format!("{} {}", case.group, case.shape),
            );
        }
    }
}

/// The benchmark's blocker: spectrum index 3 of
/// `inputs/derived/sub_profile_uk222_first600.mzML`, a centroided MS2 spectrum
/// whose 11 peaks span 596 Th, with the benchmark's 1 Th element. The element
/// is one sample wide, so the source top-hat zeroes every peak but the last.
#[test]
fn the_uk222_reproducer_keeps_its_last_peak() {
    let shapes = shapes();
    let (kind, mz, intensity) = &shapes["reproducer_uk222_s3"];
    let mut produced = spectrum(*kind, mz, intensity);
    MorphologicalFilter::new(Method::TopHat, Element::Thomson(1.0))
        .unwrap()
        .filter_spectrum(&mut produced)
        .unwrap();
    let mut expected = vec![0.0; 10];
    expected.push(intensity[10]);
    // The C++ tool wrote 882.464233 there, which is this binary32 value.
    assert_eq!(expected[10], 882.464_23_f32, "the last input peak");
    assert_same(&values(&produced), &expected, "reproducer");
    // Erosion and dilation zero that sample instead, and the direct-window
    // methods keep every sample.
    for (method, expected) in [
        (Method::Erosion, {
            let mut v = intensity.clone();
            v[10] = 0.0;
            v
        }),
        (Method::ErosionSimple, intensity.clone()),
    ] {
        let mut produced = spectrum(*kind, mz, intensity);
        MorphologicalFilter::new(method, Element::Thomson(1.0))
            .unwrap()
            .filter_spectrum(&mut produced)
            .unwrap();
        assert_same(&values(&produced), &expected, "reproducer");
    }
}

// ---- the executed sweep ----------------------------------------------------

/// The driver's 31-bit LCG: small integers, exact in binary32.
fn sweep_input(n: usize) -> Vec<f32> {
    let mut x = n as u64;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(1103515245).wrapping_add(12345) % 2147483648;
            ((x >> 16) % 11) as i64 as f32 - 5.0
        })
        .collect()
}

fn clipped(values: &[f32], length: usize, maximum: bool) -> Vec<f32> {
    let radius = length / 2;
    (0..values.len())
        .map(|center| {
            let start = center.saturating_sub(radius);
            let end = (center + radius + 1).min(values.len());
            values[start..end]
                .iter()
                .copied()
                .reduce(|a, b| if maximum { a.max(b) } else { a.min(b) })
                .unwrap()
        })
        .collect()
}

/// The clipped-window result the driver compared against, in `f32` arithmetic.
fn reference(method: Method, values: &[f32], length: usize) -> Vec<f32> {
    let erosion = |v: &[f32]| clipped(v, length, false);
    let dilation = |v: &[f32]| clipped(v, length, true);
    let subtract = |a: Vec<f32>, b: Vec<f32>| a.iter().zip(b).map(|(x, y)| x - y).collect();
    match method {
        Method::Identity => values.to_vec(),
        Method::Erosion | Method::ErosionSimple => erosion(values),
        Method::Dilation | Method::DilationSimple => dilation(values),
        Method::Opening => dilation(&erosion(values)),
        Method::Closing => erosion(&dilation(values)),
        Method::Gradient => subtract(dilation(values), erosion(values)),
        Method::TopHat => subtract(values.to_vec(), dilation(&erosion(values))),
        Method::BottomHat => subtract(values.to_vec(), erosion(&dilation(values))),
    }
}

const METHODS: [Method; 10] = [
    Method::Identity,
    Method::Erosion,
    Method::Dilation,
    Method::Opening,
    Method::Closing,
    Method::Gradient,
    Method::TopHat,
    Method::BottomHat,
    Method::ErosionSimple,
    Method::DilationSimple,
];

/// Every place where the Release build departs from a clipped-window filter,
/// and nowhere else: the executed sweep of 51,780 cases over 3.3 million
/// samples recorded exactly these samples, and this port reproduces them.
#[test]
fn executed_sweep_mismatches_are_reproduced_exactly() {
    let mut expected: BTreeMap<(String, String, usize, usize, usize), f32> = BTreeMap::new();
    for line in SWEEP.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        expected.insert(
            (
                f[1].to_string(),
                f[2].to_string(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
            ),
            f[6].parse().unwrap(),
        );
    }
    assert_eq!(expected.len(), 528, "the recorded mismatches");

    let mut produced: BTreeMap<(&str, &str, usize, usize, usize), f32> = BTreeMap::new();
    fn record(
        into: &mut BTreeMap<(&'static str, &'static str, usize, usize, usize), f32>,
        api: &'static str,
        method: Method,
        n: usize,
        length: usize,
        out: &[f32],
    ) {
        let effective = if length % 2 == 0 { length + 1 } else { length };
        let input = sweep_input(n);
        // Source `filter(MSSpectrum&)` returns before touching a spectrum with
        // fewer than two peaks, where `filterRange` filters it.
        let reference = if api == "spectrum" && n <= 1 {
            input
        } else {
            reference(method, &input, effective)
        };
        for (index, (got, want)) in out.iter().zip(&reference).enumerate() {
            if got.to_bits() != want.to_bits() {
                into.insert((api, source_name(method), n, length, index), *got);
            }
        }
    }
    let sizes = (0..=48).chain([97, 128, 255, 256, 1000, 4097]);
    for n in sizes {
        let input = sweep_input(n);
        let lengths: Vec<usize> = if n <= 48 {
            (1..=2 * n + 3).collect()
        } else {
            vec![1, 2, 3, 5, 7, 9, 15, 31, 63, 99, 255, 257, 1001, 4097, 4099]
        };
        for method in METHODS {
            for &length in &lengths {
                let filter = MorphologicalFilter {
                    method,
                    structuring_element: Element::DataPoints(length),
                };
                let mut s = spectrum(
                    SpectrumType::Profile,
                    &(0..n).map(|i| i as f64).collect::<Vec<_>>(),
                    &input,
                );
                filter.filter_spectrum(&mut s).unwrap();
                record(&mut produced, "spectrum", method, n, length, &values(&s));
                if length % 2 == 1 {
                    record(
                        &mut produced,
                        "range",
                        method,
                        n,
                        length,
                        &filter.filter_range(&input).unwrap(),
                    );
                }
            }
        }
    }
    for (key, value) in &expected {
        let borrowed = (key.0.as_str(), key.1.as_str(), key.2, key.3, key.4);
        let got = produced
            .get(&borrowed)
            .unwrap_or_else(|| panic!("missing mismatch {borrowed:?}, C++ gave {value}"));
        assert_eq!(got.to_bits(), value.to_bits(), "{borrowed:?}");
    }
    for key in produced.keys() {
        let owned = (key.0.to_string(), key.1.to_string(), key.2, key.3, key.4);
        assert!(
            expected.contains_key(&owned),
            "{key:?} differs from the clipped reference, but C++ did not"
        );
    }
}

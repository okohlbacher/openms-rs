// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! Tier-1 differential for `ANALYSIS/ID/PScore.h` and `ANALYSIS/ID/AScore.h`.
//!
//! `tests/data/pscore_ascore_oracle_cases.tsv` lists windowed ranks, peak-level
//! spectra, both `computePScore` overloads and `AScore::compute` over the five
//! upstream DTA inputs; `tests/data/pscore_ascore_oracle.tsv` is what the
//! Release libOpenMS (`bc9cc12-c19e494-174b576`, Linux x86_64, run on
//! ibminode06) printed for them (driver and generator in `../oracle/I3b-scores`,
//! hashes in `tests/data/pscore_ascore_provenance.json`).
//!
//! Ranks, level spectra, theoretical sizes, sequences and metadata keys must
//! match exactly. Scores are compared to 1e-9 relative: both builds evaluate
//! `-10 log10` of the same binomial upper tail, the source through Boost's
//! incomplete beta and this port by an explicit log-domain sum
//! (`concept::math_functions::binomial_cdf_complement` documents about 2.5e-14
//! relative agreement), so bit equality is not the contract.

use openms::MSSpectrum;
use openms::analysis::id::ascore::AScore;
use openms::analysis::id::pscore::PScore;
use openms::chemistry::AASequence;
use openms::chemistry::theoretical::TheoreticalSpectrumGenerator;
use openms::format::dta;
use openms::identification::PeptideHit;
use openms::metadata::MetaValueData;
use openms::param::{Param, ParamValue};
use std::collections::BTreeMap;

fn spectrum(name: &str) -> MSSpectrum {
    let bytes: &[u8] = match name {
        "PScore_test.dta" => include_bytes!("data/PScore_test.dta"),
        "Ascore_test_input1.dta" => include_bytes!("data/Ascore_test_input1.dta"),
        "Ascore_test_input2.dta" => include_bytes!("data/Ascore_test_input2.dta"),
        "Ascore_test_input3.dta" => include_bytes!("data/Ascore_test_input3.dta"),
        "Ascore_test_input4.dta" => include_bytes!("data/Ascore_test_input4.dta"),
        other => panic!("unknown DTA {other}"),
    };
    dta::read(bytes).unwrap()
}

fn ranks(spec: &MSSpectrum, window: f64) -> Vec<usize> {
    let mz: Vec<f64> = spec.peaks.iter().map(|p| p.mz).collect();
    let intensity: Vec<f64> = spec.peaks.iter().map(|p| f64::from(p.intensity)).collect();
    PScore::calculate_intensity_rank_in_mz_window(&mz, &intensity, window).unwrap()
}

fn close(actual: f64, expected: f64, what: &str) {
    if expected.is_infinite() || actual.is_infinite() {
        assert_eq!(actual, expected, "{what}");
        return;
    }
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= 1e-9 * scale,
        "{what}: {actual} vs {expected}"
    );
}

fn number(text: &str) -> f64 {
    match text {
        "inf" => f64::INFINITY,
        "-inf" => f64::NEG_INFINITY,
        _ => text.parse().unwrap(),
    }
}

/// The oracle's records for one case id, in output order.
fn oracle() -> BTreeMap<String, Vec<(String, String)>> {
    let mut map: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    for line in include_str!("data/pscore_ascore_oracle.tsv").lines() {
        let mut fields = line.splitn(3, '\t');
        let id = fields.next().unwrap().to_owned();
        let key = fields.next().unwrap().to_owned();
        let value = fields.next().unwrap_or("").to_owned();
        map.entry(id).or_default().push((key, value));
    }
    map
}

#[test]
fn release_oracle_scores() {
    let expected = oracle();
    let mut cases = 0;
    let mut values = 0;
    for line in include_str!("data/pscore_ascore_oracle_cases.tsv").lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (kind, id) = (f[0], f[1]);
        let records = &expected[id];
        cases += 1;
        match kind {
            "ranks" => {
                let got = ranks(&spectrum(f[2]), number(f[3]));
                let text: Vec<String> = got.iter().map(usize::to_string).collect();
                assert_eq!(records, &vec![("ranks".to_owned(), text.join(","))], "{id}");
                values += got.len();
            }
            "levels" => {
                let spec = spectrum(f[2]);
                let r = ranks(&spec, number(f[3]));
                let levels = PScore::calculate_peak_level_spectra(
                    &spec,
                    &r,
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                )
                .unwrap();
                assert_eq!(levels.len(), records.len(), "{id}");
                for ((level, s), (key, value)) in levels.iter().zip(records) {
                    assert_eq!(key, &format!("level{level}"), "{id}");
                    let mz: Vec<f64> = value.split(',').map(number).collect();
                    let got: Vec<f64> = s.peaks.iter().map(|p| p.mz).collect();
                    assert_eq!(got, mz, "{id} level {level}");
                    values += mz.len();
                }
            }
            "pscore" => {
                let spec = spectrum(f[2]);
                let generator = TheoreticalSpectrumGenerator {
                    add_first_prefix_ion: f[4] == "1",
                    ..TheoreticalSpectrumGenerator::default()
                };
                let peptide = AASequence::parse(f[3]).unwrap();
                let theo = generator.generate(&peptide, 1, 1, None).unwrap();
                let theo2 = generator.generate(&peptide, 1, 2, None).unwrap();
                let r = ranks(&spec, number(f[5]));
                let pls = PScore::calculate_peak_level_spectra(
                    &spec,
                    &r,
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                )
                .unwrap();
                let (tol, ppm, window) = (number(f[8]), f[9] == "1", number(f[10]));
                let single = PScore::compute_pscore(tol, ppm, &pls, &theo, window).unwrap();
                let best =
                    PScore::compute_pscore_best_of(tol, ppm, &pls, &[theo.clone(), theo2], window)
                        .unwrap();
                let get = |k: &str| &records.iter().find(|(key, _)| key == k).unwrap().1;
                assert_eq!(get("theo_size"), &theo.len().to_string(), "{id}");
                close(single, number(get("single")), id);
                close(best, number(get("best_of")), id);
                values += 3;
            }
            "ascore" => {
                let mut param = Param::new();
                param
                    .set_value(
                        "fragment_mass_tolerance",
                        ParamValue::Float(number(f[4])),
                        "",
                        &[],
                    )
                    .unwrap();
                param
                    .set_value(
                        "fragment_mass_unit",
                        ParamValue::String(f[5].into()),
                        "",
                        &[],
                    )
                    .unwrap();
                let ascore = AScore::from_param(&param).unwrap();
                let hit = PeptideHit::new(1.0, 1, 1, AASequence::parse(f[3]).unwrap()).unwrap();
                let mut spec = spectrum(f[2]);
                let result = ascore.compute(&hit, &mut spec);
                if records[0].0 == "error" {
                    assert!(result.is_err(), "{id}: the source threw {}", records[0].1);
                    values += 1;
                    continue;
                }
                let out = result.unwrap();
                let mut expected_meta = BTreeMap::new();
                for (key, value) in records {
                    match key.as_str() {
                        "score" => close(out.score, number(value), id),
                        "sequence" => assert_eq!(&out.sequence.to_string(), value, "{id}"),
                        other => {
                            let name = other.strip_prefix("meta:").unwrap();
                            expected_meta.insert(name.to_owned(), value.clone());
                        }
                    }
                    values += 1;
                }
                let got_keys: Vec<&String> = out.metadata.keys().collect();
                let want_keys: Vec<&String> = expected_meta.keys().collect();
                assert_eq!(got_keys, want_keys, "{id}");
                for (key, value) in &expected_meta {
                    match out.metadata[key].data() {
                        MetaValueData::Float(v) => close(*v, number(value), id),
                        MetaValueData::String(s) => assert_eq!(s, value, "{id} {key}"),
                        other => panic!("{id} {key}: unexpected {other:?}"),
                    }
                }
            }
            other => panic!("unknown kind {other}"),
        }
    }
    assert_eq!(cases, expected.len());
    assert!(values > 1000, "compared {values} values");
}

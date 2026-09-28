// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `FEATUREFINDER/Biosaur2Algorithm.h` (`src/analysis/featurefinder/biosaur2.rs`).
//!
//! Evidence (`tests/data/biosaur2_provenance.json`):
//!
//! - **tier 1, executed C++**: the Linux x86_64 **Release** build
//!   `openms4-release-bc9cc12-c19e494-174b576` on ibminode06, through the oracle
//!   driver `../oracle/FF4-biosaur2/driver.cpp`, which reaches the private
//!   helpers with `#define private public` and prints every double as its
//!   IEEE-754 bit pattern. Every case ran twice in fresh processes and
//!   reproduced byte for byte. `tests/data/biosaur2/oracle.tsv.gz` holds, per
//!   case, the input digest, the status, every hill (summary bits and a digest
//!   of its per-point vectors), every peptide feature, every `FeatureMap` entry
//!   (position, intensity, charge, qualities, meta values, hull points), the
//!   processed experiment (digest and ranges) and both TSV exports;
//!   `helpers_oracle.tsv` holds the helper probes, the Boost averagine table
//!   and the parameter defaults. The 47 cases cover the upstream class test's
//!   inputs, randomised LC-MS runs (ties, zero intensities, doublet peaks),
//!   every parameter branch (`iuse`, `nm`, `hrttol` 0, both calibrations,
//!   `convex_hulls`, `hvf`/`ivf`/`minlh`, charge ranges, windows), profile, TOF,
//!   FAIMS (with and without per-peak ion mobility), PASEF/TIMS ion mobility
//!   (full, partial, per-spectrum drift time) and five upstream mzML files;
//! - **tier 3, source review**: the literals of `Biosaur2Algorithm_test.cpp`
//!   (core `bc9cc12`), each cited;
//! - **tier 4, native**: the refusals, the atomicity of a failed run, the
//!   thread-count determinism and the independent check of the averagine table.

#![cfg(feature = "mzml")]

use openms::analysis::featurefinder::biosaur2::{
    Biosaur2Algorithm, Biosaur2LogLevel, ConvexHullMode, Hill, IsotopeCandidate, hills, isotopes,
    output,
};
use openms::concept::parallel::Threads;
use openms::concept::unique_id::UniqueIdGenerator;
use openms::kernel::faims_helper::FaimsHelper;
use openms::kernel::ranges::MSDim;
use openms::kernel::{DataArray, MSChromatogram, SpectrumType};
use openms::metadata::{DriftTimeUnit, MetaValueData, ProcessingAction};
use openms::param::ParamValue;
use openms::{Error, FeatureMap, MSExperiment, MSSpectrum, Peak1D};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn fixture_bytes(name: &str) -> Vec<u8> {
    let bytes = std::fs::read(data_dir().join("biosaur2").join(name)).unwrap();
    if name.ends_with(".gz") {
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(bytes.as_slice())
            .read_to_end(&mut out)
            .unwrap();
        out
    } else {
        bytes
    }
}

fn fixture(name: &str) -> String {
    String::from_utf8(fixture_bytes(name)).unwrap()
}

fn bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

/// The FNV-1a digest the oracle driver computes.
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Self(1_469_598_103_934_665_603)
    }
    fn bytes(&mut self, data: &[u8]) {
        for &byte in data {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(1_099_511_628_211);
        }
    }
    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }
    fn f64(&mut self, value: f64) {
        self.bytes(&value.to_le_bytes());
    }
    fn f32(&mut self, value: f32) {
        self.bytes(&value.to_le_bytes());
    }
    fn str(&mut self, value: &str) {
        self.u64(value.len() as u64);
        self.bytes(value.as_bytes());
    }
    fn hex(&self) -> String {
        format!("{:016x}", self.0)
    }
}

fn unit_code(unit: DriftTimeUnit) -> u64 {
    match unit {
        DriftTimeUnit::None => 0,
        DriftTimeUnit::Millisecond => 1,
        DriftTimeUnit::InverseReducedMobility => 2,
        DriftTimeUnit::FaimsCompensationVoltage => 3,
        DriftTimeUnit::CollisionCrossSection => 4,
    }
}

/// `MSSpectrum::getType(false)`: the stored type, or centroid when a
/// processing record says the data was picked.
fn type_code(spectrum: &MSSpectrum) -> u64 {
    match spectrum.spectrum_type {
        SpectrumType::Centroid => 1,
        SpectrumType::Profile => 2,
        SpectrumType::Unknown => {
            if spectrum
                .data_processing
                .iter()
                .any(|dp| dp.actions.contains(&ProcessingAction::PeakPicking))
            {
                1
            } else {
                0
            }
        }
    }
}

fn spectra_digest(experiment: &MSExperiment) -> String {
    let mut f = Fnv::new();
    f.u64(experiment.spectra.len() as u64);
    for s in &experiment.spectra {
        f.f64(s.rt);
        f.u64(u64::from(s.ms_level));
        f.f64(s.drift_time);
        f.u64(unit_code(s.drift_time_unit));
        f.u64(type_code(s));
        f.u64(s.peaks.len() as u64);
        for p in &s.peaks {
            f.f64(p.mz);
            f.f32(p.intensity);
        }
        f.u64(s.float_data_arrays.len() as u64);
        for a in &s.float_data_arrays {
            f.str(&a.name);
            f.u64(a.data.len() as u64);
            for &v in &a.data {
                f.f32(v);
            }
        }
    }
    f.hex()
}

fn hill_digest(h: &Hill) -> String {
    let mut f = Fnv::new();
    f.u64(h.scan_indices.len() as u64);
    for &v in &h.scan_indices {
        f.u64(v as u64);
    }
    f.u64(h.peak_indices.len() as u64);
    for &v in &h.peak_indices {
        f.u64(v as u64);
    }
    for values in [
        &h.mz_values,
        &h.intensities,
        &h.rt_values,
        &h.drift_times,
        &h.ion_mobilities,
    ] {
        f.u64(values.len() as u64);
        for &v in values {
            f.f64(v);
        }
    }
    f.hex()
}

/// The oracle's plain-text experiment format (`../oracle/FF4-biosaur2/gen_cases.py`).
fn read_exp(text: &str) -> MSExperiment {
    let mut experiment = MSExperiment::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.first() {
            Some(&"S") => {
                let unit = match fields[4] {
                    "none" => DriftTimeUnit::None,
                    "ms" => DriftTimeUnit::Millisecond,
                    "1/K0" => DriftTimeUnit::InverseReducedMobility,
                    "FAIMS_CV" => DriftTimeUnit::FaimsCompensationVoltage,
                    "CCS" => DriftTimeUnit::CollisionCrossSection,
                    other => panic!("unit {other}"),
                };
                let spectrum_type = match fields[5] {
                    "centroid" => SpectrumType::Centroid,
                    "profile" => SpectrumType::Profile,
                    _ => SpectrumType::Unknown,
                };
                let n: usize = fields[6].parse().unwrap();
                let im_name = fields[7];
                let mut spectrum = MSSpectrum {
                    rt: fields[1].parse().unwrap(),
                    ms_level: fields[2].parse().unwrap(),
                    drift_time: fields[3].parse().unwrap(),
                    drift_time_unit: unit,
                    spectrum_type,
                    ..MSSpectrum::default()
                };
                let mut im = Vec::new();
                for _ in 0..n {
                    let values: Vec<f64> = lines
                        .next()
                        .unwrap()
                        .split_whitespace()
                        .map(|v| v.parse().unwrap())
                        .collect();
                    spectrum
                        .peaks
                        .push(Peak1D::new(values[0], values[1] as f32));
                    if im_name != "-" {
                        im.push(values[2] as f32);
                    }
                }
                if im_name != "-" {
                    spectrum
                        .float_data_arrays
                        .push(DataArray::new(im_name.replace('_', " "), im));
                }
                experiment.spectra.push(spectrum);
            }
            Some(&"C") => experiment.chromatograms.push(MSChromatogram {
                native_id: fields[1].to_owned(),
                ..MSChromatogram::default()
            }),
            Some(&"D") => {
                experiment.settings.date_time =
                    openms::data_structures::DateTime::parse(fields[1]).unwrap();
            }
            _ => {}
        }
    }
    experiment
}

fn mzml_path(name: &str) -> PathBuf {
    let dir = match name {
        "FeatureFinderCentroided_1_input.mzML"
        | "FAIMS_test_data.mzML"
        | "FAIMS_CV-60C_V-45_Interleaved.mzML" => "mzml_mobility",
        "IM_FAIMS_test.mzML" => "faims_helper",
        "FeatureFinderAlgorithmPicked.mzML" => "feature_finder_picked",
        other => panic!("unknown mzML {other}"),
    };
    data_dir().join(dir).join(name)
}

fn load_input(input: &str) -> MSExperiment {
    if let Some(name) = input.strip_prefix("mzml:") {
        let file = std::fs::File::open(mzml_path(name)).unwrap();
        openms::format::mzml::read(std::io::BufReader::new(file)).unwrap()
    } else {
        read_exp(&fixture(&format!("{input}.gz")))
    }
}

/// `(input, parameters)` of every case of `cases.tsv`.
fn cases() -> &'static BTreeMap<String, (String, String)> {
    static CASES: OnceLock<BTreeMap<String, (String, String)>> = OnceLock::new();
    CASES.get_or_init(|| {
        fixture("cases.tsv")
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                (
                    f[0].to_owned(),
                    (f[1].to_owned(), f.get(2).copied().unwrap_or("").to_owned()),
                )
            })
            .collect()
    })
}

/// The oracle records, by case, then by kind.
type Records = BTreeMap<String, BTreeMap<String, Vec<Vec<String>>>>;

fn oracle() -> &'static Records {
    static ORACLE: OnceLock<Records> = OnceLock::new();
    ORACLE.get_or_init(|| {
        let mut records = Records::new();
        for line in fixture("oracle.tsv.gz").lines() {
            let f: Vec<String> = line.split('\t').map(str::to_owned).collect();
            let (kind, case) = (f[0].clone(), f[1].clone());
            records
                .entry(case)
                .or_default()
                .entry(kind)
                .or_default()
                .push(f);
        }
        records
    })
}

/// Apply `k=v;...` with each value typed as its default is (the driver's
/// `setParameters` path).
fn configured(parameters: &str) -> Biosaur2Algorithm {
    let mut param = Biosaur2Algorithm::defaults().unwrap();
    for kv in parameters.split(';').filter(|kv| !kv.is_empty()) {
        let (key, value) = kv.split_once('=').unwrap();
        let typed = match param.value(key).unwrap() {
            ParamValue::Float(_) => ParamValue::Float(value.parse().unwrap()),
            ParamValue::Integer(_) => ParamValue::Integer(value.parse().unwrap()),
            _ => ParamValue::String(value.to_owned()),
        };
        param.set_value(key, typed, "", &[]).unwrap();
    }
    Biosaur2Algorithm::from_param(&param).unwrap()
}

fn meta_text(feature: &openms::Feature) -> String {
    if feature.metadata.is_empty() {
        return "-".into();
    }
    let mut out = String::new();
    for (key, value) in &feature.metadata {
        out.push_str(key);
        out.push('=');
        match &value.data {
            MetaValueData::Float(v) => out.push_str(&format!("d:{:016x}", v.to_bits())),
            MetaValueData::Integer(v) => out.push_str(&format!("i:{v}")),
            other => out.push_str(&format!("s:{other:?}")),
        }
        out.push(';');
    }
    out
}

fn hull_text(feature: &openms::Feature) -> String {
    let mut out = String::new();
    for hull in &feature.convex_hulls {
        for p in hull.hull_points() {
            out.push_str(&format!("{:016x}/{:016x},", p.rt.to_bits(), p.mz.to_bits()));
        }
        out.push('|');
    }
    if out.is_empty() { "-".into() } else { out }
}

fn f64_field(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

/// Run one oracle case and compare everything the oracle recorded.
fn check_case(name: &str) {
    let (input, parameters) = &cases()[name];
    let expected = &oracle()[name];
    let experiment = load_input(input);
    let input_record = &expected["input"][0];
    assert_eq!(
        spectra_digest(&experiment),
        input_record[2],
        "{name}: the Rust reader does not load the input the C++ loaded"
    );
    assert_eq!(experiment.spectra.len().to_string(), input_record[3]);

    let mut algorithm = configured(parameters);
    algorithm.set_ms_data(experiment.clone());
    let mut generator = UniqueIdGenerator::from_seed(7);
    let status = &expected["status"][0];
    let output = match algorithm.run_with_details(&mut generator) {
        Ok(output) => {
            assert_eq!(
                status[2], "ok",
                "{name}: the Release build failed: {}",
                status[3]
            );
            output
        }
        Err(error) => {
            assert_eq!(status[2], "error", "{name}: {error}");
            let Error::InvalidValue(message) = &error else {
                panic!("{name}: unexpected error kind {error:?}");
            };
            assert_eq!(format!("InvalidValue: {message}"), status[3], "{name}");
            // Native: a failed run leaves the stored data unchanged (the
            // source has already removed the non-MS1 spectra).
            assert_eq!(algorithm.ms_data(), &experiment, "{name}");
            return;
        }
    };

    // Hills.
    let hills_expected = expected.get("hill").map_or(&[][..], Vec::as_slice);
    assert_eq!(
        output.hills.len(),
        hills_expected.len(),
        "{name}: hill count"
    );
    for (hill, record) in output.hills.iter().zip(hills_expected) {
        let got = [
            hill.hill_idx.to_string(),
            hill.length.to_string(),
            f64_field(hill.mz_weighted_mean),
            f64_field(hill.rt_start),
            f64_field(hill.rt_end),
            f64_field(hill.rt_apex),
            f64_field(hill.intensity_apex),
            f64_field(hill.intensity_sum),
            f64_field(hill.drift_time_median),
            f64_field(hill.ion_mobility_median),
            hill_digest(hill),
        ];
        assert_eq!(got.as_slice(), &record[3..14], "{name}: hill {}", record[2]);
    }

    // Peptide features.
    let pfs_expected = expected.get("pf").map_or(&[][..], Vec::as_slice);
    assert_eq!(
        output.peptide_features.len(),
        pfs_expected.len(),
        "{name}: peptide feature count"
    );
    for (pf, record) in output.peptide_features.iter().zip(pfs_expected) {
        let mut iso = String::new();
        for c in &pf.isotopes {
            iso.push_str(&format!(
                "{}:{}:{}:{};",
                c.hill_idx,
                c.isotope_number,
                f64_field(c.mass_diff_ppm),
                f64_field(c.cos_corr)
            ));
        }
        if iso.is_empty() {
            iso.push('-');
        }
        let got = [
            f64_field(pf.mz),
            f64_field(pf.rt_start),
            f64_field(pf.rt_end),
            f64_field(pf.rt_apex),
            f64_field(pf.intensity_apex),
            f64_field(pf.intensity_sum),
            pf.charge.to_string(),
            pf.n_isotopes.to_string(),
            pf.n_scans.to_string(),
            f64_field(pf.mass_calib),
            f64_field(pf.drift_time),
            f64_field(pf.ion_mobility),
            pf.mono_hill_idx.to_string(),
            iso,
        ];
        assert_eq!(
            got.as_slice(),
            &record[3..17],
            "{name}: peptide feature {}",
            record[2]
        );
    }

    // FeatureMap.
    let map = &output.feature_map;
    let feats_expected = expected.get("feat").map_or(&[][..], Vec::as_slice);
    assert_eq!(
        map.features.len(),
        feats_expected.len(),
        "{name}: feature count"
    );
    for (feature, record) in map.features.iter().zip(feats_expected) {
        let got = [
            f64_field(feature.rt),
            f64_field(feature.mz),
            format!("{:08x}", feature.intensity.to_bits()),
            feature.charge.to_string(),
            format!("{:08x}", feature.quality.to_bits()),
            format!("{:08x}", feature.quality_rt.to_bits()),
            format!("{:08x}", feature.quality_mz.to_bits()),
            format!("{:08x}", feature.width.to_bits()),
            meta_text(feature),
            hull_text(feature),
            feature.subordinates.len().to_string(),
        ];
        assert_eq!(
            got.as_slice(),
            &record[3..14],
            "{name}: feature {}",
            record[2]
        );
        assert_ne!(feature.unique_id, 0, "{name}: feature without unique id");
    }
    let fmap = &expected["fmap"][0];
    assert_eq!(map.features.len().to_string(), fmap[2]);
    assert_eq!(map.protein_identifications.len().to_string(), fmap[3]);
    assert_eq!(
        map.unassigned_peptide_identifications.len().to_string(),
        fmap[4]
    );
    assert_eq!(map.data_processing.len().to_string(), fmap[5]);
    assert_eq!(map.metadata.is_empty(), fmap[6] == "-");
    assert_ne!(map.unique_id, 0);

    // The processed experiment.
    let exp = &expected["exp"][0];
    let back = algorithm.ms_data();
    assert_eq!(
        back.spectra.len().to_string(),
        exp[2],
        "{name}: spectra kept"
    );
    assert_eq!(
        back.chromatograms.len().to_string(),
        exp[3],
        "{name}: chromatograms kept"
    );
    assert_eq!(spectra_digest(back), exp[4], "{name}: processed spectra");
    assert_eq!(
        back.settings.date_time.iso_string(),
        exp[5],
        "{name}: settings"
    );
    assert_eq!(back.chromatograms, experiment.chromatograms, "{name}");
    if exp[6] != "-" {
        let ranges = back.combined_range_manager().unwrap();
        let got = [
            f64_field(ranges.min(MSDim::Rt).unwrap()),
            f64_field(ranges.max(MSDim::Rt).unwrap()),
            f64_field(ranges.min(MSDim::Mz).unwrap()),
            f64_field(ranges.max(MSDim::Mz).unwrap()),
        ];
        assert_eq!(got.as_slice(), &exp[6..10], "{name}: ranges");
    }

    // The two TSV exports.
    let tsv: Vec<String> = expected.get("tsvf").map_or(Vec::new(), |rows| {
        rows.iter().map(|r| r[2..].join("\t")).collect()
    });
    let text = output::tsv_text(&output.peptide_features);
    assert_eq!(text.lines().collect::<Vec<_>>(), tsv, "{name}: feature TSV");
    let hills_text = output::hills_text(&output.hills);
    let mut f = Fnv::new();
    f.bytes(hills_text.as_bytes());
    let tsvh = &expected["tsvh"][0];
    assert_eq!(
        hills_text.len().to_string(),
        tsvh[2],
        "{name}: hill TSV size"
    );
    assert_eq!(f.hex(), tsvh[3], "{name}: hill TSV");

    // No warning is expected on this data but the RT-gating diagnostics.
    for message in &output.messages {
        if message.level == Biosaur2LogLevel::Warning {
            assert!(
                !message.text.contains("Biosaur2 isotope debug"),
                "{name}: {}",
                message.text
            );
        }
    }
}

macro_rules! oracle_cases {
    ($($test:ident => $case:literal,)*) => {
        $(
            #[test]
            fn $test() {
                check_case($case);
            }
        )*
        #[test]
        fn every_oracle_case_has_a_test() {
            let mut tested = vec![$($case),*];
            tested.sort_unstable();
            let all: Vec<&str> = cases().keys().map(String::as_str).collect();
            assert_eq!(tested, all);
        }
    };
}

oracle_cases! {
    release_ct_run_basic => "ct_run_basic",
    release_ct_run_hills => "ct_run_hills",
    release_ct_faims => "ct_faims",
    release_ct_faims_profile => "ct_faims_profile",
    release_ct_profile => "ct_profile",
    release_ct_nonfaims => "ct_nonfaims",
    release_err_no_ms1 => "err_no_ms1",
    release_err_empty => "err_empty",
    release_no_peaks_in_range => "no_peaks_in_range",
    release_syn11 => "syn11",
    release_syn12 => "syn12",
    release_syn13 => "syn13",
    release_syn_big => "syn_big",
    release_syn_big_iuse_all => "syn_big_iuse_all",
    release_syn_big_iuse2_nm => "syn_big_iuse2_nm",
    release_syn_big_hrttol0 => "syn_big_hrttol0",
    release_syn_big_noisocalib => "syn_big_noisocalib",
    release_syn_big_hillcalib => "syn_big_hillcalib",
    release_syn_big_masstraces => "syn_big_masstraces",
    release_syn_big_hvf_ivf => "syn_big_hvf_ivf",
    release_syn_big_c26 => "syn_big_c26",
    release_syn_big_c18 => "syn_big_c18",
    release_syn_big_mini => "syn_big_mini",
    release_syn_ties => "syn_ties",
    release_syn_zero => "syn_zero",
    release_syn_profile => "syn_profile",
    release_syn_tof => "syn_tof",
    release_syn_tof_off => "syn_tof_off",
    release_syn_faims => "syn_faims",
    release_syn_faims_nomerge => "syn_faims_nomerge",
    release_syn_im => "syn_im",
    release_syn_im_notol => "syn_im_notol",
    release_syn_im_pasef => "syn_im_pasef",
    release_syn_im_masstraces => "syn_im_masstraces",
    release_syn_im1 => "syn_im1",
    release_syn_im_partial => "syn_im_partial",
    release_syn_calib => "syn_calib",
    release_syn_calib_hill => "syn_calib_hill",
    release_syn_faims_im => "syn_faims_im",
    release_syn_drift => "syn_drift",
    release_ffc1_topp => "ffc1_topp",
    release_ffc1_default => "ffc1_default",
    release_ffc1_wide => "ffc1_wide",
    release_ffap => "ffap",
    release_faims_test => "faims_test",
    release_faims_interleaved => "faims_interleaved",
    release_im_faims => "im_faims",
}

/// The helper probes of `helpers.tsv` against `helpers_oracle.tsv`.
fn helper_records() -> &'static BTreeMap<(String, usize), Vec<String>> {
    static RECORDS: OnceLock<BTreeMap<(String, usize), Vec<String>>> = OnceLock::new();
    RECORDS.get_or_init(|| {
        fixture("helpers_oracle.tsv")
            .lines()
            .filter_map(|line| {
                let f: Vec<String> = line.split('\t').map(str::to_owned).collect();
                if f[0] == "default" || f[0] == "table" {
                    return None;
                }
                Some(((f[0].clone(), f[1].parse().unwrap()), f[2..].to_vec()))
            })
            .collect()
    })
}

fn vec_of(text: &str) -> Vec<f64> {
    if text == "-" || text.is_empty() {
        Vec::new()
    } else {
        text.split(',').map(|v| v.parse().unwrap()).collect()
    }
}

fn bits_list(values: &[f64]) -> String {
    values
        .iter()
        .map(|&v| format!("{},", f64_field(v)))
        .collect()
}

#[test]
fn private_helpers_match_release_build() {
    let probes = fixture("helpers.tsv");
    let mut checked = 0;
    for (index, line) in probes.lines().enumerate() {
        let lineno = index + 1;
        let f: Vec<&str> = line.split('\t').collect();
        let expected = &helper_records()[&(f[0].to_owned(), lineno)];
        let got: Vec<String> = match f[0] {
            "averagine" => {
                let (theor, max_pos) =
                    isotopes::compute_averagine(f[1].parse().unwrap(), f[2].parse().unwrap());
                vec![max_pos.to_string(), bits_list(&theor)]
            }
            "meanfilter" => vec![bits_list(&hills::mean_filter(
                &vec_of(f[2]),
                f[1].parse().unwrap(),
            ))],
            "median" => vec![f64_field(hills::calculate_median(&vec_of(f[1])).unwrap())],
            "calib" => {
                let (shift, sigma) = hills::calibrate_mass(&vec_of(f[1]), 0.05).unwrap();
                vec![f64_field(shift), f64_field(sigma)]
            }
            "coscarbon" => {
                let (cor, pos) = isotopes::checking_cos_correlation_for_carbon(
                    &vec_of(f[1]),
                    &vec_of(f[2]),
                    f[3].parse().unwrap(),
                );
                vec![f64_field(cor), pos.to_string()]
            }
            "cosine" => {
                let s1: Vec<usize> = vec_of(f[2]).iter().map(|&v| v as usize).collect();
                let s2: Vec<usize> = vec_of(f[4]).iter().map(|&v| v as usize).collect();
                vec![
                    f64_field(isotopes::cosine_correlation(
                        &vec_of(f[1]),
                        &s1,
                        &vec_of(f[3]),
                        &s2,
                    )),
                    f64_field(isotopes::cosine_correlation_1d(
                        &vec_of(f[1]),
                        &vec_of(f[3]),
                    )),
                ]
            }
            "ppm" => vec![f64_field(hills::calculate_ppm(
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
            ))],
            "split" => {
                let trace = vec_of(f[3]);
                let hill = Hill {
                    scan_indices: (0..trace.len()).map(|i| i + 3).collect(),
                    peak_indices: (0..trace.len()).collect(),
                    mz_values: (0..trace.len()).map(|i| 500.0 + i as f64 * 1e-5).collect(),
                    intensities: trace.clone(),
                    rt_values: (0..trace.len()).map(|i| i as f64 * 1.5).collect(),
                    drift_times: vec![-1.0; trace.len()],
                    ion_mobilities: vec![-1.0; trace.len()],
                    length: trace.len(),
                    hill_idx: 7,
                    ..Hill::default()
                };
                let other = Hill {
                    hill_idx: 11,
                    ..hill.clone()
                };
                let split = hills::split_hills(
                    &[hill, other],
                    f[1].parse().unwrap(),
                    f[2].parse().unwrap(),
                )
                .unwrap();
                let mut text = String::new();
                for s in &split {
                    text.push_str(&format!(
                        "{}:{}:{}:{}:{}:{}:{}:{}:{}:{};",
                        s.hill_idx,
                        s.length,
                        s.scan_indices.first().copied().unwrap_or(0),
                        f64_field(s.mz_weighted_mean),
                        f64_field(s.rt_start),
                        f64_field(s.rt_end),
                        f64_field(s.rt_apex),
                        f64_field(s.intensity_apex),
                        f64_field(s.intensity_sum),
                        hill_digest(s)
                    ));
                }
                vec![split.len().to_string(), text]
            }
            "isovalley" => {
                let trace = vec_of(f[2]);
                let hills: Vec<Hill> = trace
                    .iter()
                    .enumerate()
                    .map(|(i, &apex)| Hill {
                        hill_idx: 100 + i,
                        intensity_apex: apex,
                        ..Hill::default()
                    })
                    .collect();
                let iso: Vec<IsotopeCandidate> = (0..trace.len())
                    .map(|i| IsotopeCandidate {
                        hill_idx: 100 + trace.len() - 1 - i,
                        isotope_number: i,
                        ..Default::default()
                    })
                    .collect();
                vec![
                    isotopes::check_isotope_valley_split(&iso, &hills, f[1].parse().unwrap())
                        .unwrap()
                        .to_string(),
                ]
            }
            other => panic!("unknown probe {other}"),
        };
        let expected: Vec<String> = expected
            .iter()
            .filter(|field| !field.is_empty() || f[0] == "split")
            .cloned()
            .collect();
        let got: Vec<String> = got
            .into_iter()
            .filter(|field| !field.is_empty() || f[0] == "split")
            .collect();
        assert_eq!(got, expected, "helpers.tsv line {lineno}: {}", f[0]);
        checked += 1;
    }
    assert_eq!(checked, 378);
}

#[test]
fn averagine_table_is_the_binomial_model() {
    // Tier 4: the embedded Release-build table against an independent
    // evaluation of the source's model with statrs.
    use statrs::distribution::{Binomial, Discrete};
    let table: Vec<(usize, u32, usize, Vec<f64>)> = fixture("helpers_oracle.tsv")
        .lines()
        .filter(|line| line.starts_with("table\t"))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].split(',')
                    .filter(|v| !v.is_empty())
                    .map(bits)
                    .collect(),
            )
        })
        .collect();
    assert_eq!(table.len(), 200);
    for (bin, n_c, _, probs) in &table {
        let (embedded_n_c, embedded) = isotopes::averagine_bin(*bin).unwrap();
        assert_eq!(embedded_n_c, *n_c);
        assert_eq!(embedded.as_slice(), probs.as_slice(), "bin {bin}");
        let expected_n_c = ((*bin as f64 * 100.0) / 111.1254 * 4.9384).round().max(1.0) as u64;
        assert_eq!(u64::from(*n_c), expected_n_c);
        let binomial = Binomial::new(0.0107, expected_n_c).unwrap();
        let raw: Vec<f64> = (0..10u64)
            .map(|k| {
                if k <= expected_n_c {
                    binomial.pmf(k)
                } else {
                    0.0
                }
            })
            .collect();
        let sum: f64 = raw.iter().sum();
        for (k, (&p, &r)) in probs.iter().zip(&raw).enumerate() {
            let independent = r / sum;
            assert!(
                (p - independent).abs() <= 1e-12 * independent.max(1e-300) + 1e-300,
                "bin {bin} isotope {k}: {p} vs {independent}"
            );
        }
    }
}

#[test]
fn defaults_match_release_build() {
    let defaults = Biosaur2Algorithm::defaults().unwrap();
    let mut names = Vec::new();
    for line in fixture("helpers_oracle.tsv")
        .lines()
        .filter(|l| l.starts_with("default\t"))
    {
        let f: Vec<&str> = line.split('\t').collect();
        let key = f[2];
        if names.contains(&key) {
            continue; // the Param iterator reports the last entry twice
        }
        names.push(key);
        let entry = defaults.entry(key).unwrap();
        let (type_code, value) = match &entry.value {
            ParamValue::String(v) => ("0", v.clone()),
            ParamValue::Integer(v) => ("1", v.to_string()),
            ParamValue::Float(v) => ("2", format!("{v:?}")),
            other => panic!("{key}: {other:?}"),
        };
        assert_eq!(type_code, f[3], "{key}: value type");
        assert_eq!(value, f[4], "{key}: value");
        assert_eq!(entry.min_float, bits(f[5]), "{key}: min_float");
        assert_eq!(entry.max_float, bits(f[6]), "{key}: max_float");
        assert_eq!(entry.min_int.to_string(), f[7], "{key}: min_int");
        assert_eq!(entry.max_int.to_string(), f[8], "{key}: max_int");
        let valid: String = entry
            .valid_strings
            .iter()
            .map(|s| format!("{s},"))
            .collect();
        assert_eq!(valid, f[9], "{key}: valid strings");
        let tags: String = entry.tags.iter().map(|s| format!("{s},")).collect();
        assert_eq!(tags, f[10], "{key}: tags");
        assert_eq!(entry.description, f[11], "{key}: description");
    }
    let order: Vec<String> = defaults
        .root()
        .entries
        .iter()
        .map(|e| e.name.clone())
        .collect();
    assert_eq!(order, names, "declaration order");
    assert_eq!(
        Biosaur2Algorithm::from_param(&defaults).unwrap(),
        Biosaur2Algorithm::new()
    );
    let round_trip = configured("htol=4.5;minlh=5;nm=true;convex_hulls=mass_traces;iuse=-1");
    let param = round_trip.to_param().unwrap();
    assert_eq!(Biosaur2Algorithm::from_param(&param).unwrap(), round_trip);
    assert_eq!(round_trip.convex_hulls, ConvexHullMode::MassTraces);
}

#[test]
fn parameter_restrictions_are_refused() {
    // The source's checkDefaults: out-of-range values throw InvalidParameter.
    for bad in [
        "hvf=0.5",
        "ivf=0.9",
        "mini=-1.0",
        "minlh=0",
        "cmin=0",
        "iuse=-2",
        "tof=yes",
        "convex_hulls=box",
    ] {
        let mut param = Biosaur2Algorithm::defaults().unwrap();
        let (key, value) = bad.split_once('=').unwrap();
        let typed = match param.value(key).unwrap() {
            ParamValue::Float(_) => ParamValue::Float(value.parse().unwrap()),
            ParamValue::Integer(_) => ParamValue::Integer(value.parse().unwrap()),
            _ => ParamValue::String(value.to_owned()),
        };
        param.set_value(key, typed, "", &[]).unwrap();
        assert!(Biosaur2Algorithm::from_param(&param).is_err(), "{bad}");
    }
    // Native: the public fields are checked by run().
    let mut algorithm = Biosaur2Algorithm::new();
    algorithm.hvf = 0.5;
    algorithm.set_ms_data(read_exp(&fixture("ct_run_basic.exp.gz")));
    let before = algorithm.ms_data().clone();
    assert!(algorithm.run(&mut UniqueIdGenerator::from_seed(1)).is_err());
    assert_eq!(algorithm.ms_data(), &before);
    algorithm.hvf = 1.3;
    algorithm.htol = f64::NAN;
    assert!(algorithm.run(&mut UniqueIdGenerator::from_seed(1)).is_err());
}

// ----- Biosaur2Algorithm_test.cpp (core bc9cc12), tier 3 -----

fn class_test_run_input() -> MSExperiment {
    // Biosaur2Algorithm_test.cpp:59-92
    let mut experiment = MSExperiment::new();
    for i in 0..10 {
        let mut spectrum = MSSpectrum {
            rt: f64::from(i),
            ms_level: 1,
            ..MSSpectrum::default()
        };
        spectrum.peaks = vec![
            Peak1D::new(1000.0, 10000.0),
            Peak1D::new(1001.003355, 6000.0),
            Peak1D::new(1000.993355, 100.0),
        ];
        experiment.spectra.push(spectrum);
    }
    experiment
}

#[test]
fn class_test_constructor() {
    // Biosaur2Algorithm_test.cpp:35-46
    let algorithm = Biosaur2Algorithm::new();
    assert!(algorithm.ms_data().is_empty());
}

#[test]
fn class_test_run_feature_map() {
    // Biosaur2Algorithm_test.cpp:48-121: itol 20, minmz 900, maxmz 1100.
    let mut algorithm = configured("itol=20.0;minmz=900.0;maxmz=1100.0");
    algorithm.set_ms_data(class_test_run_input());
    let map = algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    assert_eq!(map.len(), 1); // :106
    let feature = &map.features[0];
    assert!((feature.mz - 1000.0).abs() <= 1e-5 * 1000.0); // :111
    assert_eq!(feature.charge, 1); // :112
    assert_eq!(
        feature.metadata["n_isotopes"].data,
        MetaValueData::Integer(2)
    ); // :118
}

#[test]
fn class_test_set_ms_data() {
    // Biosaur2Algorithm_test.cpp:123-148: both overloads store the data.
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(MSSpectrum::default());
    let mut algorithm = Biosaur2Algorithm::new();
    algorithm.set_ms_data(experiment.clone());
    assert_eq!(algorithm.ms_data().len(), 1); // :132
    let mut moved = experiment;
    algorithm.set_ms_data(std::mem::take(&mut moved));
    assert_eq!(algorithm.ms_data().len(), 1); // :145
    assert_eq!(moved.len(), 0); // :146
}

#[test]
fn class_test_get_ms_data() {
    // Biosaur2Algorithm_test.cpp:150-170
    let mut algorithm = Biosaur2Algorithm::new();
    assert_eq!(algorithm.ms_data_mut().len(), 0); // :154
    algorithm.ms_data_mut().spectra.push(MSSpectrum::default());
    assert_eq!(algorithm.ms_data().len(), 1); // :160
    assert_eq!(Biosaur2Algorithm::new().ms_data().len(), 0); // :168
}

#[test]
fn class_test_run_with_hills() {
    // Biosaur2Algorithm_test.cpp:172-208
    let mut algorithm = configured("minmz=400.0;maxmz=600.0;mini=100.0");
    let mut experiment = MSExperiment::new();
    for i in 0..5 {
        let mut spectrum = MSSpectrum {
            rt: f64::from(i),
            ms_level: 1,
            ..MSSpectrum::default()
        };
        spectrum.peaks.push(Peak1D::new(500.0, 1000.0));
        experiment.spectra.push(spectrum);
    }
    algorithm.set_ms_data(experiment);
    let output = algorithm
        .run_with_details(&mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    assert!(!output.hills.is_empty()); // :204
}

fn class_test_faims_input(profile: bool) -> MSExperiment {
    // Biosaur2Algorithm_test.cpp:217-255 and :303-332
    let mut experiment = MSExperiment::new();
    let cvs = [-45.0, -65.0];
    for i in 0..10 {
        for (g, &cv) in cvs.iter().enumerate() {
            let offset = g as f64;
            let mut spectrum = MSSpectrum {
                rt: f64::from(i) * 2.0 + offset,
                ms_level: 1,
                drift_time: cv,
                drift_time_unit: DriftTimeUnit::FaimsCompensationVoltage,
                ..MSSpectrum::default()
            };
            if profile {
                spectrum.spectrum_type = SpectrumType::Profile;
                for k in -3..=3 {
                    let kf = f64::from(k);
                    spectrum.peaks.push(Peak1D::new(
                        500.0 + offset + kf * 0.01,
                        10000.0f32 * ((-0.5 * kf * kf).exp() as f32),
                    ));
                }
            } else {
                spectrum.peaks.push(Peak1D::new(500.0 + offset, 10000.0));
                spectrum
                    .peaks
                    .push(Peak1D::new(500.0 + offset + 1.003355, 5000.0));
            }
            experiment.spectra.push(spectrum);
        }
    }
    if !profile {
        let mut ms2 = MSSpectrum {
            rt: 1.5,
            ms_level: 2,
            drift_time: cvs[0],
            drift_time_unit: DriftTimeUnit::FaimsCompensationVoltage,
            ..MSSpectrum::default()
        };
        ms2.peaks.push(Peak1D::new(300.0, 100.0));
        experiment.spectra.push(ms2);
    }
    experiment.chromatograms.push(MSChromatogram {
        native_id: "TIC".into(),
        ..MSChromatogram::default()
    });
    experiment
}

#[test]
fn class_test_faims_path_preserves_ms_data() {
    // Biosaur2Algorithm_test.cpp:210-296
    let experiment = class_test_faims_input(false);
    assert_eq!(
        FaimsHelper::get_compensation_voltages(&experiment)
            .unwrap()
            .voltages
            .len(),
        2
    ); // :257
    let mut algorithm = configured("minmz=400.0;maxmz=600.0;mini=100.0");
    algorithm.set_ms_data(experiment);
    algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    let back = algorithm.ms_data();
    assert_eq!(back.len(), 20); // :273
    assert_eq!(back.chromatograms.len(), 1); // :274
    assert_eq!(
        FaimsHelper::get_compensation_voltages(back)
            .unwrap()
            .voltages
            .len(),
        2
    ); // :278
    assert!(back.spectra.windows(2).all(|w| w[0].rt <= w[1].rt)); // :279-284
    assert_eq!(back.spectra[0].rt, 0.0); // :285
    assert_eq!(back.spectra[back.len() - 1].rt, 19.0); // :286
    let ranges = back.combined_range_manager().unwrap();
    assert_eq!(ranges.min(MSDim::Rt).unwrap(), 0.0); // :292
    assert_eq!(ranges.max(MSDim::Rt).unwrap(), 19.0); // :293
    assert!((ranges.max(MSDim::Mz).unwrap() - 502.003355).abs() <= 1e-5 * 502.003355); // :294
}

#[test]
fn class_test_faims_profile_path_keeps_chromatograms_and_settings() {
    // Biosaur2Algorithm_test.cpp:298-353
    let mut experiment = class_test_faims_input(true);
    experiment.settings.date_time =
        openms::data_structures::DateTime::parse("2019-09-07T09:40:04").unwrap();
    let mut algorithm = configured("profile=true;minmz=400.0;maxmz=600.0;mini=100.0");
    algorithm.set_ms_data(experiment);
    algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    let back = algorithm.ms_data();
    assert_eq!(back.len(), 20); // :347
    assert_eq!(back.chromatograms.len(), 1); // :348
    assert_eq!(back.settings.date_time.iso_string(), "2019-09-07T09:40:04"); // :349
    assert_eq!(
        FaimsHelper::get_compensation_voltages(back)
            .unwrap()
            .voltages
            .len(),
        2
    ); // :351
}

#[test]
fn class_test_profile_path_refreshes_ranges() {
    // Biosaur2Algorithm_test.cpp:355-408: ranges are computed on demand here,
    // so they always describe the returned spectra.
    let mut experiment = MSExperiment::new();
    for i in 0..10 {
        let mut spectrum = MSSpectrum {
            rt: f64::from(i),
            ms_level: 1,
            spectrum_type: SpectrumType::Profile,
            ..MSSpectrum::default()
        };
        for k in -3..=3 {
            let kf = f64::from(k);
            spectrum.peaks.push(Peak1D::new(
                500.0 + kf * 0.01,
                10000.0f32 * ((-0.5 * kf * kf).exp() as f32),
            ));
        }
        experiment.spectra.push(spectrum);
    }
    let mut algorithm = configured("profile=true;minmz=400.0;maxmz=600.0;mini=100.0");
    algorithm.set_ms_data(experiment);
    algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    let back = algorithm.ms_data();
    assert_eq!(back.len(), 10); // :393
    assert!(
        FaimsHelper::get_compensation_voltages(back)
            .unwrap()
            .voltages
            .is_empty()
    ); // :394
    let ranges = back.combined_range_manager().unwrap();
    assert_eq!(ranges.min(MSDim::Rt).unwrap(), 0.0); // :395
    assert_eq!(ranges.max(MSDim::Rt).unwrap(), 9.0); // :396
}

#[test]
fn class_test_non_faims_path_preserves_ms_data() {
    // Biosaur2Algorithm_test.cpp:410-443
    let mut experiment = MSExperiment::new();
    for i in 0..10 {
        let mut spectrum = MSSpectrum {
            rt: f64::from(i),
            ms_level: 1,
            ..MSSpectrum::default()
        };
        spectrum.peaks.push(Peak1D::new(500.0, 10000.0));
        experiment.spectra.push(spectrum);
    }
    assert!(
        FaimsHelper::get_compensation_voltages(&experiment)
            .unwrap()
            .voltages
            .is_empty()
    ); // :426
    let mut algorithm = configured("minmz=400.0;maxmz=600.0;mini=100.0");
    algorithm.set_ms_data(experiment);
    algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    assert_eq!(algorithm.ms_data().len(), 10); // :438
    let ranges = algorithm.ms_data().combined_range_manager().unwrap();
    assert_eq!(ranges.min(MSDim::Rt).unwrap(), 0.0); // :440
    assert_eq!(ranges.max(MSDim::Rt).unwrap(), 9.0); // :441
}

// ----- native (tier 4) -----

#[test]
fn faims_groups_are_bit_identical_across_thread_counts() {
    let (input, parameters) = &cases()["syn_faims"];
    let experiment = load_input(input);
    let mut runs = Vec::new();
    for threads in [Threads::serial(), Threads::from_cli(4), Threads::all()] {
        let mut algorithm = configured(parameters);
        algorithm.threads = threads;
        algorithm.set_ms_data(experiment.clone());
        let output = algorithm
            .run_with_details(&mut UniqueIdGenerator::from_seed(3))
            .unwrap();
        runs.push((output, algorithm.ms_data().clone()));
    }
    assert_eq!(runs[0], runs[1]);
    assert_eq!(runs[0], runs[2]);
}

#[test]
fn a_failed_run_leaves_the_stored_data_unchanged() {
    // An ion-mobility array shorter than the peak list fails in the hill
    // linking, after profile/TOF preprocessing would have run.
    let mut experiment = read_exp(&fixture("syn_im1.exp.gz"));
    experiment.spectra[3].float_data_arrays[0].data.truncate(2);
    let mut algorithm = configured("tof=true");
    algorithm.set_ms_data(experiment.clone());
    let error = algorithm
        .run(&mut UniqueIdGenerator::from_seed(1))
        .unwrap_err();
    assert!(
        matches!(&error, Error::InvalidValue(m) if m.contains("Ion mobility array shorter than peak list.")),
        "{error}"
    );
    assert_eq!(algorithm.ms_data(), &experiment);
}

#[test]
fn non_finite_peaks_are_refused() {
    let mut experiment = class_test_run_input();
    experiment.spectra[4].peaks[1].intensity = f32::NAN;
    let mut algorithm = configured("itol=20.0;minmz=900.0;maxmz=1100.0");
    algorithm.set_ms_data(experiment.clone());
    assert!(algorithm.run(&mut UniqueIdGenerator::from_seed(1)).is_err());
    assert_eq!(algorithm.ms_data(), &experiment);
    // A non-finite peak outside the window is filtered as in the source.
    let mut experiment = class_test_run_input();
    experiment.spectra[4]
        .peaks
        .push(Peak1D::new(2000.0, f32::INFINITY));
    let mut algorithm = configured("itol=20.0;minmz=900.0;maxmz=1100.0");
    algorithm.set_ms_data(experiment);
    assert_eq!(
        algorithm
            .run(&mut UniqueIdGenerator::from_seed(1))
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn the_isotope_combination_ceiling_is_preflighted() {
    let (input, parameters) = &cases()["syn_im"];
    let mut algorithm = configured(parameters);
    algorithm.limits.max_isotope_combinations = 1;
    let experiment = load_input(input);
    algorithm.set_ms_data(experiment.clone());
    let error = algorithm
        .run(&mut UniqueIdGenerator::from_seed(1))
        .unwrap_err();
    assert!(
        error.to_string().contains("isotope combinations"),
        "{error}"
    );
    assert_eq!(algorithm.ms_data(), &experiment);
}

#[test]
fn tsv_writers_write_files() {
    let (input, parameters) = &cases()["ffc1_topp"];
    let mut algorithm = configured(parameters);
    algorithm.set_ms_data(load_input(input));
    let output = algorithm
        .run_with_details(&mut UniqueIdGenerator::from_seed(1))
        .unwrap();
    let dir = std::env::temp_dir().join(format!("biosaur2_tsv_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let features = dir.join("features.tsv");
    let hills_path = dir.join("hills.tsv");
    let message = output::write_tsv(&output.peptide_features, &features).unwrap();
    assert_eq!(
        message.text,
        format!(
            "Wrote {} features to TSV file: {}",
            output.peptide_features.len(),
            features.display()
        )
    );
    output::write_hills(&output.hills, &hills_path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&features).unwrap(),
        output::tsv_text(&output.peptide_features)
    );
    assert_eq!(
        std::fs::read_to_string(&hills_path).unwrap(),
        output::hills_text(&output.hills)
    );
    assert!(output::write_tsv(&output.peptide_features, &dir.join("missing/x.tsv")).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn feature_map_is_returned_by_the_simple_run() {
    let mut algorithm = configured("itol=20.0;minmz=900.0;maxmz=1100.0");
    algorithm.set_ms_data(class_test_run_input());
    let map: FeatureMap = algorithm.run(&mut UniqueIdGenerator::from_seed(1)).unwrap();
    assert_eq!(map.protein_identifications.len(), 1);
}

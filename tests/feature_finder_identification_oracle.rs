// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Executed differential (evidence tier 1) for `FeatureFinderIdentificationAlgorithm`
//! and `FeatureFinderAlgorithmMetaboIdent` against the C++ Release build.
//!
//! The fixtures under `tests/data/feature_finder_identification/oracle/` are the
//! output of the driver `../oracle/FF3-ffid/driver.cpp`, built with
//! `x86_64-conda-linux-gnu-g++` `-O3 -DNDEBUG` against the Linux x86-64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
//! ibminode06 with byte-identical output
//! (`tests/data/feature_finder_identification_provenance.json`). Each record is
//! `key<TAB>value`, doubles as `%.17g`; this test recomputes every record with
//! the port (the driver's construction of inputs is mirrored here statement
//! for statement) and compares keys and values. Numbers compare as `f64`,
//! exactly; strings compare as strings.
//!
//! The cases beyond the upstream TOPP tests exercise what those leave
//! uncovered: seeds, offset peptides, decoys, the score cutoff, the `CalcMass`
//! shift, FAIMS data, ion mobility, mapping tolerance, a minimum isotope
//! probability, adducts, isotope distributions given in the table, RT ranges,
//! overlapping targets and the targets the table contract skips.
#![cfg(all(feature = "mzml", feature = "featurexml", feature = "idxml"))]

#[path = "support/ffid_inputs.rs"]
mod ffid_inputs;

use openms::analysis::featurefinder::feature_finder_identification::FeatureFinderIdentificationAlgorithm;
use openms::analysis::featurefinder::feature_finder_metabo_ident::{
    FeatureFinderAlgorithmMetaboIdent, FeatureFinderMetaboIdentCompound,
};
use openms::analysis::openswath::chromatogram_extractor_algorithm::{
    ChromatogramExtractorAlgorithm, ExtractionCoordinates,
};
use openms::analysis::openswath::data_structures::Chromatogram;
use openms::analysis::openswath::spectrum_access::SpectrumAccess;
use openms::analysis::openswath::spectrum_access_openms::SimpleOpenMSSpectraFactory;
use openms::analysis::openswath::targeted_experiment::TargetedExperiment;
use openms::chemistry::{AASequence, EmpiricalFormula, PROTON_MASS_U};
use openms::concept::UniqueIdGenerator;
use openms::format::file_handler::FileHandler;
use openms::format::file_types::FileType;
use openms::format::peak_options::PeakFileOptions;
use openms::identification::{PeptideIdentification, ProteinIdentification};
use openms::kernel::{DataArray, Feature, FeatureMap, MSExperiment};
use openms::metadata::{DriftTimeUnit, MetaInfo, MetaValue, MetaValueData};
use openms::param::{Param, ParamValue};
use openms::system::file::TempDir;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;

const TEST_SEED: u64 = 19_991_231_235_959;

fn data(name: &str) -> PathBuf {
    ffid_inputs::data(name)
}

// ---------------------------------------------------------------- records

enum Val {
    F(f64),
    S(String),
}

#[derive(Default)]
struct Out {
    records: BTreeMap<String, Val>,
}

impl Out {
    fn p(&mut self, key: &str, value: f64) {
        self.records.insert(key.to_owned(), Val::F(value));
    }
    fn pi(&mut self, key: &str, value: i64) {
        self.p(key, value as f64);
    }
    fn ps(&mut self, key: &str, value: &str) {
        self.records
            .insert(key.to_owned(), Val::S(value.to_owned()));
    }
}

fn dump_meta(out: &mut Out, pre: &str, meta: &MetaInfo) {
    for (name, value) in meta {
        let key = format!("{pre}.{name}");
        match value.data() {
            MetaValueData::Float(f) => out.p(&key, *f),
            MetaValueData::Integer(i) => out.pi(&key, *i),
            MetaValueData::String(s) => out.ps(&key, s),
            MetaValueData::StringList(list) => {
                out.pi(&format!("{key}.n"), list.len() as i64);
                for (i, s) in list.iter().enumerate() {
                    out.ps(&format!("{key}.{i}"), s);
                }
            }
            MetaValueData::FloatList(list) => {
                out.pi(&format!("{key}.n"), list.len() as i64);
                for (i, v) in list.iter().enumerate() {
                    out.p(&format!("{key}.{i}"), *v);
                }
            }
            MetaValueData::IntegerList(list) => {
                out.pi(&format!("{key}.n"), list.len() as i64);
                for (i, v) in list.iter().enumerate() {
                    out.pi(&format!("{key}.{i}"), *v);
                }
            }
            MetaValueData::Empty => out.ps(&key, "<empty>"),
        }
    }
}

fn dump_peptide_id(out: &mut Out, p: &str, id: &PeptideIdentification) {
    out.ps(&format!("{p}.identifier"), &id.identifier);
    out.p(&format!("{p}.rt"), id.rt.unwrap_or(f64::NAN));
    out.p(&format!("{p}.mz"), id.mz.unwrap_or(f64::NAN));
    out.ps(&format!("{p}.score_type"), &id.score_type);
    out.pi(
        &format!("{p}.higher_better"),
        i64::from(id.higher_score_better),
    );
    dump_meta(out, &format!("{p}.meta"), &id.metadata);
    out.pi(&format!("{p}.nhits"), id.hits.len() as i64);
    for (j, hit) in id.hits.iter().enumerate() {
        let q = format!("{p}.hit{j}");
        out.p(&format!("{q}.score"), hit.score);
        out.pi(&format!("{q}.charge"), i64::from(hit.charge));
        out.ps(&format!("{q}.sequence"), &hit.sequence.to_string());
        dump_meta(out, &format!("{q}.meta"), &hit.metadata);
        out.pi(&format!("{q}.nev"), hit.evidences.len() as i64);
        for (k, evidence) in hit.evidences.iter().enumerate() {
            out.ps(&format!("{q}.ev{k}"), &evidence.protein_accession);
        }
    }
}

fn dump_feature(out: &mut Out, pre: &str, f: &Feature, top: bool) {
    out.p(&format!("{pre}.rt"), f.rt);
    out.p(&format!("{pre}.mz"), f.mz);
    out.p(&format!("{pre}.intensity"), f64::from(f.intensity));
    out.p(&format!("{pre}.quality"), f64::from(f.quality));
    out.p(&format!("{pre}.quality_rt"), f64::from(f.quality_rt));
    out.p(&format!("{pre}.quality_mz"), f64::from(f.quality_mz));
    out.p(&format!("{pre}.width"), f64::from(f.width));
    out.pi(&format!("{pre}.charge"), i64::from(f.charge));
    out.pi(&format!("{pre}.nhulls"), f.convex_hulls.len() as i64);
    for (h, hull) in f.convex_hulls.iter().enumerate() {
        let points = hull.hull_points();
        out.pi(&format!("{pre}.hull{h}.n"), points.len() as i64);
        for (p, point) in points.iter().enumerate() {
            out.p(&format!("{pre}.hull{h}.p{p}.rt"), point.rt);
            out.p(&format!("{pre}.hull{h}.p{p}.mz"), point.mz);
        }
    }
    dump_meta(out, &format!("{pre}.meta"), &f.metadata);
    out.pi(
        &format!("{pre}.npepid"),
        f.peptide_identifications.len() as i64,
    );
    for (i, id) in f.peptide_identifications.iter().enumerate() {
        dump_peptide_id(out, &format!("{pre}.pepid{i}"), id);
    }
    out.pi(&format!("{pre}.nsub"), f.subordinates.len() as i64);
    if top {
        for (i, sub) in f.subordinates.iter().enumerate() {
            dump_feature(out, &format!("{pre}.sub{i}"), sub, false);
        }
    }
}

fn dump_map(out: &mut Out, pre: &str, fm: &FeatureMap) {
    out.pi(&format!("{pre}.n"), fm.features.len() as i64);
    for (i, f) in fm.features.iter().enumerate() {
        dump_feature(out, &format!("{pre}.f{i}"), f, true);
    }
    let un = &fm.unassigned_peptide_identifications;
    out.pi(&format!("{pre}.nunassigned"), un.len() as i64);
    for (i, id) in un.iter().enumerate() {
        dump_peptide_id(out, &format!("{pre}.un{i}"), id);
    }
    let prot = &fm.protein_identifications;
    out.pi(&format!("{pre}.nprot"), prot.len() as i64);
    for (i, p) in prot.iter().enumerate() {
        out.ps(&format!("{pre}.prot{i}.identifier"), &p.identifier);
        out.pi(&format!("{pre}.prot{i}.nhits"), p.hits.len() as i64);
        for (h, hit) in p.hits.iter().enumerate() {
            out.ps(&format!("{pre}.prot{i}.hit{h}"), &hit.accession);
        }
    }
    let run = fm.primary_ms_run_path().unwrap();
    out.pi(&format!("{pre}.nrun"), run.len() as i64);
    for (i, r) in run.iter().enumerate() {
        out.ps(&format!("{pre}.run{i}"), r);
    }
    out.ps(
        &format!("{pre}.has_unique_id"),
        if fm.unique_id != 0 { "true" } else { "false" },
    );
}

fn dump_uids(out: &mut Out, pre: &str, fm: &FeatureMap) {
    out.pi(
        &format!("{pre}.map_uid_valid"),
        i64::from(fm.unique_id != 0),
    );
    for (i, f) in fm.features.iter().enumerate() {
        out.ps(&format!("{pre}.f{i}.uid"), &f.unique_id.to_string());
    }
}

fn dump_library(out: &mut Out, pre: &str, lib: &TargetedExperiment) {
    out.pi(&format!("{pre}.npeptides"), lib.peptides().len() as i64);
    for (i, p) in lib.peptides().iter().enumerate() {
        let q = format!("{pre}.pep{i}");
        out.ps(&format!("{q}.id"), &p.id);
        out.ps(&format!("{q}.sequence"), &p.sequence);
        out.ps(&format!("{q}.group"), p.peptide_group_label());
        out.pi(
            &format!("{q}.charge"),
            i64::from(p.charge_state().unwrap_or(0)),
        );
        out.p(&format!("{q}.drift"), p.drift_time());
        out.pi(&format!("{q}.nrts"), p.rts.len() as i64);
        for (r, rt) in p.rts.iter().enumerate() {
            out.p(&format!("{q}.rt{r}"), rt.rt().unwrap_or(0.0));
            out.pi(
                &format!("{q}.rt{r}.type"),
                i64::from(rt.retention_time_type as i8),
            );
            out.pi(
                &format!("{q}.rt{r}.unit"),
                i64::from(rt.retention_time_unit as i8),
            );
        }
        out.pi(&format!("{q}.nprot"), p.protein_refs.len() as i64);
        for (r, prot) in p.protein_refs.iter().enumerate() {
            out.ps(&format!("{q}.prot{r}"), prot);
        }
    }
    out.pi(&format!("{pre}.nproteins"), lib.proteins().len() as i64);
    for (i, p) in lib.proteins().iter().enumerate() {
        out.ps(&format!("{pre}.protein{i}"), &p.id);
    }
    out.pi(&format!("{pre}.ncompounds"), lib.compounds().len() as i64);
    for (i, c) in lib.compounds().iter().enumerate() {
        let q = format!("{pre}.comp{i}");
        out.ps(&format!("{q}.id"), &c.id);
        out.ps(&format!("{q}.formula"), &c.molecular_formula);
        out.p(&format!("{q}.mass"), c.theoretical_mass);
        out.pi(
            &format!("{q}.charge"),
            i64::from(c.charge_state().unwrap_or(0)),
        );
        out.p(&format!("{q}.drift"), c.drift_time());
        dump_meta(out, &format!("{q}.meta"), &c.cv_terms.metadata);
        out.pi(&format!("{q}.nrts"), c.rts.len() as i64);
        for (r, rt) in c.rts.iter().enumerate() {
            out.p(&format!("{q}.rt{r}"), rt.rt().unwrap_or(0.0));
            out.pi(
                &format!("{q}.rt{r}.type"),
                i64::from(rt.retention_time_type as i8),
            );
            out.pi(
                &format!("{q}.rt{r}.unit"),
                i64::from(rt.retention_time_unit as i8),
            );
        }
    }
    out.pi(
        &format!("{pre}.ntransitions"),
        lib.transitions().len() as i64,
    );
    for (i, t) in lib.transitions().iter().enumerate() {
        let q = format!("{pre}.tr{i}");
        out.ps(&format!("{q}.native_id"), t.native_id());
        out.ps(&format!("{q}.peptide_ref"), t.peptide_ref());
        out.ps(&format!("{q}.compound_ref"), t.compound_ref());
        out.p(&format!("{q}.precursor_mz"), t.precursor_mz());
        out.p(&format!("{q}.product_mz"), t.product_mz());
        out.p(&format!("{q}.library_intensity"), t.library_intensity());
    }
}

// ---------------------------------------------------------------- comparing

fn oracle_records(case: &str) -> BTreeMap<String, String> {
    let path = data(&format!("oracle/{case}.stdout"));
    let text = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("oracle record file {}: {e}", path.display()));
    let mut records = BTreeMap::new();
    for line in text.lines() {
        if let Some((key, value)) = line.split_once('\t') {
            if !value.contains('\t') && key != "case" {
                records.insert(key.to_owned(), value.to_owned());
            }
        }
    }
    records
}

/// A double meta value as the source's writer prints it and its reader reads
/// it back (`Internal::NumericFormatting::appendNumeric`): fixed notation with
/// fifteen decimals for zero and magnitudes in [1e-2, 1e4), the shortest
/// round-trip text, that is the value itself, for every other magnitude.
fn as_written_meta(value: f64) -> f64 {
    let magnitude = value.abs();
    if value.is_finite() && (magnitude == 0.0 || (1e-2..1e4).contains(&magnitude)) {
        format!("{value:.15}").parse::<f64>().unwrap_or(value)
    } else {
        value
    }
}

/// The members of a feature that are `float`s in the source.
const FLOAT_MEMBERS: &[&str] = &[".quality", ".quality_rt", ".quality_mz", ".width", ".intensity"];

/// A `float` as the source's writer prints it and its reader reads it back:
/// fixed notation with six decimals for zero and magnitudes in [1e-2, 1e4),
/// the shortest round-trip text, that is the value itself, otherwise.
fn as_written_float(value: f32) -> f64 {
    let magnitude = value.abs();
    if value.is_finite() && (magnitude == 0.0 || (1e-2..1e4).contains(&magnitude)) {
        f64::from(format!("{value:.6}").parse::<f32>().unwrap_or(value))
    } else {
        f64::from(value)
    }
}

/// The candidates as the source's file holds them. The source writes the
/// candidates before the protein identifications are set, and its writer
/// leaves out a peptide identification whose run is not in the file, so a
/// candidate carries none; the in-memory candidates of the port keep them.
fn candidates_as_written(candidates: &FeatureMap) -> FeatureMap {
    let mut written = candidates.clone();
    for feature in &mut written.features {
        feature.peptide_identifications.clear();
        for sub in &mut feature.subordinates {
            sub.peptide_identifications.clear();
        }
    }
    written.unassigned_peptide_identifications.clear();
    written
}

fn oracle_stderr(case: &str) -> String {
    fs::read_to_string(data(&format!("oracle/{case}.stderr"))).unwrap()
}

/// Compare the produced records with the oracle's, key for key. Keys whose
/// name ends in one of `ignore_suffixes` are skipped on both sides.
fn compare(case: &str, produced: &Out, ignore_suffixes: &[&str]) {
    let skipped = |key: &str| ignore_suffixes.iter().any(|s| key.ends_with(s));
    let expected = oracle_records(case);
    let mut problems: Vec<String> = Vec::new();
    for (key, value) in &expected {
        if skipped(key) {
            continue;
        }
        match produced.records.get(key) {
            None => problems.push(format!("missing {key} (oracle {value})")),
            Some(Val::S(s)) => {
                if s != value {
                    problems.push(format!("{key}: '{s}' != oracle '{value}'"));
                }
            }
            Some(Val::F(x)) => match value.parse::<f64>() {
                Ok(e) => {
                    // The oracle's candidates went through the source's
                    // featureXML file; the port keeps the candidates in
                    // memory, so a double meta value is put through the text
                    // the source writes. Every other candidate value is
                    // written (and compared) in full.
                    let x = &if key.starts_with("cand.") && key.contains(".meta.") {
                        as_written_meta(*x)
                    } else if key.starts_with("cand.") && FLOAT_MEMBERS.iter().any(|m| key.ends_with(m))
                    {
                        as_written_float(*x as f32)
                    } else {
                        *x
                    };
                    let same = (x.is_nan() && e.is_nan()) || x.to_bits() == e.to_bits() || *x == e;
                    if !same {
                        problems.push(format!("{key}: {x:e} != oracle {e:e}"));
                    }
                }
                Err(_) => problems.push(format!("{key}: number {x} != oracle '{value}'")),
            },
        }
    }
    for key in produced.records.keys() {
        if !skipped(key) && !expected.contains_key(key) {
            problems.push(format!("unexpected {key}"));
        }
    }
    // the first differences of every record family, not only the first ones
    let mut shown: Vec<String> = Vec::new();
    let mut per_family: BTreeMap<String, usize> = BTreeMap::new();
    for problem in &problems {
        let family = problem
            .split('.')
            .take(2)
            .collect::<Vec<_>>()
            .join(".")
            .chars()
            .take(24)
            .collect::<String>();
        let count = per_family.entry(family).or_insert(0);
        *count += 1;
        if *count <= 6 {
            shown.push(problem.clone());
        }
    }
    assert!(
        problems.is_empty(),
        "{case}: {} difference(s) in {} family(ies), first of each:\n{}",
        problems.len(),
        per_family.len(),
        shown.join("\n")
    );
}

// ---------------------------------------------------------------- inputs

fn load_ms1(name: &str) -> MSExperiment {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    FileHandler::load_experiment_with_options(data(name), &[FileType::MzMl], &options).unwrap()
}

/// Alternate spectra get the compensation voltage -45 and -65 V.
fn add_faims(experiment: &mut MSExperiment) {
    for (i, spectrum) in experiment.spectra.iter_mut().enumerate() {
        spectrum.drift_time = if i % 2 == 0 { -45.0 } else { -65.0 };
        spectrum.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    }
}

fn frac(x: f64) -> f64 {
    x - x.floor()
}

/// A 1/K0 value per peak: 0.6 + 0.4 * frac(mz / 50).
fn im_of_mz(mz: f64) -> f64 {
    0.6 + 0.4 * frac(mz / 50.0)
}

fn add_per_peak_im(experiment: &mut MSExperiment) {
    for spectrum in &mut experiment.spectra {
        let values: Vec<f32> = spectrum
            .peaks
            .iter()
            .map(|peak| im_of_mz(peak.mz) as f32)
            .collect();
        spectrum.float_data_arrays.push(DataArray::new(
            "raw inverse reduced ion mobility array",
            values,
        ));
    }
}

fn seeds_from_fixture(case: &str) -> FeatureMap {
    FileHandler::load_feature_map(
        data(&format!("oracle/{case}.seeds.featureXML")),
        &[FileType::FeatureXml],
    )
    .unwrap()
}

/// `AASequence::getMZ(charge)` of the driver: the weight in the source's order.
fn source_mz(sequence: &AASequence, charge: i32) -> f64 {
    openms::analysis::featurefinder::feature_finder_identification::source_mz(sequence, charge)
        .unwrap()
}

fn meta_float(value: f64) -> MetaValue {
    MetaValue::try_from(value).unwrap()
}

#[derive(Default, Clone, Copy)]
struct FfidOptions {
    seeds: bool,
    faims: bool,
    faims_ids: bool,
    im: bool,
    decoys: bool,
    calcmass: bool,
}

fn float(value: f64) -> ParamValue {
    ParamValue::Float(value)
}
fn int(value: i64) -> ParamValue {
    ParamValue::Integer(value)
}
fn string(value: &str) -> ParamValue {
    ParamValue::String(value.to_owned())
}

fn ffid_base() -> Vec<(&'static str, ParamValue)> {
    vec![
        ("extract:mz_window", float(0.1)),
        ("detect:peak_width", float(60.0)),
        ("model:type", string("none")),
    ]
}

fn with(
    mut base: Vec<(&'static str, ParamValue)>,
    extra: &[(&'static str, ParamValue)],
) -> Vec<(&'static str, ParamValue)> {
    for (key, value) in extra {
        match base.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value.clone(),
            None => base.push((key, value.clone())),
        }
    }
    base
}

fn apply(params: &[(&'static str, ParamValue)], mut param: Param) -> Param {
    for (key, value) in params {
        param.set_value(key, value.clone(), "", &[]).unwrap();
    }
    param
}

fn run_ffid(case: &str, params: &[(&'static str, ParamValue)], options: FfidOptions) -> Out {
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let mut experiment = load_ms1("FeatureFinderIdentification_1_input.mzML.gz");
    if options.faims {
        add_faims(&mut experiment);
    }
    if options.im {
        add_per_peak_im(&mut experiment);
    }
    let document = ffid_inputs::load_identifications(&work);
    let proteins: Vec<ProteinIdentification> = document.protein_identifications;
    let mut peptides: Vec<PeptideIdentification> = document.peptide_identifications;
    for (j, id) in peptides.iter_mut().enumerate() {
        if id.hits.is_empty() {
            continue;
        }
        if options.faims_ids {
            if j % 3 == 0 {
                id.metadata.insert("FAIMS_CV".to_owned(), meta_float(-45.0));
            } else if j % 3 == 1 {
                id.metadata.insert("FAIMS_CV".to_owned(), meta_float(-65.0));
            }
        }
        if options.im && j % 4 != 3 {
            let mz = id.mz.unwrap_or(f64::NAN);
            id.metadata
                .insert("IM".to_owned(), meta_float(im_of_mz(mz)));
        }
        if options.decoys && j % 5 == 0 {
            id.hits[0]
                .metadata
                .insert("target_decoy".to_owned(), MetaValue::from("decoy"));
        }
        if options.calcmass {
            let hit = &mut id.hits[0];
            let z = hit.charge;
            if z != 0 && j % 7 == 0 {
                let value = source_mz(&hit.sequence, z) + 15.9949146 / f64::from(z);
                hit.metadata
                    .insert("CalcMass".to_owned(), meta_float(value));
            }
            if z != 0 && j % 7 == 3 {
                let value = source_mz(&hit.sequence, z) + 0.3 / f64::from(z);
                hit.metadata
                    .insert("CalcMass".to_owned(), meta_float(value));
            }
        }
    }
    let seeds = if options.seeds {
        seeds_from_fixture(case)
    } else {
        FeatureMap::new()
    };

    let mut param = apply(
        params,
        FeatureFinderIdentificationAlgorithm::defaults().unwrap(),
    );
    let candidates = work.path().join("candidates.featureXML");
    let wants_candidates = params.iter().any(|(k, v)| {
        *k == "candidates_out" && matches!(v, ParamValue::String(s) if s == "__out__")
    });
    if wants_candidates {
        param
            .set_value(
                "candidates_out",
                string(&candidates.to_string_lossy()),
                "",
                &[],
            )
            .unwrap();
    }
    let mut algorithm = FeatureFinderIdentificationAlgorithm::from_param(&param).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(TEST_SEED));
    algorithm.set_capture_candidates(wants_candidates);
    algorithm.set_ms_data(experiment);
    let mut features = FeatureMap::new();
    algorithm
        .run(
            peptides,
            &proteins,
            &mut features,
            &seeds,
            "FeatureFinderIdentification_1_input.mzML",
        )
        .unwrap();

    let mut out = Out::default();
    dump_map(&mut out, "fm", &features);
    dump_uids(&mut out, "fm", &features);
    dump_library(&mut out, "lib", algorithm.library());
    out.pi(
        "chrom.n",
        algorithm.chromatograms().chromatograms.len() as i64,
    );
    out.pi("msdata.n", algorithm.ms_data().spectra.len() as i64);
    if wants_candidates {
        // The candidates of the run itself. The file the run wrote is read back
        // too and must hold the same features; the native writer prints a
        // float meta value with the precision of a float (the writer's own
        // behaviour, recorded in the lane report), so the values compared with
        // the Release build are the in-memory ones.
        let kept = algorithm.candidates().expect("captured candidates");
        let loaded = ffid_inputs::load_feature_map(&candidates);
        assert_eq!(kept.features.len(), loaded.features.len());
        dump_map(&mut out, "cand", &candidates_as_written(kept));
    }
    out
}

// ---------------------------------------------------------------- MetaboIdent

fn list_of(text: &str) -> Vec<f64> {
    text.split(',')
        .filter(|part| !part.trim().is_empty())
        .map(|part| part.trim().parse::<f64>().unwrap())
        .collect()
}

fn read_table() -> Vec<FeatureFinderMetaboIdentCompound> {
    let content = fs::read_to_string(data("FeatureFinderMetaboIdent_1_input.tsv")).unwrap();
    let mut names = BTreeSet::new();
    let mut table = Vec::new();
    for line in content.lines().skip(1) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 7 || parts[0].is_empty() || !names.insert(parts[0].to_owned()) {
            continue;
        }
        table.push(FeatureFinderMetaboIdentCompound::new(
            parts[0],
            parts[1],
            parts[2].parse::<f64>().unwrap(),
            list_of(parts[3]).into_iter().map(|c| c as i32).collect(),
            list_of(parts[4]),
            list_of(parts[5]),
            list_of(parts[6]),
            Vec::new(),
            "",
        ));
    }
    table
}

fn variants_table(im: bool) -> Vec<FeatureFinderMetaboIdentCompound> {
    let base = read_table();
    let mut table = Vec::new();
    for (i, b) in base.iter().enumerate() {
        let mut charges = b.charges().to_vec();
        let mut rt_ranges = b.rt_ranges().to_vec();
        let mut iso = b.isotope_distribution().to_vec();
        let mut formula = b.formula().to_owned();
        let mut mass = b.mass();
        let mut adduct = String::new();
        let mut ims = Vec::new();
        if i < 3 {
            charges = vec![1, 2];
        }
        if i == 3 {
            adduct = "M+Na;1+".to_owned();
        }
        if i == 4 {
            formula = String::new();
            mass = 244.0695;
        }
        if i == 5 {
            iso = vec![1.0, 0.11, 0.02];
        }
        if i == 6 {
            rt_ranges = vec![30.0];
        }
        if i == 7 {
            charges = vec![1, -1];
        }
        if i == 2 {
            adduct = "M+K;1+".to_owned();
        }
        if im {
            // the 1/K0 the peaks of the [M+H]+ ion carry
            let m = if mass > 0.0 {
                mass
            } else {
                EmpiricalFormula::parse(&formula).unwrap().mono_mass()
            };
            ims = vec![im_of_mz(m + PROTON_MASS_U)];
        }
        table.push(FeatureFinderMetaboIdentCompound::new(
            b.name(),
            formula,
            mass,
            charges,
            b.rts().to_vec(),
            rt_ranges,
            iso,
            ims,
            adduct,
        ));
    }
    let compound = |name: &str,
                    formula: &str,
                    mass: f64,
                    charges: Vec<i32>,
                    rts: Vec<f64>,
                    rt_ranges: Vec<f64>,
                    ims: Vec<f64>,
                    adduct: &str| {
        FeatureFinderMetaboIdentCompound::new(
            name,
            formula,
            mass,
            charges,
            rts,
            rt_ranges,
            Vec::new(),
            ims,
            adduct,
        )
    };
    // violating the contract: skipped with a message
    table.push(compound(
        "no_mass_or_formula",
        "",
        0.0,
        vec![1],
        vec![100.0],
        vec![],
        vec![],
        "",
    ));
    table.push(compound(
        "no_rt",
        "C6H12O6",
        0.0,
        vec![1],
        vec![],
        vec![],
        vec![],
        "",
    ));
    table.push(compound(
        "zero_charge",
        "C6H12O6",
        0.0,
        vec![0, 1],
        vec![150.0],
        vec![],
        vec![],
        "",
    ));
    table.push(compound(
        "bad_rt_ranges",
        "C6H12O6",
        0.0,
        vec![1],
        vec![150.0, 160.0],
        vec![1.0, 2.0, 3.0],
        vec![],
        "",
    ));
    table.push(compound(
        "bad_im",
        "C6H12O6",
        0.0,
        vec![1],
        vec![150.0, 160.0],
        vec![],
        vec![0.7, 0.8, 0.9],
        "",
    ));
    table.push(compound(
        "bad_adduct",
        "C6H12O6",
        0.0,
        vec![1],
        vec![150.0],
        vec![],
        vec![],
        "M?H;1+",
    ));
    table.push(compound(
        "bad_adduct_mass",
        "",
        180.0634,
        vec![2],
        vec![150.0],
        vec![],
        vec![],
        "oops",
    ));
    table
}

fn mi_base() -> Vec<(&'static str, ParamValue)> {
    vec![
        ("extract:mz_window", float(5.0)),
        ("extract:rt_window", float(20.0)),
        ("detect:peak_width", float(3.0)),
    ]
}

fn run_metabo(
    case: &str,
    params: &[(&'static str, ParamValue)],
    table_kind: u8,
    faims: bool,
    im: bool,
) -> Out {
    let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let mut experiment = load_ms1("FeatureFinderMetaboIdent_1_input.mzML.gz");
    if faims {
        add_faims(&mut experiment);
    }
    if im {
        add_per_peak_im(&mut experiment);
    }
    let table = if table_kind == 0 {
        read_table()
    } else {
        variants_table(im)
    };
    let mut param = apply(
        params,
        FeatureFinderAlgorithmMetaboIdent::defaults().unwrap(),
    );
    let candidates = work.path().join("candidates.featureXML");
    let wants_candidates = params.iter().any(|(k, v)| {
        *k == "candidates_out" && matches!(v, ParamValue::String(s) if s == "__out__")
    });
    if wants_candidates {
        param
            .set_value(
                "candidates_out",
                string(&candidates.to_string_lossy()),
                "",
                &[],
            )
            .unwrap();
    }
    let mut algorithm = FeatureFinderAlgorithmMetaboIdent::from_param(&param).unwrap();
    algorithm.set_unique_id_generator(UniqueIdGenerator::from_seed(TEST_SEED));
    algorithm.set_capture_candidates(wants_candidates);
    algorithm.set_ms_data(experiment);
    let mut features = FeatureMap::new();
    algorithm
        .run(
            &table,
            &mut features,
            "FeatureFinderMetaboIdent_1_input.mzML",
        )
        .unwrap();

    let _ = case;
    let mut out = Out::default();
    dump_map(&mut out, "fm", &features);
    dump_uids(&mut out, "fm", &features);
    dump_library(&mut out, "lib", algorithm.library());
    out.pi(
        "chrom.n",
        algorithm.chromatograms().chromatograms.len() as i64,
    );
    out.pi("nshared", algorithm.n_shared() as i64);
    let points = algorithm.transformations().data_points();
    out.pi("trafo.n", points.len() as i64);
    for (i, point) in points.iter().enumerate() {
        out.p(&format!("trafo.p{i}.first"), point.x);
        out.p(&format!("trafo.p{i}.second"), point.y);
        out.ps(&format!("trafo.p{i}.note"), &point.note);
    }
    if wants_candidates {
        // The candidates of the run itself. The file the run wrote is read back
        // too and must hold the same features; the native writer prints a
        // float meta value with the precision of a float (the writer's own
        // behaviour, recorded in the lane report), so the values compared with
        // the Release build are the in-memory ones.
        let kept = algorithm.candidates().expect("captured candidates");
        let loaded = ffid_inputs::load_feature_map(&candidates);
        assert_eq!(kept.features.len(), loaded.features.len());
        dump_map(&mut out, "cand", &candidates_as_written(kept));
    }
    out
}

// ---------------------------------------------------------------- the cases

/// The unique IDs are compared in `unique_ids_follow_the_source_stream`.
const IGNORE: &[&str] = &[".uid"];

/// The identification run identifiers of the FeatureFinderIdentification cases
/// are generated by the idXML reader (`<engine>_<date>_<random>` in the source),
/// so only their count and structure are comparable.
const FFID_IGNORE: &[&str] = &[".uid", ".identifier"];

#[test]
fn ffid_none() {
    let out = run_ffid("ffid_none", &ffid_base(), FfidOptions::default());
    compare("ffid_none", &out, FFID_IGNORE);
}

#[test]
fn ffid_symmetric() {
    let params = with(ffid_base(), &[("model:type", string("symmetric"))]);
    let out = run_ffid("ffid_symmetric", &params, FfidOptions::default());
    compare("ffid_symmetric", &out, FFID_IGNORE);
}

#[test]
fn ffid_asymmetric() {
    let params = with(ffid_base(), &[("model:type", string("asymmetric"))]);
    let out = run_ffid("ffid_asymmetric", &params, FfidOptions::default());
    compare("ffid_asymmetric", &out, FFID_IGNORE);
}

#[test]
fn ffid_batch10_candidates() {
    let params = with(
        ffid_base(),
        &[
            ("extract:batch_size", int(10)),
            ("candidates_out", string("__out__")),
        ],
    );
    let out = run_ffid("ffid_batch10_candidates", &params, FfidOptions::default());
    compare("ffid_batch10_candidates", &out, FFID_IGNORE);
}

#[test]
fn ffid_mapping_tolerance() {
    let params = with(
        ffid_base(),
        &[
            ("detect:mapping_tolerance", float(0.3)),
            ("extract:rt_window", float(120.0)),
        ],
    );
    let out = run_ffid("ffid_mapping_tolerance", &params, FfidOptions::default());
    compare("ffid_mapping_tolerance", &out, FFID_IGNORE);
}

#[test]
fn ffid_pmin() {
    let params = with(ffid_base(), &[("extract:isotope_pmin", float(0.05))]);
    let out = run_ffid("ffid_pmin", &params, FfidOptions::default());
    compare("ffid_pmin", &out, FFID_IGNORE);
}

#[test]
fn ffid_ppm() {
    let params = with(
        ffid_base(),
        &[
            ("extract:mz_window", float(20.0)),
            ("extract:n_isotopes", int(3)),
        ],
    );
    let out = run_ffid("ffid_ppm", &params, FfidOptions::default());
    compare("ffid_ppm", &out, FFID_IGNORE);
}

#[test]
fn ffid_seeds() {
    let options = FfidOptions {
        seeds: true,
        ..FfidOptions::default()
    };
    let out = run_ffid("ffid_seeds", &ffid_base(), options);
    compare("ffid_seeds", &out, FFID_IGNORE);
}

#[test]
fn ffid_seeds_offset() {
    let options = FfidOptions {
        seeds: true,
        ..FfidOptions::default()
    };
    let params = with(ffid_base(), &[("add_mass_offset_peptides", float(5.0))]);
    let out = run_ffid("ffid_seeds_offset", &params, options);
    compare("ffid_seeds_offset", &out, FFID_IGNORE);
}

#[test]
fn ffid_seeds_apex_off() {
    let options = FfidOptions {
        seeds: true,
        ..FfidOptions::default()
    };
    let params = with(ffid_base(), &[("seed_apex_rt_tolerance", float(0.0))]);
    let out = run_ffid("ffid_seeds_apex_off", &params, options);
    compare("ffid_seeds_apex_off", &out, FFID_IGNORE);
}

#[test]
fn ffid_decoys() {
    let options = FfidOptions {
        decoys: true,
        ..FfidOptions::default()
    };
    let out = run_ffid("ffid_decoys", &ffid_base(), options);
    compare("ffid_decoys", &out, FFID_IGNORE);
}

#[test]
fn ffid_decoys_quantified() {
    let options = FfidOptions {
        decoys: true,
        ..FfidOptions::default()
    };
    let params = with(ffid_base(), &[("quantify_decoys", string("true"))]);
    let out = run_ffid("ffid_decoys_quantified", &params, options);
    compare("ffid_decoys_quantified", &out, FFID_IGNORE);
}

#[test]
fn ffid_psm_cutoff_string_is_read_as_the_number_it_spells() {
    // The Release build reads a string `min_psm_cutoff` through the `double`
    // member of the value union (the address of the string), so a cutoff of
    // "1e30" acts as a number near zero there and sets all but eight features
    // aside (oracle case `ffid_psm_cutoff_string`, the defect candidate of the
    // support document). The port reads the number: with q-values (lower is
    // better) a cutoff of 1e30 drops nothing, and the run is the one without a
    // cutoff.
    let defective = oracle_records("ffid_psm_cutoff_string");
    let plain = oracle_records("ffid_none");
    assert_eq!(defective["fm.n"], "8");
    assert_eq!(plain["fm.n"], "23");
    let params = with(ffid_base(), &[("min_psm_cutoff", string("1e30"))]);
    let out = run_ffid("ffid_none", &params, FfidOptions::default());
    compare("ffid_none", &out, FFID_IGNORE);
}

#[test]
fn ffid_calcmass() {
    let options = FfidOptions {
        calcmass: true,
        ..FfidOptions::default()
    };
    let out = run_ffid("ffid_calcmass", &ffid_base(), options);
    compare("ffid_calcmass", &out, FFID_IGNORE);
}

#[test]
fn ffid_faims() {
    let options = FfidOptions {
        faims: true,
        ..FfidOptions::default()
    };
    let out = run_ffid("ffid_faims", &ffid_base(), options);
    compare("ffid_faims", &out, FFID_IGNORE);
}

#[test]
fn ffid_faims_ids() {
    let options = FfidOptions {
        faims: true,
        faims_ids: true,
        ..FfidOptions::default()
    };
    let params = with(ffid_base(), &[("faims:merge_features", string("false"))]);
    let out = run_ffid("ffid_faims_ids", &params, options);
    compare("ffid_faims_ids", &out, FFID_IGNORE);
}

#[test]
fn ffid_im() {
    let options = FfidOptions {
        im: true,
        ..FfidOptions::default()
    };
    let out = run_ffid("ffid_im", &ffid_base(), options);
    compare("ffid_im", &out, FFID_IGNORE);
}

#[test]
fn ffid_im_off() {
    let options = FfidOptions {
        im: true,
        ..FfidOptions::default()
    };
    let params = with(ffid_base(), &[("extract:IM_window", float(0.0))]);
    let out = run_ffid("ffid_im_off", &params, options);
    compare("ffid_im_off", &out, FFID_IGNORE);
}

#[test]
fn mi_tsv() {
    let out = run_metabo("mi_tsv", &mi_base(), 0, false, false);
    compare("mi_tsv", &out, IGNORE);
}

#[test]
fn mi_none_candidates() {
    let params = with(
        mi_base(),
        &[
            ("model:type", string("none")),
            ("candidates_out", string("__out__")),
        ],
    );
    let out = run_metabo("mi_none_candidates", &params, 0, false, false);
    compare("mi_none_candidates", &out, IGNORE);
}

#[test]
fn mi_asymmetric() {
    let params = with(mi_base(), &[("model:type", string("asymmetric"))]);
    let out = run_metabo("mi_asymmetric", &params, 0, false, false);
    compare("mi_asymmetric", &out, IGNORE);
}

#[test]
fn mi_pmin() {
    let params = with(mi_base(), &[("extract:isotope_pmin", float(0.1))]);
    let out = run_metabo("mi_pmin", &params, 1, false, false);
    compare("mi_pmin", &out, IGNORE);
}

#[test]
fn mi_variants() {
    let out = run_metabo("mi_variants", &mi_base(), 1, false, false);
    compare("mi_variants", &out, IGNORE);
}

#[test]
fn mi_variants_none() {
    let params = with(mi_base(), &[("model:type", string("none"))]);
    let out = run_metabo("mi_variants_none", &params, 1, false, false);
    compare("mi_variants_none", &out, IGNORE);
}

#[test]
fn mi_faims() {
    let out = run_metabo("mi_faims", &mi_base(), 0, true, false);
    compare("mi_faims", &out, IGNORE);
}

#[test]
fn mi_im() {
    let params = with(mi_base(), &[("extract:im_window", float(0.05))]);
    let out = run_metabo("mi_im", &params, 1, false, true);
    compare("mi_im", &out, IGNORE);
}

#[test]
fn mi_im_off() {
    let out = run_metabo("mi_im_off", &mi_base(), 1, false, true);
    compare("mi_im_off", &out, IGNORE);
}

#[test]
fn mi_no_candidates_is_the_sources_empty_range_error() {
    // The source's `FeatureOverlapFilter::filter` reads the m/z range of the
    // map, which throws for an empty one: a run that finds no candidate at all
    // ends with `Exception::InvalidRange` ("Empty or uninitalized range
    // object") instead of an empty result (oracle case `mi_no_candidates`,
    // exit code 2).
    let table = vec![FeatureFinderMetaboIdentCompound::new(
        "no_signal",
        "C6H12O6",
        0.0,
        vec![1],
        vec![50.0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        "",
    )];
    let param = apply(
        &mi_base(),
        FeatureFinderAlgorithmMetaboIdent::defaults().unwrap(),
    );
    let mut algorithm = FeatureFinderAlgorithmMetaboIdent::from_param(&param).unwrap();
    algorithm.set_ms_data(load_ms1("FeatureFinderMetaboIdent_1_input.mzML.gz"));
    let mut features = FeatureMap::new();
    let error = algorithm
        .run(
            &table,
            &mut features,
            "FeatureFinderMetaboIdent_1_input.mzML",
        )
        .unwrap_err();
    assert!(
        matches!(error, openms::Error::InvalidRange(_)),
        "unexpected error {error:?}"
    );
    assert!(features.features.is_empty());
    let stderr = oracle_stderr("mi_no_candidates");
    assert!(
        stderr.contains("Empty or uninitalized range object"),
        "{stderr}"
    );
}

#[test]
fn mi_rt_window_default() {
    let params = vec![
        ("extract:mz_window", float(5.0)),
        ("detect:peak_width", float(3.0)),
    ];
    let out = run_metabo("mi_rt_window_default", &params, 0, false, false);
    compare("mi_rt_window_default", &out, IGNORE);
}

#[test]
fn unique_ids_follow_the_source_stream() {
    // the generator is seeded as under `-test`, and the draws happen in the
    // order of the source: offset peptides, then the features, then the map
    let out = run_ffid("ffid_none", &ffid_base(), FfidOptions::default());
    let expected = oracle_records("ffid_none");
    let mut checked = 0;
    for (key, value) in &expected {
        if key.ends_with(".uid") {
            match out.records.get(key) {
                Some(Val::S(s)) => assert_eq!(s, value, "{key}"),
                _ => panic!("missing {key}"),
            }
            checked += 1;
        }
    }
    assert!(checked > 0);
}

// ---------------------------------------------------------------- extraction

fn xic_access() -> openms::analysis::openswath::spectrum_access::SpectrumAccessPtr {
    let mut experiment = load_ms1("FeatureFinderIdentification_1_input.mzML.gz");
    add_per_peak_im(&mut experiment);
    SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr(std::sync::Arc::new(experiment))
        .unwrap()
}

/// Monoisotopic m/z and charge of the 23 features of the `ffid_none` case.
const XIC_FEATURES: [(f64, f64); 23] = [
    (461.74765303542097, 2.0),
    (395.72652290782094, 2.0),
    (455.75528053542098, 2.0),
    (569.75261839302107, 2.0),
    (443.71126790782091, 2.0),
    (487.73253447162097, 2.0),
    (325.49078180333771, 3.0),
    (554.26060201207099, 2.0),
    (646.30468413227095, 2.0),
    (431.20554824377092, 3.0),
    (379.71510077282096, 2.0),
    (532.24874658327099, 2.0),
    (404.20341232807101, 2.0),
    (300.19552859760432, 3.0),
    (368.86210990473768, 3.0),
    (449.74438990042097, 2.0),
    (300.1653520892043, 3.0),
    (395.2394634875709, 2.0),
    (501.79513472682095, 2.0),
    (358.17457648633768, 3.0),
    (421.75835453542095, 2.0),
    (613.30659575872096, 2.0),
    (464.25036251947103, 2.0),
];

/// Four isotope traces of every feature, one sorted walk, with and without the
/// ion-mobility window; `mixed` takes the ion mobility from every third
/// feature and gives every other one an RT window, as the driver does.
fn xic_all(mixed: bool) -> Out {
    let access = xic_access();
    let mut coords: Vec<ExtractionCoordinates> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    for (i, &(mz, z)) in XIC_FEATURES.iter().enumerate() {
        for k in 0..4 {
            let mut c = ExtractionCoordinates {
                mz: mz + f64::from(k) * 1.0033548378 / z,
                ion_mobility: im_of_mz(mz),
                mz_precursor: mz,
                rt_start: -1.0,
                rt_end: -1.0,
                id: format!("f{i}_iso{k}"),
            };
            if mixed {
                if i % 3 == 0 {
                    c.ion_mobility = -1.0;
                }
                if i % 2 == 1 {
                    c.rt_start = 1800.0 + 20.0 * i as f64;
                    c.rt_end = c.rt_start + 400.0;
                }
            }
            coords.push(c);
        }
    }
    coords.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
    let mut out = Out::default();
    out.pi("ncoords", coords.len() as i64);
    for mode in 0..2 {
        let mut chroms: Vec<Chromatogram> =
            (0..coords.len()).map(|_| Chromatogram::new()).collect();
        ChromatogramExtractorAlgorithm::new()
            .extract_chromatograms(
                access.as_ref(),
                &mut chroms,
                &coords,
                0.1,
                false,
                if mode == 0 { 0.0 } else { 0.06 },
                "tophat",
            )
            .unwrap();
        for (coord, chrom) in coords.iter().zip(&chroms) {
            let pre = format!("mode{mode}.{}", coord.id);
            out.p(&format!("{pre}.mz"), coord.mz);
            out.pi(&format!("{pre}.n"), chrom.time_array().data.len() as i64);
            out.p(&format!("{pre}.sum"), chrom.intensity_array().data.iter().sum());
        }
    }
    out
}

/// The extractor with and without the ion-mobility window on the four isotope
/// traces of AEFVEVTK/2, with the points of every trace.
#[test]
fn xic_im() {
    let access = xic_access();
    let base_mz = 461.7476530354209;
    let im = im_of_mz(base_mz);
    let mut out = Out::default();
    out.p("im", im);
    out.pi("nspec", access.nr_spectra() as i64);
    for mode in 0..2 {
        let coords: Vec<ExtractionCoordinates> = (0..4)
            .map(|k| ExtractionCoordinates {
                mz: base_mz + f64::from(k) * 1.0033548378 / 2.0,
                ion_mobility: im,
                mz_precursor: base_mz,
                rt_start: 1900.0,
                rt_end: 2150.0,
                id: format!("iso{k}"),
            })
            .collect();
        let mut chroms: Vec<Chromatogram> = (0..4).map(|_| Chromatogram::new()).collect();
        ChromatogramExtractorAlgorithm::new()
            .extract_chromatograms(
                access.as_ref(),
                &mut chroms,
                &coords,
                0.1,
                false,
                if mode == 0 { 0.0 } else { 0.06 },
                "tophat",
            )
            .unwrap();
        for (k, chrom) in chroms.iter().enumerate() {
            let pre = format!("mode{mode}.c{k}");
            let times = &chrom.time_array().data;
            let values = &chrom.intensity_array().data;
            out.pi(&format!("{pre}.n"), times.len() as i64);
            out.p(&format!("{pre}.sum"), values.iter().sum());
            for (i, (t, v)) in times.iter().zip(values).enumerate() {
                out.p(&format!("{pre}.t{i}"), *t);
                out.p(&format!("{pre}.v{i}"), *v);
            }
        }
    }
    compare("xic_im", &out, &[]);
}

/// Every coordinate carries ion mobility: the window removes only the isotope
/// traces of one feature, which lie across the 1/K0 wrap of the test data.
#[test]
fn xic_im_all() {
    let out = xic_all(false);
    compare("xic_im_all", &out, &[]);
}

/// Coordinates with and without ion mobility in one walk: the source advances
/// its ion-mobility iterator only in the ion-mobility overload, so every
/// ion-mobility coordinate after the first m/z-only one reads the ion mobility
/// of an earlier peak. The Release build does it, and so does the port.
#[test]
fn xic_im_mixed() {
    let out = xic_all(true);
    compare("xic_im_mixed", &out, &[]);
}

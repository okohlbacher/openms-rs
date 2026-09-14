// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `FileHandler::get_type`, `load_experiment_with_options` and
//! `load_feature_map_with_options`.
//!
//! Type detection, the FeatureFinderCentroided loading options, the intensity
//! range edges and the feature-file options are compared with
//! `data/mzml_mobility/a3_format_io_oracle.tsv` (oracle-generated, tier 1
//! executed differential against product-sdk; `../oracle/a3-format-io/`). The
//! DTA2D, DTA and MGF option routing is source review of `FileHandler.cpp:869-936`
//! and `DTA2DFile.h:208-243`. See `docs/MZML_MOBILITY_SUPPORT.md`.

use openms::format::{FileHandler, FileType, PeakFileOptions};
use openms::kernel::NumericRange;
use openms::{Error, MSExperiment, MSSpectrum, Peak1D};

const ORACLE: &str = include_str!("data/mzml_mobility/a3_format_io_oracle.tsv");

struct Directory(std::path::PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "openms-rs-type-detection-{}-{id}",
            std::process::id()
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn fixture(name: &str) -> String {
    format!(
        "{}/tests/data/mzml_mobility/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// Oracle records of one kind, split into fields; `label` filters the second field.
fn rows(record: &str, label: Option<&str>) -> Vec<Vec<&'static str>> {
    ORACLE
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[0] == record && label.is_none_or(|l| fields.get(1) == Some(&l)))
        .collect()
}

/// Parse a C `printf("%a")` hexadecimal float exactly.
#[cfg(feature = "mzml")]
fn hex(text: &str) -> f64 {
    let (negative, body) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let body = body.strip_prefix("0x").expect("hexadecimal float");
    let (mantissa, exponent) = body.split_once('p').expect("binary exponent");
    let exponent: i32 = exponent.parse().expect("decimal exponent");
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = u64::from_str_radix(&format!("{whole}{fraction}"), 16).expect("hex digits");
    let scale = exponent - 4 * i32::try_from(fraction.len()).expect("short fraction");
    let magnitude = digits as f64 * 2f64.powi(scale);
    if negative { -magnitude } else { magnitude }
}

#[test]
fn get_type_matches_the_oracle_for_names_directories_and_content() {
    let dir = Directory::new();
    let root = dir.0.to_str().unwrap().to_owned();
    for name in ["FileInfo_5_input.mzDat", "FileInfo_8_input.notype"] {
        std::fs::copy(fixture(name), dir.0.join(name)).unwrap();
    }
    std::fs::create_dir(dir.0.join("sample.d")).unwrap();
    std::fs::write(dir.0.join("sample.d/analysis.tdf"), b"").unwrap();
    std::fs::write(dir.0.join("empty.unknownext"), b"").unwrap();
    std::fs::copy(
        fixture("spectrum_representation.mzML"),
        dir.0.join("mzml_content.unknownext"),
    )
    .unwrap();
    // The oracle gzips spectrum_representation.mzML; any gzip mzML sniffs alike.
    #[cfg(feature = "mzml")]
    {
        let written = dir.0.join("gz_mzml.mzML.gz");
        let experiment = MSExperiment {
            spectra: vec![MSSpectrum {
                native_id: "scan=1".into(),
                peaks: vec![Peak1D::new(100.0, 1.0)],
                ..Default::default()
            }],
            ..Default::default()
        };
        FileHandler::store_experiment(&written, &experiment, None).unwrap();
        std::fs::rename(&written, dir.0.join("gz_mzml.bin")).unwrap();
    }
    let records = rows("type", None);
    assert_eq!(records.len(), 13);
    for record in records {
        let (name, expected) = (record[1], record[2]);
        if name == "gz_mzml.bin" && !cfg!(feature = "mzml") {
            continue;
        }
        let actual = FileHandler::get_type(format!("{root}/{name}"));
        if let Some(error) = expected.strip_prefix("error:") {
            assert_eq!(error, "FileNotFound");
            assert!(
                matches!(&actual, Err(Error::Io(e)) if e.kind() == std::io::ErrorKind::NotFound),
                "{name}: {actual:?}"
            );
        } else {
            assert_eq!(actual.unwrap(), FileType::from_name(expected), "{name}");
        }
    }
    // Unknown content stays unknown, and an unknown extension is resolved by content.
    assert_eq!(
        FileHandler::get_type(format!("{root}/FileInfo_8_input.notype")).unwrap(),
        FileType::Unknown
    );
    assert_eq!(
        FileHandler::get_type(format!("{root}/FileInfo_5_input.mzDat")).unwrap(),
        FileType::MzData
    );
}

#[test]
fn get_type_reads_no_file_for_known_or_directory_names() {
    // A known extension, including behind separators or compression suffixes,
    // and every Bruker TDF name are decided lexically.
    for (name, expected) in [
        ("/nonexistent/peaks.mzML", FileType::MzMl),
        ("/nonexistent/folder.featureXML//", FileType::FeatureXml),
        ("/nonexistent/back.mzML\\\\", FileType::MzMl),
        ("/nonexistent/features.featureXML.bz2", FileType::FeatureXml),
        ("/nonexistent/sample.d", FileType::Unknown),
        ("/nonexistent/sample.d/", FileType::Unknown),
        ("/nonexistent/sample.D.zip", FileType::Unknown),
    ] {
        assert_eq!(FileHandler::get_type(name).unwrap(), expected, "{name}");
    }
    assert!(FileHandler::get_type("/nonexistent/file.unknownext").is_err());
}

#[cfg(feature = "mzml")]
fn feature_finder_options(min: f64) -> PeakFileOptions {
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_intensity_range(NumericRange { min, max: f64::MAX });
    options
}

#[cfg(feature = "mzml")]
#[test]
fn feature_finder_options_load_the_ffc_1_input() {
    let path = fixture("FeatureFinderCentroided_1_input.mzML");
    let allowed = [FileType::MzMl, FileType::Raw];
    // MS1 spectra with intensities in [f64::MIN_POSITIVE, f64::MAX).
    let intended = FileHandler::load_experiment_with_options(
        &path,
        &allowed,
        &feature_finder_options(f64::MIN_POSITIVE),
    )
    .unwrap();
    assert_eq!(intended.spectra.len(), 112);
    assert_eq!(
        intended.spectra.iter().map(MSSpectrum::len).sum::<usize>(),
        3084
    );
    assert!(
        intended
            .settings
            .document
            .loaded_file_path
            .ends_with("FeatureFinderCentroided_1_input.mzML")
    );
    // FeatureFinderCentroided.cpp:184 passes std::numeric_limits<DPosition<1>>::min().
    // DPosition has no numeric_limits specialisation, so that is DPosition() = 0
    // and the executed range is [0, DBL_MAX); this input has no zero peak, so
    // both ranges agree here.
    let executed =
        FileHandler::load_experiment_with_options(&path, &allowed, &feature_finder_options(0.0))
            .unwrap();
    let totals = rows("filtered", Some("FeatureFinderCentroided_1_input"));
    assert_eq!(executed.spectra.len().to_string(), totals[0][2]);
    assert_eq!(
        executed
            .spectra
            .iter()
            .map(MSSpectrum::len)
            .sum::<usize>()
            .to_string(),
        totals[0][3]
    );
    let per_spectrum = rows("filtered_spectrum", Some("FeatureFinderCentroided_1_input"));
    assert_eq!(per_spectrum.len(), executed.spectra.len());
    for (spectrum, row) in executed.spectra.iter().zip(per_spectrum) {
        assert_eq!(spectrum.native_id, row[3]);
        assert_eq!(spectrum.len().to_string(), row[4]);
    }
    // A detected type outside the allowed list: the source throws ParseError.
    assert_eq!(
        rows("experiment_disallowed", Some("spectrum_representation"))[0][2],
        "Parse Error"
    );
    assert!(matches!(
        FileHandler::load_experiment_with_options(
            fixture("spectrum_representation.mzML"),
            &[FileType::Mgf],
            &PeakFileOptions::default()
        ),
        Err(Error::InvalidValue(_))
    ));
}

#[cfg(feature = "mzml")]
#[test]
fn intensity_range_edges_match_the_oracle() {
    type Point = (usize, u64, u32);
    fn points(experiment: &MSExperiment) -> Vec<Point> {
        experiment
            .spectra
            .iter()
            .enumerate()
            .flat_map(|(index, s)| {
                s.peaks
                    .iter()
                    .map(move |p| (index, p.mz.to_bits(), p.intensity.to_bits()))
            })
            .collect()
    }
    fn oracle(record: &str) -> Vec<Point> {
        rows(record, Some("intensity_bounds"))
            .iter()
            .map(|r| {
                (
                    r[2].parse().unwrap(),
                    hex(r[4]).to_bits(),
                    (hex(r[5]) as f32).to_bits(),
                )
            })
            .collect()
    }
    let path = fixture("intensity_bounds.mzML");
    let plain =
        FileHandler::load_experiment_with_options(&path, &[], &PeakFileOptions::default()).unwrap();
    assert_eq!(points(&plain), oracle("spectrum_peak"));
    // Source-executed range [0, DBL_MAX): 0.0 and a subnormal stay, -1 goes.
    let executed =
        FileHandler::load_experiment_with_options(&path, &[], &feature_finder_options(0.0))
            .unwrap();
    assert_eq!(executed.spectra.len(), 1);
    assert_eq!(points(&executed), oracle("filtered_detail_peak"));
    // The range the source comment intends also drops 0.0 and 1e-310; the
    // minimum itself, f64::MIN_POSITIVE, is kept and stored as f32 zero.
    let intended = FileHandler::load_experiment_with_options(
        &path,
        &[],
        &feature_finder_options(f64::MIN_POSITIVE),
    )
    .unwrap();
    assert_eq!(
        intended.spectra[0]
            .peaks
            .iter()
            .map(|p| (p.mz, p.intensity))
            .collect::<Vec<_>>(),
        [(103.0, 0.0), (104.0, 1.0), (105.0, 1e30)]
    );
}

#[test]
fn dta2d_receives_only_the_three_ranges_the_source_reads() {
    let dir = Directory::new();
    let path = dir.0.join("peaks.dta2d");
    std::fs::write(
        &path,
        "#SEC\tMZ\tINT\n10\t100\t5\n10\t200\t50\n20\t100\t7\n",
    )
    .unwrap();
    let load = |options: &PeakFileOptions| {
        FileHandler::load_experiment_with_options(&path, &[FileType::Dta2d], options).unwrap()
    };
    let summary = |experiment: &MSExperiment| {
        experiment
            .spectra
            .iter()
            .map(|s| (s.rt, s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>()))
            .collect::<Vec<_>>()
    };
    let all = summary(&load(&PeakFileOptions::default()));
    assert_eq!(all, [(10.0, vec![100.0, 200.0]), (20.0, vec![100.0])]);
    let mut rt = PeakFileOptions::default();
    rt.set_rt_range(NumericRange {
        min: 15.0,
        max: 25.0,
    });
    assert_eq!(summary(&load(&rt)), [(20.0, vec![100.0])]);
    // DTA2DFile.h:235-243 drops a last spectrum the m/z range leaves empty.
    let mut mz = PeakFileOptions::default();
    mz.set_mz_range(NumericRange {
        min: 150.0,
        max: 250.0,
    });
    assert_eq!(summary(&load(&mz)), [(10.0, vec![200.0])]);
    // Half-open: an intensity equal to the maximum is dropped.
    let mut intensity = PeakFileOptions::default();
    intensity.set_intensity_range(NumericRange {
        min: 6.0,
        max: 50.0,
    });
    assert_eq!(summary(&load(&intensity)), [(20.0, vec![100.0])]);
    // DTA2DFile ignores MS levels.
    let mut levels = PeakFileOptions::default();
    levels.add_ms_level(2).unwrap();
    assert_eq!(summary(&load(&levels)), all);
}

#[test]
fn dta_and_mgf_ignore_peak_file_options_like_the_source() {
    let dir = Directory::new();
    let mut options = PeakFileOptions::default();
    options.add_ms_level(1).unwrap();
    options.set_mz_range(NumericRange { min: 0.0, max: 1.0 });
    options.set_intensity_range(NumericRange {
        min: 1000.0,
        max: 2000.0,
    });
    let mgf = dir.0.join("peaks.mgf");
    let experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            ms_level: 2,
            peaks: vec![Peak1D::new(100.5, 7.0)],
            ..Default::default()
        }],
        ..Default::default()
    };
    FileHandler::store_experiment(&mgf, &experiment, None).unwrap();
    let dta = dir.0.join("peaks.dta");
    std::fs::write(&dta, "501.0 1\n100.5 7\n").unwrap();
    for (path, kind) in [(mgf, FileType::Mgf), (dta, FileType::Dta)] {
        let loaded = FileHandler::load_experiment_with_options(&path, &[kind], &options).unwrap();
        assert_eq!(loaded.spectra.len(), 1, "{kind:?}");
        assert_eq!(
            loaded.spectra[0].peaks, experiment.spectra[0].peaks,
            "{kind:?}"
        );
    }
}

#[cfg(feature = "featurexml")]
#[test]
fn feature_file_options_match_the_oracle() {
    use openms::format::featurexml::FeatureFileOptions;
    use openms::kernel::FeatureMap;
    fn counts(map: &FeatureMap) -> [String; 3] {
        [
            map.features.len().to_string(),
            map.features
                .iter()
                .map(|f| f.convex_hulls.len())
                .sum::<usize>()
                .to_string(),
            map.features
                .iter()
                .map(|f| f.subordinates.len())
                .sum::<usize>()
                .to_string(),
        ]
    }
    let inputs = [
        (
            "FeatureFinderCentroided_1_1_output",
            fixture("FeatureFinderCentroided_1_1_output.featureXML"),
        ),
        (
            "featurexml_source_1",
            format!(
                "{}/tests/data/featurexml_source_1.featureXML",
                env!("CARGO_MANIFEST_DIR")
            ),
        ),
    ];
    for (label, path) in inputs {
        let oracle = rows("features", Some(label));
        for (mode, options) in [
            ("default", FeatureFileOptions::default()),
            (
                "no_hulls_no_subordinates",
                FeatureFileOptions {
                    load_convex_hulls: false,
                    load_subordinates: false,
                    ..Default::default()
                },
            ),
        ] {
            let row = oracle.iter().find(|r| r[2] == mode).expect("oracle mode");
            let map = FileHandler::load_feature_map_with_options(&path, &options).unwrap();
            assert_eq!(counts(&map), [row[3], row[4], row[5]], "{label} {mode}");
            assert_eq!(map.loaded_file_type, FileType::FeatureXml);
        }
        // The source's allowed-type check; the native call has no list to fail.
        assert_eq!(
            rows("features_disallowed", Some(label))[0][2],
            "InvalidFileType"
        );
    }
    assert!(matches!(
        FileHandler::load_feature_map_with_options(
            fixture("FeatureFinderCentroided_1_input.mzML"),
            &FeatureFileOptions::default()
        ),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn options_loader_rejects_types_without_an_experiment_adapter() {
    assert!(matches!(
        FileHandler::load_experiment_with_options(
            fixture("FeatureFinderCentroided_1_1_output.featureXML"),
            &[],
            &PeakFileOptions::default()
        ),
        Err(Error::Unsupported(_))
    ));
}

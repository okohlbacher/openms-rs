// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::{dta, dta2d, ms2};
use openms::{MSExperiment, MSSpectrum, Peak1D, Precursor};
use std::io::{self, Cursor, Write};
const MS2: &[u8] = include_bytes!("data/text_peak_lists/MS2File_test_spectra.ms2");
const DTA: [&[u8]; 3] = [
    include_bytes!("data/text_peak_lists/DTA2DFile_test_1.dta2d"),
    include_bytes!("data/text_peak_lists/DTA2DFile_test_2.dta2d"),
    include_bytes!("data/text_peak_lists/DTA2DFile_test_3.dta2d"),
];
fn peak_spectrum(rt: f64) -> MSSpectrum {
    MSSpectrum {
        rt,
        peaks: vec![Peak1D::new(123.456789012345, 42.25)],
        ..Default::default()
    }
}
fn experiment(s: MSSpectrum) -> MSExperiment {
    MSExperiment {
        spectra: vec![s],
        ..Default::default()
    }
}

#[test]
fn ms2_literal_source_fixture_ignores_charge_and_analysis_lines() {
    let e = ms2::read(MS2).unwrap();
    assert_eq!(e.len(), 2);
    for (i, mz) in [444.44, 555.555].into_iter().enumerate() {
        let s = &e.spectra[i];
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.native_id, format!("index={i}"));
        assert_eq!(s.rt, -1.);
        assert_eq!(s.precursors, [Precursor::new(mz, 0)]);
        assert_eq!(
            s.peaks,
            vec![
                Peak1D::new(1., 2.3),
                Peak1D::new(2., 3.4),
                Peak1D::new(3., 4.5),
                Peak1D::new(6., 9.)
            ]
        );
    }
    assert!(e.settings.metadata.is_empty());
    let mut text = Vec::new();
    ms2::write(&mut text, &e).unwrap();
    assert_eq!(ms2::read(text.as_slice()).unwrap(), e);
}

#[test]
fn ms2_source_scan_tokens_empty_scans_and_discarded_pre_scan_peaks() {
    let e=ms2::read(&b"1 2\nH ignored\nS first last -3\nS arbitrary ignored +4\nI arbitrary\nZ nonsense\nD anything\n 5\t -6 \r\n"[..]).unwrap();
    assert_eq!(e.len(), 2);
    assert!(e.spectra[0].peaks.is_empty());
    assert_eq!(e.spectra[0].precursors[0].mz, -3.);
    assert_eq!(e.spectra[1].peaks, [Peak1D::new(5., -6.)]);
    assert!(ms2::read(&b"1 2\nH ignored\n"[..]).unwrap().is_empty());
    for data in [
        "S 1 2\n",
        "S 1 2 3 4\n",
        "S 1 1 4\n1 2 3\n",
        "bad peak\nS 1 1 1\n",
    ] {
        assert!(ms2::read(data.as_bytes()).is_err(), "{data}");
    }
}

#[test]
fn dta2d_all_three_source_fixtures_preserve_values_and_minutes() {
    let first = dta2d::read(DTA[0]).unwrap();
    let second = dta2d::read(DTA[1]).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.len(), 9);
    assert_eq!(
        first.spectra.iter().map(|s| s.peaks.len()).sum::<usize>(),
        11
    );
    let mz = [
        230.02, 231.51, 139.42, 149.93, 169.65, 189.30, 202.28, 207.82, 219.72,
    ];
    let intensity = [
        47218.89_f32,
        89935.22,
        318.52,
        61870.99,
        62074.22,
        53737.85,
        49410.25,
        17038.71,
        73629.98,
    ];
    for i in 0..9 {
        let s = &first.spectra[i];
        assert_eq!(s.native_id, format!("index={i}"));
        assert_eq!(
            s.rt,
            [
                4711.1, 4711.2, 4711.3, 4711.4, 4711.5, 4711.6, 4711.7, 4711.8, 4711.9
            ][i]
        );
        assert_eq!(s.peaks[0], Peak1D::new(mz[i], intensity[i]));
    }
    let minutes = dta2d::read(DTA[2]).unwrap();
    for (s, expected) in minutes
        .spectra
        .iter()
        .zip([282666., 282672., 282678., 282684., 282690.])
    {
        assert_eq!(s.rt, expected);
    }
    let mut text = Vec::new();
    dta2d::write(&mut text, &first).unwrap();
    assert_eq!(dta2d::read(text.as_slice()).unwrap(), first);
}

#[test]
fn dta2d_source_filters_are_half_open_and_ids_keep_gaps() {
    let cases = [
        (
            dta2d::ReadOptions {
                rt_range: Some(4711.15..4711.45),
                ..Default::default()
            },
            vec![1, 2, 3],
        ),
        (
            dta2d::ReadOptions {
                mz_range: Some(150.0..220.0),
                ..Default::default()
            },
            vec![4, 5, 6, 7, 8],
        ),
        (
            dta2d::ReadOptions {
                intensity_range: Some(30000.0..70000.0),
                ..Default::default()
            },
            vec![0, 3, 4, 5, 6],
        ),
    ];
    for (options, ids) in cases {
        let e = dta2d::read_with_options(DTA[0], &options).unwrap();
        assert_eq!(
            e.spectra
                .iter()
                .map(|s| s.native_id.clone())
                .collect::<Vec<_>>(),
            ids.iter().map(|i| format!("index={i}")).collect::<Vec<_>>()
        );
    }
    let options = dta2d::ReadOptions {
        rt_range: Some(1.0..2.0),
        mz_range: Some(10.0..20.0),
        intensity_range: Some(100.0..200.0),
        ..Default::default()
    };
    let e = dta2d::read_with_options(&b"1 10 100\n1 20 100\n1 10 200\n2 10 100\n"[..], &options)
        .unwrap();
    assert_eq!(e.len(), 1);
    assert_eq!(e.spectra[0].peaks, [Peak1D::new(10., 100.)]);
    assert!(
        dta2d::read_with_options(
            &b""[..],
            &dta2d::ReadOptions {
                rt_range: Some(2.0..1.0),
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn dta2d_header_delimiter_and_sticky_minute_conventions() {
    let e = dta2d::read(&b"# Int Mz Min ignored\n2 3 1\n#INT MZ SEC\n4 5 2\n"[..]).unwrap();
    assert_eq!(e.spectra[0].rt, 60.);
    assert_eq!(e.spectra[1].rt, 120.);
    assert_eq!(e.spectra[1].peaks, [Peak1D::new(5., 4.)]);
    let e = dta2d::read(&b"#INTENSITY\tMASS-TO-CHARGE\tRETENTION_TIME\n-1\t-2\t-3\n"[..]).unwrap();
    assert_eq!(e.spectra[0].rt, -3.);
    assert_eq!(e.spectra[0].peaks, [Peak1D::new(-2., -1.)]);
    for text in [
        "#SEC MZ\n",
        "#SEC SEC INT\n",
        "#SEC  MZ INT\n",
        "1  2 3\n",
        "1\t2 3\n",
        "# comment\n",
    ] {
        assert!(dta2d::read(text.as_bytes()).is_err(), "{text}");
    }
}

#[test]
fn dta2d_anchor_tolerance_initial_sentinel_and_group_limits() {
    let e = dta2d::read(&b"0 1 1\n0.0001 2 2\n0.00015 3 3\n"[..]).unwrap();
    assert_eq!(e.len(), 2);
    assert_eq!(e.spectra[0].rt, 0.);
    assert_eq!(e.spectra[0].peaks.len(), 2);
    let e = dta2d::read(&b"-0.99995 1 1\n0 2 2\n"[..]).unwrap();
    assert_eq!(e.spectra[0].rt, -1.);
    assert!(e.spectra[0].native_id.is_empty());
    assert_eq!(e.spectra[1].native_id, "index=0");
    let options = dta2d::ReadOptions {
        limits: ms2::Limits {
            max_spectra: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(dta2d::read_with_options(&b"-1 1 1\n0 2 2\n"[..], &options).is_err());
}

#[test]
fn direct_f32_rounding_and_nonfinite_or_underflow_errors() {
    let value = "1.0000000596046447753906250000000000000001";
    let text = format!("S 1 1 1\n2 {value}\n");
    assert_eq!(
        ms2::read(text.as_bytes()).unwrap().spectra[0].peaks[0]
            .intensity
            .to_bits(),
        0x3f800001
    );
    let text = format!("1 2 {value}\n");
    assert_eq!(
        dta2d::read(text.as_bytes()).unwrap().spectra[0].peaks[0]
            .intensity
            .to_bits(),
        0x3f800001
    );
    for value in ["nan", "inf", "1e1000", "1e-1000"] {
        assert!(ms2::read(format!("S 1 1 {value}\n").as_bytes()).is_err());
        assert!(dta2d::read(format!("1 2 {value}\n").as_bytes()).is_err());
    }
    assert!(dta2d::read(&b"#MIN MZ INT\n1e308 2 3\n"[..]).is_err());
}

#[test]
fn limits_include_ignored_lines_and_filtered_data_and_replacement_is_atomic() {
    let e = experiment(peak_spectrum(42.));
    for options in [
        ms2::Limits {
            max_line_bytes: 3,
            ..Default::default()
        },
        ms2::Limits {
            max_bytes: 5,
            ..Default::default()
        },
        ms2::Limits {
            max_lines: 1,
            ..Default::default()
        },
        ms2::Limits {
            max_peaks: 0,
            ..Default::default()
        },
        ms2::Limits {
            max_spectra: 0,
            ..Default::default()
        },
    ] {
        let mut target = e.clone();
        assert!(
            ms2::read_into_with_options(&b"S 1 1 4\n1 2\n"[..], &mut target, &options).is_err()
        );
        assert_eq!(target, e);
    }
    let mut target = e.clone();
    assert!(dta2d::read_into(&b"1 2 3\n4 bad 6\n"[..], &mut target).is_err());
    assert_eq!(target, e);
    let options = dta2d::ReadOptions {
        mz_range: Some(0.0..1.0),
        limits: ms2::Limits {
            max_peaks: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(dta2d::read_with_options(&b"1 2 3\n1 2 4\n"[..], &options).is_err());
    assert!(
        ms2::read_with_options(
            &b"H ignored\nH ignored\n"[..],
            &ms2::Limits {
                max_lines: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn writers_reject_unrepresentable_data_before_output() {
    let mut cases = Vec::new();
    let mut s = peak_spectrum(1.);
    s.name = "name".into();
    cases.push(experiment(s));
    let mut s = peak_spectrum(1.);
    s.native_id = "scan=9".into();
    cases.push(experiment(s));
    cases.push(experiment(MSSpectrum::default()));
    cases.push(MSExperiment {
        spectra: vec![peak_spectrum(1.), peak_spectrum(1.00001)],
        ..Default::default()
    });
    cases.push(experiment(peak_spectrum(-0.99995)));
    let mut s = peak_spectrum(1.);
    s.precursors.push(Precursor::new(500., 2));
    cases.push(experiment(s));
    for e in cases {
        let mut bytes = vec![9];
        assert!(dta2d::write(&mut bytes, &e).is_err());
        assert_eq!(bytes, [9]);
    }
    let mut e = ms2::read(MS2).unwrap();
    e.spectra[1].precursors[0].charge = 2;
    let mut bytes = vec![9];
    assert!(ms2::write(&mut bytes, &e).is_err());
    assert_eq!(bytes, [9]);
    let e = experiment(peak_spectrum(1.));
    assert!(
        dta2d::write_with_options(
            &mut bytes,
            &e,
            &dta2d::WriteOptions::from(ms2::Limits {
                max_output_bytes: 2,
                ..Default::default()
            })
        )
        .is_err()
    );
    assert_eq!(bytes, [9]);
}

#[test]
fn tic_matches_source_ms1_only_f32_accumulation_and_explicit_projection() {
    let e = dta2d::read(DTA[0]).unwrap();
    let mut text = Vec::new();
    dta2d::write_tic(&mut text, &e).unwrap();
    let tic = dta2d::read(text.as_slice()).unwrap();
    assert_eq!(tic.len(), 9);
    assert_eq!(
        tic.spectra[0].peaks,
        [Peak1D::new(0., f32::from_bits(0x480a546b))]
    );
    let mut s = peak_spectrum(1.);
    s.peaks = vec![
        Peak1D::new(f64::NAN, 16777216.),
        Peak1D::new(0., 1.),
        Peak1D::new(0., 1.),
    ];
    let mut higher = peak_spectrum(2.);
    higher.ms_level = 2;
    higher.peaks[0].intensity = f32::NAN;
    let e = MSExperiment {
        spectra: vec![
            s,
            higher,
            MSSpectrum {
                rt: 3.,
                ..Default::default()
            },
        ],
        settings: openms::metadata::ExperimentalSettings {
            metadata: [("discarded".into(), "for explicit TIC".into())].into(),
            ..Default::default()
        },
        ..Default::default()
    };
    text.clear();
    dta2d::write_tic(&mut text, &e).unwrap();
    let result = dta2d::read(text.as_slice()).unwrap();
    assert_eq!(result.len(), 2);
    assert_eq!(result.spectra[0].peaks[0].intensity, 16777216.);
    assert_eq!(result.spectra[1].peaks[0].intensity, 0.);
    let mut s = peak_spectrum(1.);
    s.peaks = vec![Peak1D::new(0., f32::MAX), Peak1D::new(0., f32::MAX)];
    text.clear();
    assert!(dta2d::write_tic(&mut text, &experiment(s)).is_err());
    assert!(text.is_empty());
}

struct FlushFailure;
impl Write for FlushFailure {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("flush"))
    }
}
#[test]
fn flush_errors_and_existing_file_preservation() {
    let e = dta2d::read(DTA[0]).unwrap();
    assert!(dta2d::write(FlushFailure, &e).is_err());
    assert!(dta2d::write_tic(FlushFailure, &e).is_err());
    assert!(ms2::write(FlushFailure, &ms2::read(MS2).unwrap()).is_err());
    let path = std::env::temp_dir().join(format!(
        "openms-flat-preflight-{}.dta2d",
        std::process::id()
    ));
    std::fs::write(&path, b"existing").unwrap();
    assert!(dta2d::store(&path, &experiment(MSSpectrum::default())).is_err());
    assert_eq!(std::fs::read(&path).unwrap(), b"existing");
    std::fs::remove_file(path).unwrap();
    let mut cursor = Cursor::new(Vec::new());
    ms2::write(&mut cursor, &ms2::read(MS2).unwrap()).unwrap();
    assert!(!cursor.into_inner().is_empty());
}

fn data_path(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/text_peak_lists")
        .join(name)
}

/// Tier 1: the Release build's `DTA2DFile::store` and `storeTIC` output for
/// `DTA2DFile_test_1.dta2d`, retained by `../oracle/progress-format-readers`
/// (ibminode06, three identical repetitions). Both the default writer and the
/// source mode write those bytes: every number of this input reads back.
#[test]
fn dta2d_store_writes_the_release_build_bytes() {
    let e = dta2d::load(data_path("DTA2DFile_test_1.dta2d")).unwrap();
    let store = std::fs::read(data_path("dta2d_store_release.dta2d")).unwrap();
    let tic = std::fs::read(data_path("dta2d_store_tic_release.dta2d")).unwrap();
    for options in [
        dta2d::WriteOptions::default(),
        dta2d::WriteOptions::source(),
    ] {
        let mut bytes = Vec::new();
        dta2d::write_with_options(&mut bytes, &e, &options).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            String::from_utf8(store.clone()).unwrap()
        );
        let mut bytes = Vec::new();
        dta2d::write_tic_with_options(&mut bytes, &e, &options).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            String::from_utf8(tic.clone()).unwrap()
        );
    }
    // The source text reads back as the loaded values.
    assert_eq!(dta2d::read(store.as_slice()).unwrap(), e);
    // Three literal lines of the Release file, for a reader of this test:
    // 15 fraction digits for a double, 6 fraction digits or the shortest
    // scientific text for a float.
    let text = String::from_utf8(store).unwrap();
    assert!(
        text.starts_with("#SEC\tMZ\tINT\n4711.100000000000364\t230.02000000000001\t4.721889e04\n")
    );
    assert!(text.contains("\n4711.300000000000182\t139.419999999999987\t318.519989\n"));
    assert!(
        String::from_utf8(tic)
            .unwrap()
            .contains("\n4711.100000000000364\t0\t1.4164967e05\n")
    );
}

/// The source text keeps 15 fraction digits of a double and 6 of a float in
/// `[0.01, 10000)`. The default writer refuses a value that text would change;
/// the source mode writes it as the source does.
#[test]
fn dta2d_numbers_the_source_text_would_change_are_refused_by_default() {
    let mut small_mz = peak_spectrum(1.);
    small_mz.peaks[0].mz = 1.234_567_890_123_456_7;
    let mut small_intensity = peak_spectrum(1.);
    small_intensity.peaks[0].intensity = 0.123_456_79;
    let mut rt = peak_spectrum(0.123_456_789_012_345_67);
    rt.peaks[0].mz = 100.;
    for (e, line) in [
        (experiment(small_mz), "1.0\t1.234567890123457\t42.25\n"),
        (
            experiment(small_intensity),
            "1.0\t123.456789012344998\t0.123457\n",
        ),
        (experiment(rt), "0.123456789012346\t100.0\t42.25\n"),
    ] {
        let mut bytes = vec![9];
        let error = dta2d::write(&mut bytes, &e).unwrap_err();
        assert!(matches!(error, openms::Error::Unsupported(_)), "{error:?}");
        assert_eq!(bytes, [9]);
        let mut bytes = Vec::new();
        dta2d::write_with_options(&mut bytes, &e, &dta2d::WriteOptions::source()).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            format!("#SEC\tMZ\tINT\n{line}")
        );
    }
    // The TIC applies the same check to its retention times and sums.
    let mut s = peak_spectrum(1.);
    s.peaks[0].intensity = 0.123_456_79;
    let mut bytes = Vec::new();
    assert!(dta2d::write_tic(&mut bytes, &experiment(s.clone())).is_err());
    assert!(bytes.is_empty());
    dta2d::write_tic_with_options(&mut bytes, &experiment(s), &dta2d::WriteOptions::source())
        .unwrap();
    assert_eq!(bytes, b"#SEC\tMZ\tINT\n1.0\t0\t0.123457\n");
}

/// Source `DTA2DFile::store` writes every spectrum's peaks under its retention
/// time, whatever it carries; an empty spectrum writes nothing. The source mode
/// does the same, where the default refuses each of these.
#[test]
fn dta2d_source_mode_writes_what_the_source_store_writes() {
    let mut ms2_spectrum = peak_spectrum(2.);
    ms2_spectrum.ms_level = 2;
    ms2_spectrum.precursors.push(Precursor::new(500., 2));
    ms2_spectrum.native_id = "scan=7".into();
    ms2_spectrum.name = "named".into();
    let e = MSExperiment {
        spectra: vec![
            peak_spectrum(1.),
            MSSpectrum {
                rt: 1.5,
                ..Default::default()
            },
            ms2_spectrum,
            peak_spectrum(2.00001),
        ],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    assert!(dta2d::write(&mut bytes, &e).is_err());
    assert!(bytes.is_empty());
    dta2d::write_with_options(&mut bytes, &e, &dta2d::WriteOptions::source()).unwrap();
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        "#SEC\tMZ\tINT\n\
         1.0\t123.456789012344998\t42.25\n\
         2.0\t123.456789012344998\t42.25\n\
         2.00001\t123.456789012344998\t42.25\n"
    );
    // Nonfinite values and the ceilings are refused in both modes.
    let mut nan = peak_spectrum(1.);
    nan.peaks[0].intensity = f32::NAN;
    let mut bytes = Vec::new();
    assert!(
        dta2d::write_with_options(&mut bytes, &experiment(nan), &dta2d::WriteOptions::source())
            .is_err()
    );
    let limited = dta2d::WriteOptions {
        limits: ms2::Limits {
            max_output_bytes: 20,
            ..Default::default()
        },
        ..dta2d::WriteOptions::source()
    };
    assert!(dta2d::write_with_options(&mut bytes, &e, &limited).is_err());
    assert!(bytes.is_empty());
}

/// `DTA2DFile_test.cpp:72-73` and `MS2File_test.cpp:49-50`: a path load
/// records the document's path and its content-detected type; a stream read
/// records neither.
#[test]
fn path_loads_record_the_loaded_file_path_and_type() {
    let path = data_path("DTA2DFile_test_1.dta2d");
    for e in [
        dta2d::load(&path).unwrap(),
        dta2d::load_with_options(&path, &Default::default()).unwrap(),
    ] {
        assert_eq!(e.settings.document.loaded_file_path, path.to_str().unwrap());
        assert_eq!(e.settings.document.loaded_file_type.name(), "dta2d");
        assert_eq!(e.len(), 9);
    }
    let path = data_path("MS2File_test_spectra.ms2");
    for e in [
        ms2::load(&path).unwrap(),
        ms2::load_with_options(&path, &Default::default()).unwrap(),
    ] {
        assert_eq!(e.settings.document.loaded_file_path, path.to_str().unwrap());
        assert_eq!(e.settings.document.loaded_file_type.name(), "ms2");
        assert_eq!(e.len(), 2);
    }
    let read = dta2d::read(DTA[0]).unwrap();
    assert!(read.settings.document.loaded_file_path.is_empty());
    let read = ms2::read(MS2).unwrap();
    assert!(read.settings.document.loaded_file_path.is_empty());
    // A file that cannot be opened is an I/O error, the source's FileNotFound.
    let missing = data_path("missing.dta2d");
    assert!(matches!(dta2d::load(&missing), Err(openms::Error::Io(_))));
    assert!(matches!(
        ms2::load(data_path("missing.ms2")),
        Err(openms::Error::Io(_))
    ));
}

/// `TEST_REAL_SIMILAR` with the class tests' default tolerances: similar when
/// the relative difference is at most 1e-5 or the absolute one at most 1e-5.
fn similar(actual: f64, expected: f64) -> bool {
    let absolute = (actual - expected).abs();
    absolute <= 1e-5 || absolute <= 1e-5 * expected.abs().max(actual.abs())
}

fn temp_path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("openms-dta-{}-{name}", std::process::id()))
}

/// `DTAFile_test.cpp:26-265`, both halves (the source loads the same file into
/// a templated and a plain spectrum): 25 peaks, the precursor recomputed from
/// MH+ 1745.22 and charge 3, MS level 2, and the file's base name as the
/// spectrum's name (`DTAFile.h:164`).
#[test]
fn dta_class_test_load() {
    let path = data_path("DTAFile_test.dta");
    let s = dta::load(&path).unwrap();
    assert_eq!(s.name, "DTAFile_test.dta");
    assert_eq!(s.ms_level, dta::DEFAULT_MS_LEVEL);
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.precursors.len(), 1);
    assert!(
        similar(s.precursors[0].mz, 582.40666),
        "{}",
        s.precursors[0].mz
    );
    assert!(
        similar(s.precursors[0].mz, 582.4066),
        "{}",
        s.precursors[0].mz
    );
    assert_eq!(s.precursors[0].charge, 3);
    let expected = [
        (139.42, 318.52),
        (149.93, 61870.99),
        (169.65, 62074.22),
        (189.30, 53737.85),
        (202.28, 49410.25),
        (207.82, 17038.71),
        (219.72, 73629.98),
        (230.02, 47218.89),
        (231.51, 89935.22),
        (263.88, 81685.67),
        (318.38, 11876.49),
        (362.22, 91984.30),
        (389.84, 35049.17),
        (489.86, 55259.42),
        (508.47, 63411.68),
        (521.44, 42096.49),
        (546.98, 7133.92),
        (562.72, 47540.11),
        (579.61, 63350.29),
        (600.60, 29075.85),
        (628.64, 41030.62),
        (629.66, 51153.79),
        (630.37, 31411.63),
        (652.64, 26967.56),
        (712.18, 29235.72),
    ];
    assert_eq!(s.peaks.len(), 25);
    for (peak, (mz, intensity)) in s.peaks.iter().zip(expected) {
        assert!(similar(peak.mz, mz), "{} {mz}", peak.mz);
        assert!(
            similar(f64::from(peak.intensity), intensity),
            "{} {intensity}",
            peak.intensity
        );
    }
    // The stream form is the same spectrum without the name.
    let read = dta::read(std::io::BufReader::new(std::fs::File::open(&path).unwrap())).unwrap();
    assert!(read.name.is_empty());
    assert_eq!(read.peaks, s.peaks);
    assert_eq!(read.precursors, s.precursors);
    // `TEST_EXCEPTION(Exception::FileNotFound, f1.load("data_Idontexist", s))`.
    assert!(matches!(
        dta::load("data_Idontexist"),
        Err(openms::Error::Io(_))
    ));
}

/// `DTAFile_test.cpp:266-340`: a stored spectrum loads back with its
/// precursor and peaks, with and without a precursor, in the default and in
/// the source mode (whose legacy 1.0 Da proton mass moves the charge-3
/// precursor by 7.3 mDa, still within `TEST_REAL_SIMILAR`).
#[test]
fn dta_class_test_store() {
    let peaks = vec![
        Peak1D::new(11.4, 11.5),
        Peak1D::new(12.4, 12.5),
        Peak1D::new(13.4, 13.5),
    ];
    let with_precursor = MSSpectrum {
        precursors: vec![Precursor::new(582.40666, 3)],
        peaks: peaks.clone(),
        ..Default::default()
    };
    let without = MSSpectrum {
        peaks,
        ..Default::default()
    };
    let path = temp_path("store.dta");
    for options in [dta::WriteOptions::default(), dta::WriteOptions::source()] {
        dta::store_with_options(&path, &with_precursor, &options).unwrap();
        let loaded = dta::load(&path).unwrap();
        assert_eq!(loaded.precursors.len(), 1);
        assert!(similar(loaded.precursors[0].mz, 582.40666));
        assert_eq!(loaded.precursors[0].charge, 3);
        assert_eq!(loaded.peaks.len(), 3);
        for (peak, (mz, intensity)) in
            loaded
                .peaks
                .iter()
                .zip([(11.4, 11.5), (12.4, 12.5), (13.4, 13.5)])
        {
            assert!(similar(peak.mz, mz));
            assert!(similar(f64::from(peak.intensity), intensity));
        }
        dta::store_with_options(&path, &without, &options).unwrap();
        let loaded = dta::load(&path).unwrap();
        assert_eq!(loaded.peaks.len(), 3);
        // The source writes a default precursor, "0 0", for a spectrum without one.
        assert_eq!(loaded.precursors, [Precursor::new(0., 0)]);
    }
    // The default mode is the reader's exact inverse up to the text's precision.
    dta::store(&path, &with_precursor).unwrap();
    assert!((dta::load(&path).unwrap().precursors[0].mz - 582.40666).abs() < 1e-9);
    // A refused spectrum creates no file.
    std::fs::remove_file(&path).unwrap();
    let mut two = with_precursor.clone();
    two.precursors.push(Precursor::new(100., 1));
    assert!(dta::store(&path, &two).is_err());
    assert!(!path.exists());
    dta::store_with_options(&path, &two, &dta::WriteOptions::source()).unwrap();
    std::fs::remove_file(&path).unwrap();
}

/// `DTAFile.h:87-160`: a line is trimmed and split on a tab when it has one,
/// otherwise on single spaces, and must give exactly two fields; the number
/// conversion (`StringUtils::toDouble`/`toInt32`) skips surrounding spaces,
/// tabs, CR and LF.
#[test]
fn dta_lines_split_as_the_source_splits_them() {
    let ok = |text: &str| dta::read(text.as_bytes()).unwrap();
    let s = ok("1745.22\t 3\r\n139.42 \t318.52\n\n  \t \n149.93\t61870.99 \n");
    assert_eq!(s.precursors[0].charge, 3);
    assert_eq!(
        s.peaks,
        [Peak1D::new(139.42, 318.52), Peak1D::new(149.93, 61870.99)]
    );
    assert_eq!(
        ok("500 +2\n").precursors,
        [Precursor::new(
            (500. - openms::chemistry::PROTON_MASS_U) / 2. + openms::chemistry::PROTON_MASS_U,
            2
        )]
    );
    assert_eq!(ok("500 0\n").precursors, [Precursor::new(500., 0)]);
    for bad in [
        "",
        "\n1 2\n",
        "500  2\n",
        "500 2 1\n",
        "500 2.5\n",
        "500 2\n1  2\n",
        "500 2\n1 2 3\n",
        "500 2\n1\n",
        "500 2\n1 x\n",
        "500 2\n1\t\t2\n",
    ] {
        let error = dta::read(bad.as_bytes()).unwrap_err();
        assert!(
            matches!(error, openms::Error::Parse { .. }),
            "{bad:?}: {error:?}"
        );
    }
}

/// `DTA2DFile::getOptions() = options`: a source `PeakFileOptions` gives the
/// reader its three consumed ranges (`DTA2DFile.h:211-240`) and nothing else.
#[test]
fn dta2d_read_options_from_peak_file_options() {
    use openms::format::PeakFileOptions;
    use openms::kernel::NumericRange;
    let none = dta2d::ReadOptions::from_peak_options(&PeakFileOptions::default());
    assert_eq!(none, dta2d::ReadOptions::default());
    let mut options = PeakFileOptions::default();
    options.set_rt_range(NumericRange {
        min: 4711.15,
        max: 4711.45,
    });
    options.set_mz_range(NumericRange {
        min: 150.,
        max: 250.,
    });
    options.set_intensity_range(NumericRange {
        min: 50000.,
        max: 90000.,
    });
    options.set_ms_levels(&[2]).unwrap();
    options.metadata_only = true;
    let converted = dta2d::ReadOptions::from_peak_options(&options);
    assert_eq!(
        converted,
        dta2d::ReadOptions {
            rt_range: Some(4711.15..4711.45),
            mz_range: Some(150.0..250.0),
            intensity_range: Some(50000.0..90000.0),
            ..Default::default()
        }
    );
    // The MS-level and metadata-only settings have no DTA2D behaviour.
    let e = dta2d::read_with_options(DTA[0], &converted).unwrap();
    assert!(
        e.spectra
            .iter()
            .all(|s| s.ms_level == 1 && !s.peaks.is_empty())
    );
}

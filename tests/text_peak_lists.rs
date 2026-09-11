// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::format::{dta2d, ms2};
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
            &ms2::Limits {
                max_output_bytes: 2,
                ..Default::default()
            }
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

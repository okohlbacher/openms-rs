//! `BrukerTimsFile` and the opentims layer under it (`TimsDataHandle`) on
//! synthetic timsTOF `.d` directories written by this test.
//!
//! The upstream class test `BrukerTimsFile_test.cpp` (core pin bc9cc12) runs
//! its loading sections only against two externally downloaded datasets
//! (`OPENTIMS_DDA_TEST_DATA`, `OPENTIMS_DIA_TEST_DATA`), and the Release
//! install is built without opentims, so there is neither retained C++ output
//! nor a C++ build to execute. The section that needs no data, loading a
//! missing path, is ported as `missing_directory_is_not_readable` (tier 3).
//! Every other expectation here is derived independently of the port (tier 4):
//! the synthetic frames are encoded by `tests/support/tdf_fixture.rs` from the TDF layout that
//! opentims decodes (`TimsFrame::save_to_buffs`), and the expected m/z, 1/K0
//! and intensities are computed from the opentims formulas with values chosen
//! so the arithmetic is exact or stated.
#![cfg(feature = "bruker-tdf")]

use openms::Error;
use openms::format::bruker_tims_file::{
    BrukerTimsFile, CentroidAlgo, Config, ExportMode, IM_ARRAY_NAME, TimsCalibrationStrategy,
};
use openms::format::tims_data_handle::{TdfLimits, TimsDataHandle};

#[path = "support/tdf_fixture.rs"]
mod fixture;
use fixture::{Frame, dataset, im_of, mz_of, two_frames};

fn frame_config() -> Config {
    Config {
        export_mode: ExportMode::Frame,
        ..Config::default()
    }
}

/// Class-test section `load(path, exp, config)`: a missing path is
/// `Exception::FileNotReadable`, here an I/O error of kind `NotFound`.
#[test]
fn missing_directory_is_not_readable() {
    match BrukerTimsFile::new().load("/nonexistent/path.d", &Config::default()) {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("{other:?}"),
    }
}

#[test]
fn handle_decodes_frames_as_opentims() {
    let data = dataset("handle", &two_frames(), "");
    let handle = TimsDataHandle::open(&data.dir, TdfLimits::default()).unwrap();
    assert_eq!(handle.min_frame_id(), 1);
    assert_eq!(handle.max_frame_id(), 2);
    let ms1 = handle.read_frame(1).unwrap();
    assert_eq!(ms1.scan_ids, vec![0, 0, 3, 7, 7]);
    assert_eq!(ms1.tofs, vec![500, 900, 100, 250, 251]);
    // AccumulationTime 100: correction 1, intensities unchanged.
    assert_eq!(ms1.intensities, vec![7, 3, 11, 5, 6]);
    // AccumulationTime 50: correction 2, `(uint32_t)(i * 2.0 + 0.5)`.
    assert_eq!(handle.read_frame(2).unwrap().intensities, vec![40, 80]);
    assert_eq!(handle.tof2mz().convert(&[500]), vec![mz_of(500)]);
}

#[test]
fn frame_export_loads_every_frame() {
    let data = dataset("frame", &two_frames(), "");
    let exp = BrukerTimsFile::new()
        .load(&data.dir, &frame_config())
        .unwrap();
    assert_eq!(exp.spectra.len(), 2);

    let ms1 = &exp.spectra[0];
    assert_eq!(ms1.rt, 10.0);
    assert_eq!(ms1.ms_level, 1);
    assert_eq!(ms1.native_id, "frame=1");
    assert_eq!(ms1.spectrum_type, openms::kernel::SpectrumType::Centroid);
    assert_eq!(
        ms1.drift_time_unit,
        openms::metadata::DriftTimeUnit::InverseReducedMobility
    );
    // Sorted by m/z at the end of load (the source's `sortSpectra(true)`),
    // each 1/K0 moving with its peak.
    let expected: Vec<(u32, u32, u32)> = vec![
        (100, 11, 3),
        (250, 5, 7),
        (251, 6, 7),
        (500, 7, 0),
        (900, 3, 0),
    ];
    let mzs: Vec<f64> = ms1.peaks.iter().map(|p| p.mz).collect();
    assert_eq!(
        mzs,
        expected
            .iter()
            .map(|&(tof, _, _)| mz_of(tof))
            .collect::<Vec<_>>()
    );
    let intensities: Vec<f32> = ms1.peaks.iter().map(|p| p.intensity).collect();
    assert_eq!(
        intensities,
        expected
            .iter()
            .map(|&(_, i, _)| i as f32)
            .collect::<Vec<_>>()
    );
    assert_eq!(ms1.float_data_arrays.len(), 1);
    assert_eq!(ms1.float_data_arrays[0].name, IM_ARRAY_NAME);
    assert_eq!(
        ms1.float_data_arrays[0].data,
        expected
            .iter()
            .map(|&(_, _, scan)| im_of(scan, 8) as f32)
            .collect::<Vec<_>>()
    );

    let ms2 = &exp.spectra[1];
    assert_eq!(ms2.ms_level, 2);
    assert_eq!(ms2.native_id, "frame=2");
    assert_eq!(ms2.peaks.len(), 2);

    let source = &exp.settings.source_files;
    assert_eq!(source.len(), 1);
    assert_eq!(source[0].file_type, "Bruker TDF");
    assert_eq!(source[0].native_id_type, "Bruker TDF nativeID format");
    assert_eq!(source[0].native_id_type_accession, "MS:1002818");
    assert!(data.dir.to_string_lossy().ends_with(&source[0].name));
}

#[test]
fn frame_and_rt_ranges_select_frames() {
    let data = dataset("ranges", &two_frames(), "");
    let file = BrukerTimsFile::new();
    let only_second = Config {
        frame_id_min: 2,
        ..frame_config()
    };
    let exp = file.load(&data.dir, &only_second).unwrap();
    assert_eq!(exp.spectra.len(), 1);
    assert_eq!(exp.spectra[0].native_id, "frame=2");

    let no_ms1 = Config {
        load_ms1: false,
        ..frame_config()
    };
    assert_eq!(file.load(&data.dir, &no_ms1).unwrap().spectra.len(), 1);

    let rt = Config {
        rt_max_sec: 10.2,
        ..frame_config()
    };
    let exp = file.load(&data.dir, &rt).unwrap();
    assert_eq!(exp.spectra.len(), 1);
    assert_eq!(exp.spectra[0].native_id, "frame=1");

    for empty in [
        Config {
            rt_min_sec: 20.0,
            ..frame_config()
        },
        Config {
            rt_min_sec: 11.0,
            rt_max_sec: 10.0,
            ..frame_config()
        },
        Config {
            frame_id_min: 2,
            frame_id_max: 1,
            ..frame_config()
        },
    ] {
        assert!(file.load(&data.dir, &empty).unwrap().spectra.is_empty());
    }
}

#[test]
fn rational_calibration_is_preferred_when_present() {
    // ModelType-2 coefficients of the class test; 1/K0 of scan 3 by the
    // source formula, computed here from the coefficients.
    let (c0, c1, c2, c3, c4, c6, c7) =
        (1.0, 917.0, 213.5998, 75.81729, 33.0, -0.009065829, 135.4364);
    let calibration = format!(
        "CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, C0 REAL, C1 REAL,
           C2 REAL, C3 REAL, C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL);
         INSERT INTO TimsCalibration VALUES (1, 2, {c0:?}, {c1:?}, {c2:?}, {c3:?}, {c4:?}, 1.0,
           {c6:?}, {c7:?}, 13.32608, 1663.341);"
    );
    let data = dataset("rational", &two_frames(), &calibration);
    let voltage = c2 + ((c3 - c2) / c1) * (3.0 - c4 - c0);
    let rational = (1.0 / (c6 + c7 / voltage)) as f32;
    let file = BrukerTimsFile::new();

    let exp = file.load(&data.dir, &frame_config()).unwrap();
    // Peak (tof 100, scan 3) is the smallest m/z of frame 1.
    assert_eq!(exp.spectra[0].float_data_arrays[0].data[0], rational);

    let linear = Config {
        tims_calibration_strategy: TimsCalibrationStrategy::Linear,
        ..frame_config()
    };
    let exp = file.load(&data.dir, &linear).unwrap();
    assert_eq!(
        exp.spectra[0].float_data_arrays[0].data[0],
        im_of(3, 8) as f32
    );
}

#[test]
fn calibration_strategies_that_cannot_be_met_are_refused() {
    let data = dataset("strategies", &two_frames(), "");
    let file = BrukerTimsFile::new();
    let rational = Config {
        tims_calibration_strategy: TimsCalibrationStrategy::Rational,
        ..frame_config()
    };
    assert!(matches!(file.load(&data.dir, &rational), Err(Error::Io(_))));
    let sdk = Config {
        tims_calibration_strategy: TimsCalibrationStrategy::BrukerSdk,
        ..frame_config()
    };
    assert!(matches!(
        file.load(&data.dir, &sdk),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn dda_export_is_not_ported() {
    let data = dataset("dda", &two_frames(), "");
    let file = BrukerTimsFile::new();
    assert!(matches!(
        file.load(&data.dir, &Config::default()),
        Err(Error::Unsupported(_))
    ));
    let spectrum = Config {
        export_mode: ExportMode::Spectrum,
        ..Config::default()
    };
    assert!(matches!(
        file.load(&data.dir, &spectrum),
        Err(Error::Unsupported(_))
    ));
}

const DIA_TABLES: &str = "
    CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER,
      ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);
    INSERT INTO DiaFrameMsMsWindows VALUES (1, 0, 3, 500.0, 25.0, 30.0), (1, 4, 7, 700.0, 25.0, 30.0);
    CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER PRIMARY KEY, WindowGroup INTEGER);
    INSERT INTO DiaFrameMsMsInfo VALUES (2, 1);";

#[test]
fn dia_export_splits_frames_by_window() {
    let data = dataset("dia", &two_frames(), DIA_TABLES);
    let exp = BrukerTimsFile::new()
        .load(&data.dir, &Config::default())
        .unwrap();
    // MS1 frame 1, then frame 2 once per window of its group (equal RT, kept
    // in window order by the stable sort).
    assert_eq!(exp.spectra.len(), 3);
    assert_eq!(exp.spectra[0].native_id, "frame=1");
    let first = &exp.spectra[1];
    assert_eq!(first.native_id, "frame=2 windowGroup=1 scan=0");
    assert_eq!(first.ms_level, 2);
    assert_eq!(first.precursors.len(), 1);
    assert_eq!(first.precursors[0].mz, 500.0);
    assert_eq!(first.precursors[0].isolation_window_lower_offset, 12.5);
    assert_eq!(first.precursors[0].isolation_window_upper_offset, 12.5);
    // Window scans 0..=3: 1/K0 in [im(3), im(0)]; frame 2's peak at scan 2.
    assert_eq!(first.peaks.len(), 1);
    assert_eq!(first.peaks[0].mz, mz_of(400));
    assert_eq!(first.peaks[0].intensity, 40.0);
    assert_eq!(
        first.metadata["ion mobility lower limit"].to_string(),
        im_of(3, 8).to_string()
    );
    let second = &exp.spectra[2];
    assert_eq!(second.native_id, "frame=2 windowGroup=1 scan=4");
    assert_eq!(second.peaks.len(), 1);
    assert_eq!(second.peaks[0].mz, mz_of(300));
    assert_eq!(second.float_data_arrays[0].data, vec![im_of(5, 8) as f32]);
}

#[test]
fn dia_aggregation_is_not_ported() {
    let data = dataset("dia-agg", &two_frames(), DIA_TABLES);
    let aggregated = Config {
        dia_ms2_n_neighbors: 1,
        ..Config::default()
    };
    assert!(matches!(
        BrukerTimsFile::new().load(&data.dir, &aggregated),
        Err(Error::Unsupported(_))
    ));
}

/// One ion in three consecutive scans at the same TOF, intensities 10, 30,
/// 10: one hill, summed intensity 50, apex at the middle scan.
#[test]
fn hill_based_ms1_centroiding() {
    let mut scans = vec![Vec::new(); 8];
    scans[2] = vec![(600, 10)];
    scans[3] = vec![(600, 30)];
    scans[4] = vec![(600, 10)];
    let frames = vec![Frame {
        id: 1,
        time: 1.0,
        msms_type: 0,
        accumulation_time: 100.0,
        scans,
    }];
    let data = dataset("hill", &frames, "");
    let config = Config {
        ms1_centroid_algo: CentroidAlgo::HillBased,
        expose_hill_bounds: true,
        ..frame_config()
    };
    let exp = BrukerTimsFile::new().load(&data.dir, &config).unwrap();
    let spectrum = &exp.spectra[0];
    assert_eq!(spectrum.peaks.len(), 1);
    assert_eq!(spectrum.peaks[0].mz, mz_of(600));
    assert_eq!(spectrum.peaks[0].intensity, 50.0);
    let names: Vec<&str> = spectrum
        .float_data_arrays
        .iter()
        .map(|a| a.name.as_str())
        .collect();
    assert_eq!(
        names,
        [
            IM_ARRAY_NAME,
            "im lower bound",
            "im upper bound",
            "m/z lower bound",
            "m/z upper bound"
        ]
    );
    assert_eq!(spectrum.float_data_arrays[0].data, vec![im_of(3, 8) as f32]);
    assert_eq!(spectrum.float_data_arrays[1].data, vec![im_of(4, 8) as f32]);
    assert_eq!(spectrum.float_data_arrays[2].data, vec![im_of(2, 8) as f32]);
}

#[test]
fn greedy2d_is_not_ported() {
    let data = dataset("greedy", &two_frames(), "");
    let config = Config {
        ms1_centroid_algo: CentroidAlgo::Greedy2d,
        ms1_centroid_im_pct: 3.0,
        ..frame_config()
    };
    assert!(matches!(
        BrukerTimsFile::new().load(&data.dir, &config),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn corrupt_frames_are_refused() {
    let data = dataset("corrupt", &two_frames(), "");
    let bin = data.dir.join("analysis.tdf_bin");
    let mut bytes = std::fs::read(&bin).unwrap();
    // Frame 1's scan count no longer matches Frames.NumScans.
    bytes[4] = 9;
    std::fs::write(&bin, &bytes).unwrap();
    let handle = TimsDataHandle::open(&data.dir, TdfLimits::default()).unwrap();
    assert!(matches!(handle.read_frame(1), Err(Error::Parse { .. })));
    bytes[4] = 8;
    // Garbage instead of a zstd frame.
    for byte in &mut bytes[8..16] {
        *byte = 0xAB;
    }
    std::fs::write(&bin, &bytes).unwrap();
    let handle = TimsDataHandle::open(&data.dir, TdfLimits::default()).unwrap();
    assert!(matches!(handle.read_frame(1), Err(Error::Parse { .. })));

    // A frame over the peak ceiling is refused before allocation.
    let small = TdfLimits {
        max_frame_peaks: 4,
        ..TdfLimits::default()
    };
    let data = dataset("ceiling", &two_frames(), "");
    let handle = TimsDataHandle::open(&data.dir, small).unwrap();
    assert!(matches!(handle.read_frame(1), Err(Error::InvalidValue(_))));
    assert_eq!(handle.read_frame(2).unwrap().tofs.len(), 2);
}

#[test]
fn metadata_problems_are_refused() {
    let unsupported = dataset(
        "compression",
        &two_frames(),
        "UPDATE GlobalMetadata SET Value = '1' WHERE Key = 'TimsCompressionType';",
    );
    assert!(matches!(
        TimsDataHandle::open(&unsupported.dir, TdfLimits::default()),
        Err(Error::Unsupported(_))
    ));
    let null = dataset(
        "null",
        &two_frames(),
        "UPDATE Frames SET AccumulationTime = NULL WHERE Id = 2;",
    );
    assert!(matches!(
        TimsDataHandle::open(&null.dir, TdfLimits::default()),
        Err(Error::Parse { .. })
    ));
    let no_range = dataset(
        "range",
        &two_frames(),
        "DELETE FROM GlobalMetadata WHERE Key = 'MzAcqRangeUpper';",
    );
    assert!(matches!(
        TimsDataHandle::open(&no_range.dir, TdfLimits::default()),
        Err(Error::InvalidValue(_))
    ));
}

/// `transform` announces the SQL frame counts and feeds the sorted spectra.
#[test]
fn transform_feeds_a_consumer() {
    use openms::interfaces::MSDataConsumer;
    use openms::metadata::ExperimentalSettings;
    use openms::{MSChromatogram, MSSpectrum};
    use std::ops::ControlFlow;

    #[derive(Default)]
    struct Collect {
        expected: usize,
        files: usize,
        ids: Vec<String>,
    }
    impl MSDataConsumer for Collect {
        fn set_expected_size(&mut self, spectra: usize, _: usize) -> openms::Result<()> {
            self.expected = spectra;
            Ok(())
        }
        fn set_experimental_settings(
            &mut self,
            settings: &ExperimentalSettings,
        ) -> openms::Result<()> {
            self.files = settings.source_files.len();
            Ok(())
        }
        fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> openms::Result<ControlFlow<()>> {
            self.ids.push(s.native_id.clone());
            Ok(ControlFlow::Continue(()))
        }
        fn consume_chromatogram(
            &mut self,
            _: &mut MSChromatogram,
        ) -> openms::Result<ControlFlow<()>> {
            Ok(ControlFlow::Continue(()))
        }
    }
    let data = dataset("transform", &two_frames(), DIA_TABLES);
    let mut consumer = Collect::default();
    BrukerTimsFile::new()
        .transform(&data.dir, &mut consumer, &Config::default())
        .unwrap();
    assert_eq!(consumer.expected, 3);
    assert_eq!(consumer.files, 1);
    assert_eq!(
        consumer.ids,
        [
            "frame=1",
            "frame=2 windowGroup=1 scan=0",
            "frame=2 windowGroup=1 scan=4"
        ]
    );
}

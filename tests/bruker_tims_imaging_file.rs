//! Replication of the upstream class test `BrukerTimsImagingFile_test.cpp`
//! (OpenMS core `bc9cc12`, `src/tests/class_tests/openms/source/`), plus
//! end-to-end loads the class test does not have.
//!
//! Every section of the class test is ported (tier 3): its fake `.d` folders
//! hold only an `analysis.tdf` with `GlobalMetadata` and `MaldiFrameInfo`,
//! built here by `fake_d` exactly as the test's `createFakeD` builds them, and
//! its `load` sections exercise only the error paths. The Release install on
//! ibminode06 is built without opentims, so no C++ run can supply more.
//!
//! The positive loads are tier 4: the frames are synthesised by
//! `tests/support/tdf_fixture.rs`, and the expected grid, offsets, pixel size
//! and pixel-to-spectrum binding follow from the source's `load` (subtract the
//! per-axis minimum, bounding-box extent, `MaldiSpotSize_um`, spectra bound by
//! their `frame=<id>` native ID).
#![cfg(all(feature = "bruker-tdf", feature = "mzml"))]

use openms::Error;
use openms::format::bruker_tims_file::{self, ExportMode};
use openms::format::bruker_tims_imaging_file::{BrukerTimsImagingFile, Config, MaldiPixel};
use rusqlite::Connection;
use std::path::{Path, PathBuf};

#[path = "support/tdf_fixture.rs"]
mod fixture;
use fixture::{Frame, dataset};

/// A temporary directory, removed on drop (the class test's `File::TempDir`).
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("openms-tims-img-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The class test's `createFakeD`: `<parent>/<name>.d/analysis.tdf` with
/// `GlobalMetadata` (`MaldiApplicationType`, optionally `MaldiSpotSize_um`)
/// and `MaldiFrameInfo` rows `(Frame, X, Y[, Z])`.
fn fake_d(
    parent: &Path,
    name: &str,
    application_type: &str,
    rows: &[(i32, i32, i32, i32)],
    with_z_column: bool,
    spot_size_um: &str,
) -> PathBuf {
    let d = parent.join(format!("{name}.d"));
    std::fs::create_dir_all(&d).unwrap();
    let db = Connection::open(d.join("analysis.tdf")).unwrap();
    db.execute_batch("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);")
        .unwrap();
    if !application_type.is_empty() {
        db.execute(
            "INSERT INTO GlobalMetadata (Key, Value) VALUES ('MaldiApplicationType', ?1)",
            [application_type],
        )
        .unwrap();
    }
    if !spot_size_um.is_empty() {
        db.execute(
            "INSERT INTO GlobalMetadata (Key, Value) VALUES ('MaldiSpotSize_um', ?1)",
            [spot_size_um],
        )
        .unwrap();
    }
    if with_z_column {
        db.execute_batch(
            "CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, \
             YIndexPos INTEGER, ZIndexPos INTEGER);",
        )
        .unwrap();
    } else {
        db.execute_batch(
            "CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER);",
        )
        .unwrap();
    }
    for &(frame, x, y, z) in rows {
        if with_z_column {
            db.execute(
                "INSERT INTO MaldiFrameInfo VALUES (?1, ?2, ?3, ?4)",
                [frame, x, y, z],
            )
            .unwrap();
        } else {
            db.execute(
                "INSERT INTO MaldiFrameInfo VALUES (?1, ?2, ?3)",
                [frame, x, y],
            )
            .unwrap();
        }
    }
    d
}

/// The source's `Exception::FileNotReadable`: an I/O error of kind `NotFound`.
fn assert_not_readable<T: std::fmt::Debug>(result: Result<T, Error>) {
    match result {
        Err(Error::Io(error)) => assert_eq!(error.kind(), std::io::ErrorKind::NotFound),
        other => panic!("expected FileNotReadable, got {other:?}"),
    }
}

/// Section `static bool isImagingDataset(const String& path)`.
#[test]
fn is_imaging_dataset() {
    let tmp = TempDir::new("probe");
    let imaging = fake_d(
        &tmp.0,
        "imaging",
        "Imaging",
        &[(1, 0, 0, 0), (2, 1, 0, 0), (3, 0, 1, 0), (4, 1, 1, 0)],
        false,
        "",
    );
    assert!(BrukerTimsImagingFile::is_imaging_dataset(&imaging));
    let droplet = fake_d(
        &tmp.0,
        "droplet",
        "SingleSpectra",
        &[(1, 0, 0, 0)],
        false,
        "",
    );
    assert!(!BrukerTimsImagingFile::is_imaging_dataset(&droplet));
    assert!(!BrukerTimsImagingFile::is_imaging_dataset(
        tmp.0.join("does_not_exist.d")
    ));
}

/// Section `static String readGlobalMetadataValue(const String& d_folder, const String& key)`.
#[test]
fn read_global_metadata_value() {
    let tmp = TempDir::new("gm");
    let d = fake_d(&tmp.0, "gm", "Imaging", &[(1, 0, 0, 0)], false, "50");
    assert_eq!(
        BrukerTimsImagingFile::read_global_metadata_value(&d, "MaldiApplicationType").unwrap(),
        "Imaging"
    );
    assert_eq!(
        BrukerTimsImagingFile::read_global_metadata_value(&d, "MaldiSpotSize_um").unwrap(),
        "50"
    );
    assert_eq!(
        BrukerTimsImagingFile::read_global_metadata_value(&d, "UnknownKey").unwrap(),
        ""
    );
    assert_not_readable(BrukerTimsImagingFile::read_global_metadata_value(
        tmp.0.join("nope.d"),
        "MaldiApplicationType",
    ));
}

/// Section `static std::vector<MaldiPixel> readMaldiFrameInfo(const String& d_folder)`.
#[test]
fn read_maldi_frame_info() {
    let tmp = TempDir::new("mfi");
    let d = fake_d(
        &tmp.0,
        "grid2x2",
        "Imaging",
        &[
            (10, 10, 5, 0),
            (11, 11, 5, 0),
            (12, 10, 6, 0),
            (13, 11, 6, 0),
        ],
        false,
        "",
    );
    let rows = BrukerTimsImagingFile::read_maldi_frame_info(&d).unwrap();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0].frame_id, 10);
    assert_eq!(rows[0].x, 10);
    assert_eq!(rows[0].y, 5);
    assert_eq!(rows[3].frame_id, 13);
    assert_eq!(rows[3].x, 11);
    assert_eq!(rows[3].y, 6);

    let with_z = fake_d(
        &tmp.0,
        "withz",
        "Imaging",
        &[(1, 3, 4, 0), (2, 3, 5, 0)],
        true,
        "",
    );
    let rows_z = BrukerTimsImagingFile::read_maldi_frame_info(&with_z).unwrap();
    assert_eq!(rows_z.len(), 2);
    assert_eq!(rows_z[0].x, 3);
    assert_eq!(rows_z[0].y, 4);
    assert_eq!(rows_z[1].y, 5);

    let empty = fake_d(&tmp.0, "empty", "Imaging", &[], false, "");
    assert!(matches!(
        BrukerTimsImagingFile::read_maldi_frame_info(&empty),
        Err(Error::Parse { .. })
    ));
}

/// Section `void load(const String& path, MSImagingExperiment& exp)`: a
/// missing directory, and a directory without `analysis.tdf`, are both
/// `FileNotReadable`, never masked as `InvalidValue` by the imaging check.
#[test]
fn load_missing_inputs_are_not_readable() {
    let tmp = TempDir::new("load");
    let file = BrukerTimsImagingFile::new();
    assert_not_readable(file.load(tmp.0.join("nope.d"), &Config::default()));
    let no_tdf = tmp.0.join("no_tdf.d");
    std::fs::create_dir_all(&no_tdf).unwrap();
    assert_not_readable(file.load(&no_tdf, &Config::default()));
}

/// Section `[SingleSpectra] strict_imaging_only rejects dried-droplet datasets`.
#[test]
fn strict_imaging_only_rejects_dried_droplet() {
    let tmp = TempDir::new("rej");
    let d = fake_d(
        &tmp.0,
        "rej",
        "SingleSpectra",
        &[(1, 0, 0, 0), (2, 1, 0, 0)],
        false,
        "",
    );
    assert!(matches!(
        BrukerTimsImagingFile::new().load(&d, &Config::default()),
        Err(Error::InvalidValue(_))
    ));
}

/// Section `[Multi-section] rejects datasets with multiple distinct ZIndexPos`.
#[test]
fn multi_section_is_rejected() {
    let tmp = TempDir::new("multiz");
    let d = fake_d(
        &tmp.0,
        "multiz",
        "Imaging",
        &[(1, 0, 0, 0), (2, 1, 0, 0), (3, 0, 0, 1), (4, 1, 0, 1)],
        true,
        "",
    );
    assert!(matches!(
        BrukerTimsImagingFile::new().load(&d, &Config::default()),
        Err(Error::InvalidValue(_))
    ));
}

/// Four MS1 frames, Ids 1-4, each with one peak whose raw intensity is ten
/// times the frame Id (AccumulationTime 100: no correction).
fn imaging_frames() -> Vec<Frame> {
    (1..=4u32)
        .map(|id| {
            let mut scans = vec![Vec::new(); 4];
            scans[1] = vec![(100 * id, 10 * id)];
            Frame {
                id,
                time: f64::from(id),
                msms_type: 0,
                accumulation_time: 100.0,
                scans,
            }
        })
        .collect()
}

/// Tier 4, not in the class test: a grid at raw offset (10, 5) binds each
/// pixel, normalised to zero, to its frame's spectrum; `MaldiSpotSize_um` is
/// the pixel size; a pixel whose frame does not exist is dropped without
/// failing the load but still widens the bounding box; the inner export mode
/// and `load_ms1` are forced to `Frame` and true.
#[test]
fn load_binds_pixels_to_frame_spectra() {
    let data = dataset(
        "imaging-load",
        &imaging_frames(),
        "INSERT INTO GlobalMetadata VALUES ('MaldiApplicationType', 'Imaging'),
           ('MaldiSpotSize_um', '20');
         CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER,
           ZIndexPos INTEGER);
         INSERT INTO MaldiFrameInfo VALUES (4, 11, 6, 0), (1, 10, 5, 0), (2, 11, 5, 0),
           (3, 10, 6, 0), (9, 12, 7, 0);",
    );
    let config = Config {
        inner_config: bruker_tims_file::Config {
            export_mode: ExportMode::Spectrum,
            load_ms1: false,
            ..bruker_tims_file::Config::default()
        },
        strict_imaging_only: true,
    };
    let exp = BrukerTimsImagingFile::new()
        .load(&data.dir, &config)
        .unwrap();

    assert_eq!(exp.number_of_spectra(), 4);
    assert_eq!(exp.number_of_pixels(), 4);
    let geometry = exp.geometry();
    assert_eq!((geometry.width(), geometry.height()), (3, 3));
    assert_eq!(geometry.pixel_size_x(), 20.0);
    assert_eq!(geometry.pixel_size_y(), 20.0);
    assert_eq!(geometry.pixel_size_unit(), "micrometer");
    for (x, y, frame) in [(0, 0, 1u32), (1, 0, 2), (0, 1, 3), (1, 1, 4)] {
        let spectrum = exp.spectrum(x, y).unwrap();
        assert_eq!(spectrum.native_id, format!("frame={frame}"));
        assert_eq!(spectrum.peaks.len(), 1);
        assert_eq!(
            f64::from(spectrum.peaks[0].intensity),
            f64::from(10 * frame)
        );
    }
    assert!(!exp.has_pixel(2, 2));

    // readMaldiFrameInfo orders by frame and keeps raw coordinates.
    let rows = BrukerTimsImagingFile::read_maldi_frame_info(&data.dir).unwrap();
    let expected: Vec<MaldiPixel> = [(1, 10, 5), (2, 11, 5), (3, 10, 6), (4, 11, 6), (9, 12, 7)]
        .iter()
        .map(|&(frame_id, x, y)| MaldiPixel { frame_id, x, y })
        .collect();
    assert_eq!(rows, expected);
}

/// Tier 4: a restricted frame range loads fewer spectra, and the pixels of the
/// frames left out are dropped (the header's partial-result `@note`); without
/// a spot size in `GlobalMetadata` the pixel size is 1.
#[test]
fn restricted_frame_range_drops_pixels() {
    let data = dataset(
        "imaging-range",
        &imaging_frames(),
        "INSERT INTO GlobalMetadata VALUES ('MaldiApplicationType', 'Imaging');
         CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER);
         INSERT INTO MaldiFrameInfo VALUES (1, 0, 0), (2, 1, 0), (3, 0, 1), (4, 1, 1);",
    );
    let config = Config {
        inner_config: bruker_tims_file::Config {
            frame_id_min: 3,
            ..bruker_tims_file::Config::default()
        },
        strict_imaging_only: true,
    };
    let exp = BrukerTimsImagingFile::new()
        .load(&data.dir, &config)
        .unwrap();
    assert_eq!(exp.number_of_spectra(), 2);
    assert_eq!(exp.number_of_pixels(), 2);
    assert!(!exp.has_pixel(0, 0));
    assert_eq!(exp.spectrum(0, 1).unwrap().native_id, "frame=3");
    assert_eq!(exp.spectrum(1, 1).unwrap().native_id, "frame=4");
    assert_eq!(exp.geometry().pixel_size_x(), 1.0);
}

/// Tier 4: with `strict_imaging_only` off a non-imaging dataset gets past the
/// imaging check, and two rows at one coordinate are then `InvalidValue`, as
/// the source's `MSImagingGeometry::addPixel` throws.
#[test]
fn lenient_mode_reaches_duplicate_coordinate_check() {
    let data = dataset(
        "imaging-dup",
        &imaging_frames(),
        "INSERT INTO GlobalMetadata VALUES ('MaldiApplicationType', 'SingleSpectra');
         CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER);
         INSERT INTO MaldiFrameInfo VALUES (1, 0, 0), (2, 1, 0);",
    );
    let lenient = Config {
        strict_imaging_only: false,
        ..Config::default()
    };
    let exp = BrukerTimsImagingFile::new()
        .load(&data.dir, &lenient)
        .unwrap();
    assert_eq!(exp.number_of_pixels(), 2);

    let duplicate = dataset(
        "imaging-dup2",
        &imaging_frames(),
        "INSERT INTO GlobalMetadata VALUES ('MaldiApplicationType', 'Imaging');
         CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER);
         INSERT INTO MaldiFrameInfo VALUES (1, 0, 0), (2, 0, 0);",
    );
    assert!(matches!(
        BrukerTimsImagingFile::new().load(&duplicate.dir, &Config::default()),
        Err(Error::InvalidValue(_))
    ));
}

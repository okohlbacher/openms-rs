//! Synthetic timsTOF `.d` directories for the Bruker TDF tests: frames encoded
//! exactly as opentims decodes them (`TimsFrame::save_to_buffs`), with a
//! minimal `analysis.tdf`. Shared through `#[path]` by
//! `tests/bruker_tims_file.rs` and `tests/bruker_tims_imaging_file.rs`.
#![allow(dead_code)]

use rusqlite::Connection;
use std::path::PathBuf;

/// One synthetic frame: its metadata and the `(tof, raw intensity)` peaks of
/// each scan, TOF ascending.
pub struct Frame {
    pub id: u32,
    pub time: f64,
    pub msms_type: u32,
    pub accumulation_time: f64,
    pub scans: Vec<Vec<(u32, u32)>>,
}

/// Encode one frame packet exactly as opentims decodes it: header words with
/// twice the per-scan peak count at `s + 1`, TOF deltas from `u32::MAX`, four
/// byte planes, zstd, then the `(packet size, scan count)` header.
pub fn encode(frame: &Frame) -> Vec<u8> {
    let num_scans = frame.scans.len();
    let mut words = vec![0u32; num_scans];
    for s in 0..num_scans.saturating_sub(1) {
        words[s + 1] = 2 * frame.scans[s].len() as u32;
    }
    for scan in &frame.scans {
        let mut previous = u32::MAX;
        for &(tof, intensity) in scan {
            words.push(tof.wrapping_sub(previous));
            words.push(intensity);
            previous = tof;
        }
    }
    let mut planes = Vec::with_capacity(words.len() * 4);
    for byte in 0..4 {
        planes.extend(words.iter().map(|w| (w >> (8 * byte)) as u8));
    }
    let compressed =
        ruzstd::encoding::compress_to_vec(&planes[..], ruzstd::encoding::CompressionLevel::Fastest);
    let mut packet = Vec::new();
    packet.extend_from_slice(&(8 + compressed.len() as u32).to_le_bytes());
    packet.extend_from_slice(&(num_scans as u32).to_le_bytes());
    packet.extend_from_slice(&compressed);
    packet
}

pub const MZ_LOWER: f64 = 100.0;
pub const MZ_UPPER: f64 = 1600.0;
pub const SAMPLES: u32 = 1000;

/// `(sqrt(100) + (sqrt(1600) - sqrt(100)) / 1000 * tof)^2 = (10 + 0.03 tof)^2`.
pub fn mz_of(tof: u32) -> f64 {
    let slope = (MZ_UPPER.sqrt() - MZ_LOWER.sqrt()) / f64::from(SAMPLES);
    let value = MZ_LOWER.sqrt() + slope * f64::from(tof);
    value * value
}

/// Linear 1/K0 over `[0.5, 1.5]` with `MAX(NumScans) = scan_max`: exact in
/// binary for a power-of-two `scan_max`.
pub fn im_of(scan: u32, scan_max: u32) -> f64 {
    1.5 + (0.5 - 1.5) / f64::from(scan_max) * f64::from(scan)
}

/// A temporary `.d` directory, removed on drop.
pub struct Dataset {
    pub dir: PathBuf,
}

impl Drop for Dataset {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// Write `analysis.tdf` and `analysis.tdf_bin` for `frames`, then run `extra`
/// SQL on the database. `name` must be unique per test binary.
pub fn dataset(name: &str, frames: &[Frame], extra: &str) -> Dataset {
    let dir = std::env::temp_dir().join(format!("openms-tdf-{}-{name}.d", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let mut bin = Vec::new();
    let db = Connection::open(dir.join("analysis.tdf")).unwrap();
    db.execute_batch(
        "CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
         INSERT INTO GlobalMetadata VALUES
           ('TimsCompressionType', '2'), ('MzAcqRangeLower', '100.000000'),
           ('MzAcqRangeUpper', '1600.000000'), ('DigitizerNumSamples', '1000'),
           ('AcquisitionSoftware', 'timsControl'), ('OneOverK0AcqRangeLower', '0.5'),
           ('OneOverK0AcqRangeUpper', '1.5');
         CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, NumScans INTEGER,
           NumPeaks INTEGER, MsMsType INTEGER, AccumulationTime REAL, TimsId INTEGER,
           TimsCalibration INTEGER);",
    )
    .unwrap();
    for frame in frames {
        let packet = encode(frame);
        let peaks: usize = frame.scans.iter().map(Vec::len).sum();
        db.execute(
            "INSERT INTO Frames VALUES (?, ?, ?, ?, ?, ?, ?, 1)",
            rusqlite::params![
                frame.id,
                frame.time,
                frame.scans.len() as u32,
                peaks as u32,
                frame.msms_type,
                frame.accumulation_time,
                bin.len() as i64
            ],
        )
        .unwrap();
        bin.extend_from_slice(&packet);
    }
    db.execute_batch(extra).unwrap();
    std::fs::write(dir.join("analysis.tdf_bin"), bin).unwrap();
    Dataset { dir }
}

/// An MS1 frame (Id 1, 8 scans) and an MS2 frame (Id 2, 8 scans).
pub fn two_frames() -> Vec<Frame> {
    let mut ms1 = vec![Vec::new(); 8];
    ms1[0] = vec![(500, 7), (900, 3)];
    ms1[3] = vec![(100, 11)];
    ms1[7] = vec![(250, 5), (251, 6)];
    let mut ms2 = vec![Vec::new(); 8];
    ms2[2] = vec![(400, 20)];
    ms2[5] = vec![(300, 40)];
    vec![
        Frame {
            id: 1,
            time: 10.0,
            msms_type: 0,
            accumulation_time: 100.0,
            scans: ms1,
        },
        Frame {
            id: 2,
            time: 10.5,
            msms_type: 9,
            accumulation_time: 50.0,
            scans: ms2,
        },
    ]
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! TIMS scan-index to inverse-reduced-ion-mobility (1/K0) conversion: the
//! source `RationalScan2ImConverter` (`FORMAT/RationalScan2ImConverter.h`),
//! the opentims interface it implements (`Scan2InvIonMobilityConverter`) and
//! opentims' linear fallback (`OpenSourceScan2ImConverter`).
//!
//! A timsTOF frame records ions by TIMS scan index; the Bruker TDF reader
//! converts the index to 1/K0 (V·s/cm², "VSSC") with one of three
//! calibrations, tried in this order by the source's `AUTO` strategy: the
//! proprietary Bruker SDK, the per-frame rational model of the
//! `TimsCalibration` table ([`RationalScan2ImConverter`](crate::format::rational_scan2im_converter::RationalScan2ImConverter)), and a linear model
//! from `GlobalMetadata` ([`LinearScan2ImConverter`](crate::format::rational_scan2im_converter::LinearScan2ImConverter)). The Bruker SDK is a
//! closed-source shared library and is not ported; see
//! `docs/BRUKER_TIMS_SUPPORT.md`.
//!
//! The source header documents `RationalScan2ImConverter` as the first
//! open-source reader of the per-frame `TimsCalibration` coefficients; other
//! open-source readers (OpenTIMS, AlphaTIMS, timsrust) use the
//! frame-independent linear approximation.

use crate::{Error, Result};
use std::collections::BTreeMap;

/// Conversion between TIMS scan indices and 1/K0 for one frame at a time.
///
/// The opentims base class `Scan2InvIonMobilityConverter`. The source methods
/// write into caller-provided arrays of `size` elements; these return a new
/// vector of the input's length instead. Every conversion is total: no input
/// is refused, as in the source.
pub trait Scan2InvIonMobilityConverter {
    /// 1/K0 of each (possibly fractional) scan index of frame `frame_id`.
    ///
    /// The source `convert(frame_id, inv_ion_mobilities, const double* scans,
    /// size)`.
    fn convert(&self, frame_id: u32, scans: &[f64]) -> Vec<f64>;

    /// 1/K0 of each integer scan index of frame `frame_id`.
    ///
    /// The source `convert(frame_id, inv_ion_mobilities, const uint32_t*
    /// scans, size)`; each index is converted to `f64` first, so this equals
    /// [`Self::convert`] on the widened values.
    fn convert_indices(&self, frame_id: u32, scans: &[u32]) -> Vec<f64> {
        let widened: Vec<f64> = scans.iter().map(|&scan| f64::from(scan)).collect();
        self.convert(frame_id, &widened)
    }

    /// The nearest integer scan index of each 1/K0 value of frame `frame_id`.
    ///
    /// The source `inverse_convert(frame_id, scans, inv_ion_mobilities,
    /// size)`.
    fn inverse_convert(&self, frame_id: u32, inv_ion_mobilities: &[f64]) -> Vec<u32>;

    /// A human-readable name of the converter.
    ///
    /// The opentims base implementation returns
    /// `Scan2InvIonMobilityConverter default`; this default does too.
    fn description(&self) -> String {
        "Scan2InvIonMobilityConverter default".into()
    }
}

/// The source's rounding of a fractional scan index: `val > 0 ?
/// (uint32_t)(val + 0.5) : 0`.
///
/// NaN and non-positive values give 0, as in the source. A value of `2^32 -
/// 0.5` or more saturates at `u32::MAX`; the source's cast is undefined
/// behaviour there.
fn round_scan(value: f64) -> u32 {
    if value > 0.0 { (value + 0.5) as u32 } else { 0 }
}

/// Calibration coefficients from one row of the `TimsCalibration` table.
///
/// The source struct `RationalScan2ImConverter::Coefficients`, columns `C0` to
/// `C9`. `c5`, `c8` and `c9` are stored but unused by ModelType 2; the source
/// notes that no public documentation assigns them a meaning and keeps them
/// for forward compatibility.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Coefficients {
    /// `C0`, a scan-index offset.
    pub c0: f64,
    /// `C1`, the scan-index span of the voltage ramp.
    pub c1: f64,
    /// `C2`, the TIMS elution voltage at the start of the ramp.
    pub c2: f64,
    /// `C3`, the TIMS elution voltage at the end of the ramp.
    pub c3: f64,
    /// `C4`, a second scan-index offset.
    pub c4: f64,
    /// `C5`, unused by ModelType 2.
    pub c5: f64,
    /// `C6`, the constant term of the voltage-to-mobility model.
    pub c6: f64,
    /// `C7`, the `1/V` term of the voltage-to-mobility model.
    pub c7: f64,
    /// `C8`, unused by ModelType 2.
    pub c8: f64,
    /// `C9`, unused by ModelType 2.
    pub c9: f64,
}

/// Per-frame rational-function scan-to-1/K0 converter from the Bruker
/// `TimsCalibration` table.
///
/// The source class `RationalScan2ImConverter` (a
/// [`Scan2InvIonMobilityConverter`]). Physical model, ModelType 2:
///
/// ```text
/// V    = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)    scan index -> TIMS elution voltage
/// 1/K0 = 1 / (c6 + c7 / V)                            voltage -> inverse reduced mobility
/// ```
///
/// Limitations the source states against the Bruker SDK: no pressure
/// compensation (ambient gas pressure drift), no awareness of a recalibrated
/// state (the raw acquisition calibration is used), and only ModelType 2
/// (other types make [`try_create_rational_converter`] fall back to linear).
///
/// Immutable after construction, so concurrent conversions are safe, as the
/// source states for its `convert` and `inverse_convert`.
#[derive(Clone, Debug, PartialEq)]
pub struct RationalScan2ImConverter {
    calibrations: BTreeMap<u32, Coefficients>,
    frame_to_cal: Vec<u32>,
}

impl RationalScan2ImConverter {
    /// A converter from calibration segments and the frame-to-segment map.
    ///
    /// The source constructor `RationalScan2ImConverter(calibrations,
    /// frame_to_cal)`. `calibrations` maps a `TimsCalibration` Id to its
    /// coefficients; `frame_to_cal[frame_id]` is the calibration Id of frame
    /// `frame_id`, 1-based, so index 0 is unused.
    ///
    /// The source takes an `unordered_map`; a `BTreeMap` fixes the fallback
    /// segment (see [`Self::calibration`]).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `calibrations` is empty. The source accepts
    /// it and then dereferences `calibrations_.begin()` of an empty map on the
    /// first conversion, which is undefined behaviour.
    pub fn new(calibrations: BTreeMap<u32, Coefficients>, frame_to_cal: Vec<u32>) -> Result<Self> {
        if calibrations.is_empty() {
            return Err(Error::InvalidValue(
                "RationalScan2ImConverter needs at least one calibration segment".into(),
            ));
        }
        Ok(Self {
            calibrations,
            frame_to_cal,
        })
    }

    /// The calibration segments by `TimsCalibration` Id.
    pub fn calibrations(&self) -> &BTreeMap<u32, Coefficients> {
        &self.calibrations
    }

    /// The calibration Id of each frame, indexed by 1-based frame Id.
    pub fn frame_to_cal(&self) -> &[u32] {
        &self.frame_to_cal
    }

    /// The coefficients used for frame `frame_id`.
    ///
    /// The source private `getCalibration`: the segment that
    /// `frame_to_cal[frame_id]` names when `frame_id` is in `1..len` and the
    /// segment exists; otherwise a warning is logged (unless `frame_id` is 0)
    /// and the *first* segment is used, which the source says should not
    /// happen in valid data.
    ///
    /// The source's first segment is `unordered_map::begin()`, whose element
    /// depends on the hash table's layout when there is more than one segment;
    /// here it is the segment with the smallest Id.
    pub fn calibration(&self, frame_id: u32) -> &Coefficients {
        let found = usize::try_from(frame_id)
            .ok()
            .filter(|&index| index > 0)
            .and_then(|index| self.frame_to_cal.get(index))
            .and_then(|cal_id| self.calibrations.get(cal_id));
        if let Some(coefficients) = found {
            return coefficients;
        }
        if frame_id != 0 {
            let _ = crate::openms_log_warn!(
                "RationalScan2ImConverter: frame_id {frame_id} out of range, using first calibration"
            );
        }
        self.calibrations
            .values()
            .next()
            .expect("constructor refuses an empty calibration map")
    }

    /// The forward model for one scan index.
    ///
    /// The source private static `applyFormula`:
    /// `V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)`, then
    /// `1/K0 = 1 / (c6 + c7 / V)`. A `V` or denominator of exactly 0 is
    /// replaced by `1e-10`, the source's guard against division by zero
    /// (physically unreachable).
    pub fn apply_formula(c: &Coefficients, scan: f64) -> f64 {
        let mut voltage = c.c2 + ((c.c3 - c.c2) / c.c1) * (scan - c.c4 - c.c0);
        if voltage == 0.0 {
            voltage = 1e-10;
        }
        let mut denominator = c.c6 + c.c7 / voltage;
        if denominator == 0.0 {
            denominator = 1e-10;
        }
        1.0 / denominator
    }

    /// The inverse model for one 1/K0 value.
    ///
    /// The source private static `invertFormula`:
    /// `scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1/inv_k0 - c6) - c2)`, with
    /// the two denominators guarded like [`Self::apply_formula`].
    pub fn invert_formula(c: &Coefficients, inv_k0: f64) -> f64 {
        let reciprocal = 1.0 / inv_k0;
        let mut inner = reciprocal - c.c6;
        if inner == 0.0 {
            inner = 1e-10;
        }
        let mut slope = c.c3 - c.c2;
        if slope == 0.0 {
            slope = 1e-10;
        }
        c.c0 + c.c4 + (c.c1 / slope) * (c.c7 / inner - c.c2)
    }
}

impl Scan2InvIonMobilityConverter for RationalScan2ImConverter {
    fn convert(&self, frame_id: u32, scans: &[f64]) -> Vec<f64> {
        let c = self.calibration(frame_id);
        scans
            .iter()
            .map(|&scan| Self::apply_formula(c, scan))
            .collect()
    }

    fn inverse_convert(&self, frame_id: u32, inv_ion_mobilities: &[f64]) -> Vec<u32> {
        let c = self.calibration(frame_id);
        inv_ion_mobilities
            .iter()
            .map(|&value| round_scan(Self::invert_formula(c, value)))
            .collect()
    }

    /// `RationalScan2ImConverter (N calibration segment[s])`, the source
    /// format.
    fn description(&self) -> String {
        let count = self.calibrations.len();
        format!(
            "RationalScan2ImConverter ({count} calibration segment{})",
            if count != 1 { "s" } else { "" }
        )
    }
}

/// Frame-independent linear scan-to-1/K0 converter, opentims'
/// `OpenSourceScan2ImConverter`.
///
/// `1/K0 = intercept + slope * scan`, with `intercept = im_max` and
/// `slope = (im_min - im_max) / scan_max_index`: scan 0 maps to the upper end
/// of the acquisition range and `scan_max_index` to the lower end. The TDF
/// reader takes `im_min`/`im_max` from the `GlobalMetadata` keys
/// `OneOverK0AcqRangeLower`/`OneOverK0AcqRangeUpper` and `scan_max_index` from
/// the largest `Frames.NumScans`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearScan2ImConverter {
    intercept: f64,
    slope: f64,
}

impl LinearScan2ImConverter {
    /// The opentims constructor `OpenSourceScan2ImConverter(im_min, im_max,
    /// scan_max_index)`. A `scan_max_index` of 0 gives an infinite or NaN
    /// slope, as in the source.
    pub fn new(im_min: f64, im_max: f64, scan_max_index: u32) -> Self {
        Self {
            intercept: im_max,
            slope: (im_min - im_max) / f64::from(scan_max_index),
        }
    }

    /// 1/K0 at scan 0.
    pub fn intercept(&self) -> f64 {
        self.intercept
    }

    /// Change of 1/K0 per scan.
    pub fn slope(&self) -> f64 {
        self.slope
    }
}

impl Scan2InvIonMobilityConverter for LinearScan2ImConverter {
    fn convert(&self, _frame_id: u32, scans: &[f64]) -> Vec<f64> {
        scans
            .iter()
            .map(|&scan| self.intercept + self.slope * scan)
            .collect()
    }

    fn inverse_convert(&self, _frame_id: u32, inv_ion_mobilities: &[f64]) -> Vec<u32> {
        inv_ion_mobilities
            .iter()
            .map(|&value| round_scan((value - self.intercept) / self.slope))
            .collect()
    }

    /// `OpenSourceScan2ImConverter (linear)`, the opentims text.
    fn description(&self) -> String {
        "OpenSourceScan2ImConverter (linear)".into()
    }
}

/// Largest `Frames.Id` [`try_create_rational_converter`] accepts.
///
/// Native ceiling: the source sizes its frame-to-calibration vector by the
/// largest frame Id, so a corrupt Id near 2^32 allocates 16 GiB (and `Id + 1`
/// wraps to 0 for `u32::MAX`). A timsTOF run has about 10^5 frames; 2^24
/// frames (a 64 MiB map) is two orders of magnitude above that.
pub const MAX_FRAME_ID: u32 = 1 << 24;

/// Try to create a [`RationalScan2ImConverter`] from the `TimsCalibration`
/// table of the `analysis.tdf` in the `.d` directory `tims_dir_path`.
///
/// The source free function `tryCreateRationalConverter`. Returns `None`, so
/// the caller falls back to the linear converter, when the database cannot be
/// opened or has no `TimsCalibration` table, when a row has a ModelType other
/// than 2 (warning), when the table is empty, when `Frames` has no
/// `TimsCalibration` column (an old TDF version), when a frame's
/// `TimsCalibration` is NULL or names a missing Id (warning), and when there
/// are no frames. Otherwise the converter maps every frame Id to its
/// calibration Id; frames missing from `Frames` map to Id 0, as the source's
/// zero-filled vector does.
///
/// Column values are read with SQLite's own conversions, as SQLiteCpp's
/// `getInt`/`getDouble` do: NULL reads as 0, an integer Id keeps its low 32
/// bits, text reads as its leading number.
///
/// Native difference: a `Frames.Id` above [`MAX_FRAME_ID`] also returns
/// `None` with a warning, instead of the source's allocation of a vector that
/// large.
#[cfg(feature = "sqlite")]
pub fn try_create_rational_converter(
    tims_dir_path: impl AsRef<std::path::Path>,
) -> Option<RationalScan2ImConverter> {
    use rusqlite::types::ValueRef;
    use rusqlite::{Connection, OpenFlags};

    fn int(value: ValueRef<'_>) -> u32 {
        match value {
            ValueRef::Integer(i) => i as u32,
            ValueRef::Real(r) => (r as i64) as u32,
            ValueRef::Text(t) => leading_number(t) as i64 as u32,
            ValueRef::Null | ValueRef::Blob(_) => 0,
        }
    }
    fn double(value: ValueRef<'_>) -> f64 {
        match value {
            ValueRef::Integer(i) => i as f64,
            ValueRef::Real(r) => r,
            ValueRef::Text(t) => leading_number(t),
            ValueRef::Null | ValueRef::Blob(_) => 0.0,
        }
    }

    let tdf_path = tims_dir_path.as_ref().join("analysis.tdf");
    let db = match Connection::open_with_flags(
        &tdf_path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    ) {
        Ok(db) => db,
        Err(_) => {
            let _ = crate::openms_log_debug_nofile!(
                "TimsCalibration table not found in {}, falling back to linear converter",
                tdf_path.display()
            );
            return None;
        }
    };

    let mut calibrations = BTreeMap::new();
    let mut read_calibrations = || -> rusqlite::Result<Option<()>> {
        let mut statement = db.prepare(
            "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration",
        )?;
        let mut rows = statement.query([])?;
        while let Some(row) = rows.next()? {
            let id = int(row.get_ref(0)?);
            let model_type = int(row.get_ref(1)?) as i32;
            if model_type != 2 {
                let _ = crate::openms_log_warn!(
                    "TimsCalibration ModelType={model_type} unsupported (only ModelType=2 implemented), falling back to linear converter"
                );
                return Ok(None);
            }
            let c = |i: usize| row.get_ref(i).map(double);
            calibrations.insert(
                id,
                Coefficients {
                    c0: c(2)?,
                    c1: c(3)?,
                    c2: c(4)?,
                    c3: c(5)?,
                    c4: c(6)?,
                    c5: c(7)?,
                    c6: c(8)?,
                    c7: c(9)?,
                    c8: c(10)?,
                    c9: c(11)?,
                },
            );
        }
        Ok(Some(()))
    };
    match read_calibrations() {
        Ok(Some(())) => {}
        Ok(None) => return None,
        Err(_) => {
            let _ = crate::openms_log_debug_nofile!(
                "TimsCalibration table not found in {}, falling back to linear converter",
                tdf_path.display()
            );
            return None;
        }
    }
    if calibrations.is_empty() {
        let _ = crate::openms_log_debug_nofile!(
            "TimsCalibration table empty, falling back to linear converter"
        );
        return None;
    }

    let mut frame_to_cal: Vec<u32> = Vec::new();
    let mut read_frames = || -> rusqlite::Result<Option<()>> {
        let mut statement = db.prepare("SELECT Id, TimsCalibration FROM Frames ORDER BY Id")?;
        let mut rows = statement.query([])?;
        while let Some(row) = rows.next()? {
            let frame_id = int(row.get_ref(0)?);
            let cal = row.get_ref(1)?;
            if matches!(cal, ValueRef::Null) {
                let _ = crate::openms_log_warn!(
                    "Frame {frame_id} has NULL TimsCalibration FK, falling back to linear converter"
                );
                return Ok(None);
            }
            if frame_id > MAX_FRAME_ID {
                let _ = crate::openms_log_warn!(
                    "Frame {frame_id} exceeds the frame Id ceiling {MAX_FRAME_ID}, falling back to linear converter"
                );
                return Ok(None);
            }
            let cal_id = int(cal);
            let index = frame_id as usize;
            if frame_to_cal.len() <= index {
                frame_to_cal.resize(index + 1, 0);
            }
            frame_to_cal[index] = cal_id;
            if !calibrations.contains_key(&cal_id) {
                let _ = crate::openms_log_warn!(
                    "Frame {frame_id} references unknown TimsCalibration Id={cal_id}, falling back to linear converter"
                );
                return Ok(None);
            }
        }
        Ok(Some(()))
    };
    match read_frames() {
        Ok(Some(())) => {}
        Ok(None) => return None,
        Err(_) => {
            let _ = crate::openms_log_debug_nofile!(
                "Frames table has no TimsCalibration column, falling back to linear converter"
            );
            return None;
        }
    }
    if frame_to_cal.len() <= 1 {
        let _ =
            crate::openms_log_debug_nofile!("No frames found, falling back to linear converter");
        return None;
    }
    let count = calibrations.len();
    let _ = crate::openms_log_info!(
        "TIMS calibration: rational (TimsCalibration table, {count} segment{} for {} frames)",
        if count != 1 { "s" } else { "" },
        frame_to_cal.len() - 1
    );
    RationalScan2ImConverter::new(calibrations, frame_to_cal).ok()
}

/// SQLite's text-to-number conversion: the longest leading decimal number,
/// after leading whitespace; 0 when there is none.
#[cfg(feature = "sqlite")]
fn leading_number(text: &[u8]) -> f64 {
    let text = String::from_utf8_lossy(text);
    let text = text.trim_start();
    let mut best = 0.0;
    for (index, ch) in text.char_indices() {
        if !(ch.is_ascii_digit() || matches!(ch, '+' | '-' | '.' | 'e' | 'E')) {
            break;
        }
        if let Ok(value) = text[..index + ch.len_utf8()].parse::<f64>() {
            best = value;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_calibrations_are_refused() {
        assert!(RationalScan2ImConverter::new(BTreeMap::new(), vec![0, 1]).is_err());
    }

    #[test]
    fn fallback_uses_the_smallest_id() {
        let a = Coefficients {
            c1: 1.0,
            c6: 1.0,
            ..Coefficients::default()
        };
        let b = Coefficients {
            c1: 1.0,
            c6: 2.0,
            ..Coefficients::default()
        };
        let converter =
            RationalScan2ImConverter::new(BTreeMap::from([(7, b), (3, a)]), vec![0, 7]).unwrap();
        assert_eq!(converter.calibration(1), &b);
        assert_eq!(converter.calibration(0), &a);
        assert_eq!(converter.calibration(99), &a);
    }

    #[test]
    fn rounding_matches_the_source_cast() {
        assert_eq!(round_scan(-1.0), 0);
        assert_eq!(round_scan(f64::NAN), 0);
        assert_eq!(round_scan(2.49), 2);
        assert_eq!(round_scan(2.5), 3);
        assert_eq!(round_scan(1e20), u32::MAX);
    }

    #[test]
    fn linear_converter_maps_the_range_ends() {
        // Expected values: opentims `OpenSourceScan2ImConverter`
        // (scan2inv_ion_mobility_converter.cpp:170-199), `intercept_ + slope_ *
        // scan` with `slope_ = (im_min - im_max) / scan_max_index`, evaluated in
        // IEEE-754 doubles outside Rust (Python 3). The far end is not 0.6: the
        // slope times 1000 rounds to -1.0, and the double 1.6 minus 1.0 is the
        // double just above 0.6.
        let linear = LinearScan2ImConverter::new(0.6, 1.6, 1000);
        assert_eq!(linear.convert(1, &[0.0]), vec![1.6]);
        assert_eq!(
            linear.convert_indices(1, &[1000]),
            vec![0.600_000_000_000_000_1]
        );
        assert_eq!(linear.inverse_convert(1, &[1.6, 0.6]), vec![0, 1000]);
        // A dyadic range, where every operation is exact.
        let dyadic = LinearScan2ImConverter::new(0.5, 1.5, 1024);
        assert_eq!(
            dyadic.convert_indices(1, &[0, 512, 1024]),
            vec![1.5, 1.0, 0.5]
        );
        assert_eq!(
            dyadic.inverse_convert(1, &[1.5, 1.0, 0.5, 2.0]),
            vec![0, 512, 1024, 0]
        );
        assert_eq!(linear.description(), "OpenSourceScan2ImConverter (linear)");
    }

    #[cfg(feature = "sqlite")]
    #[test]
    fn leading_number_follows_sqlite() {
        assert_eq!(leading_number(b"  12.5abc"), 12.5);
        assert_eq!(leading_number(b"abc"), 0.0);
        assert_eq!(leading_number(b"-3"), -3.0);
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for Kroenik (Hardklör sibling) feature files.
//!
//! Ports `FORMAT/KroenikFile.h`/`.cpp`. A Kroenik file is tab-separated: the
//! first line is a header naming the columns
//!
//! ```text
//! File  First Scan  Last Scan  Num of Scans  Charge  Monoisotopic Mass
//! Base Isotope Peak  Best Intensity  Summed Intensity  First RTime
//! Last RTime  Best RTime  Best Correlation  Modifications
//! ```
//!
//! and every later line is one feature. See `docs/SEARCH_ENGINE_FILES_SUPPORT.md`
//! for the API mapping, preserved conventions and native differences.

use crate::concept::constants::PROTON_MASS_U;
use crate::data_structures::string_utils::{split, to_double, to_int32};
use crate::format::text::{ReadOptions, TextFile};
use crate::kernel::features::{Feature, FeatureMap};
use crate::kernel::geometry::{ConvexHull2D, Point2D};
use crate::metadata::MetaValue;
use crate::{Error, MSSpectrum, Result};
use std::io::BufRead;
use std::path::Path;

/// The number of tab-separated columns every Kroenik data line must have.
pub const COLUMN_COUNT: usize = 14;

/// The m/z width, in Th for charge one, the source gives every convex hull,
/// since a Kroenik file carries no m/z extent (divided by the charge).
pub const HULL_MZ_WIDTH: f64 = 3.0;

/// File adapter for Kroenik (Hardklör sibling) files.
///
/// Every property in the file becomes a [`Feature`] property, while
/// "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
/// meta values named `FirstScan`, `LastScan`, `NumOfScans` and
/// `AveragineModifications`; the monoisotopic mass is kept as `Mass`.
///
/// The m/z width of each feature's convex hull is set to 3 Th (divided by the
/// charge) for lack of a value in the file.
///
/// Kroenik files are tab (`\t`) separated.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct KroenikFile;

impl KroenikFile {
    /// A Kroenik file adapter; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load a Kroenik file into a new feature map.
    ///
    /// The first line is skipped as the header whatever it contains; an empty
    /// file yields an empty map. Each data line becomes one feature: charge
    /// from column 5, m/z as `Monoisotopic Mass / charge + proton mass`, RT
    /// from `Best RTime`, overall quality from `Best Correlation`, intensity
    /// from `Summed Intensity`, and a two-scan convex hull spanning
    /// `First RTime`..`Last RTime` by `m/z`..`m/z + 3/charge`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened (source
    /// `Exception::FileNotFound`); [`Error::Parse`] when a data line does not
    /// have exactly 14 tab-separated entries (source `Exception::ParseError`,
    /// with the same message) or when a numeric field does not convert.
    /// See [`KroenikFile::load_reader`] for the numeric rules.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<FeatureMap> {
        let input = TextFile::from_path(path, &ReadOptions::default())?;
        self.parse(&input)
    }

    /// Load Kroenik content from a reader; see [`KroenikFile::load`].
    ///
    /// # Errors
    ///
    /// As [`KroenikFile::load`]. A numeric field the source's
    /// `StringUtils::toDouble`/`toInt32` would refuse is an [`Error::Parse`]
    /// naming the line; the source lets that `Exception::ConversionError`
    /// escape unconverted. Non-finite values, and a charge of zero (which
    /// divides the mass by zero in the source and stores an infinite m/z), are
    /// refused, because a [`Feature`] with a non-finite coordinate fails every
    /// later check. On any error nothing is returned; the source has already
    /// replaced the caller's map and appended the lines before the failure.
    pub fn load_reader(&self, reader: impl BufRead) -> Result<FeatureMap> {
        let input = TextFile::from_reader(reader, &ReadOptions::default())?;
        self.parse(&input)
    }

    fn parse(&self, input: &TextFile) -> Result<FeatureMap> {
        let mut map = FeatureMap::new();
        for (index, line) in input.lines().iter().enumerate().skip(1) {
            let number = index + 1;
            let parts = split(line, '\t');
            if parts.len() != COLUMN_COUNT {
                return Err(Error::Parse {
                    line: number,
                    message: format!(
                        "Failed parsing in line {number}: missing 14 tab-separated entries (got {})\nLine was: '{line}'",
                        parts.len()
                    ),
                });
            }
            map.features.push(feature(&parts, number)?);
        }
        Ok(map)
    }

    /// Store a map as a Kroenik file: not implemented, as in the source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`], the source's `Exception::NotImplemented`.
    /// The source also prints the file name and spectrum size to standard
    /// error; the port reports them in the error message instead.
    pub fn store(&self, path: impl AsRef<Path>, spectrum: &MSSpectrum) -> Result<()> {
        Err(Error::Unsupported(format!(
            "Store() for KroenikFile not implemented. Filename was: {}, spec of size {}",
            path.as_ref().display(),
            spectrum.peaks.len()
        )))
    }
}

fn number(text: &str, line: usize) -> Result<f64> {
    to_double(text)
        .map_err(|message| Error::Parse { line, message })
        .and_then(|value| finite(value, line))
}

fn integer(text: &str, line: usize) -> Result<i32> {
    to_int32(text).map_err(|message| Error::Parse { line, message })
}

fn finite(value: f64, line: usize) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Parse {
            line,
            message: format!("non-finite value {value} in a Kroenik feature"),
        })
    }
}

fn feature(parts: &[&str], line: usize) -> Result<Feature> {
    let charge = integer(parts[4], line)?;
    let mass = number(parts[5], line)?;
    let mz = finite(mass / f64::from(charge) + PROTON_MASS_U, line)?;
    let mut feature = Feature::new(number(parts[11], line)?, mz, 0.0);
    feature.charge = charge;
    feature.quality = number(parts[12], line)? as f32;
    feature.intensity = number(parts[8], line)? as f32;
    finite(f64::from(feature.quality), line)?;
    finite(f64::from(feature.intensity), line)?;
    let first_rt = number(parts[9], line)?;
    let last_rt = number(parts[10], line)?;
    let top = finite(mz + HULL_MZ_WIDTH / f64::from(charge), line)?;
    let hull = ConvexHull2D::from_points(&[
        Point2D::new(first_rt, mz),
        Point2D::new(first_rt, top),
        Point2D::new(last_rt, top),
        Point2D::new(last_rt, mz),
        Point2D::new(first_rt, mz),
    ])?;
    feature.convex_hulls = vec![hull];
    let meta = &mut feature.metadata;
    meta.insert("Mass".into(), MetaValue::try_from(mass)?);
    meta.insert(
        "FirstScan".into(),
        MetaValue::try_from(number(parts[1], line)?)?,
    );
    meta.insert("LastScan".into(), integer(parts[2], line)?.into());
    meta.insert(
        "NumOfScans".into(),
        MetaValue::try_from(number(parts[3], line)?)?,
    );
    meta.insert("AveragineModifications".into(), parts[13].into());
    Ok(feature)
}

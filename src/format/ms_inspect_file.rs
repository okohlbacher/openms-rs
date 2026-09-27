// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for msInspect feature files.
//!
//! Ports `FORMAT/MsInspectFile.h`/`.cpp`. Lines beginning with `#` are
//! comments and are ignored. The first non-comment line is the header naming
//! the columns
//!
//! ```text
//! scan  time  mz  accurateMZ  mass  intensity  charge  chargeStates  kl
//! background  median  peaks  scanFirst  scanLast  scanCount  totalIntensity
//! sumSquaresDist  description
//! ```
//!
//! and every later line is one feature. See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::format::source_text::{split, to_double, to_int32};
use crate::format::text::{ReadOptions, TextFile};
use crate::kernel::features::{Feature, FeatureMap};
use crate::metadata::MetaValue;
use crate::{Error, MSSpectrum, Result};
use std::io::BufRead;
use std::path::Path;

/// The minimum number of tab-separated columns of an msInspect data line.
pub const MIN_COLUMN_COUNT: usize = 18;

/// File adapter for msInspect files.
///
/// Columns `time`, `mz`, `intensity`, `charge` and `kl` become the feature's
/// RT, m/z, intensity, charge and overall quality. `accurateMZ` is stored as a
/// *string* meta value, exactly as the file spells it, and `description` as a
/// string; `mass`, `background`, `median`, `totalIntensity` and
/// `sumSquaresDist` become floating meta values and `chargeStates`, `peaks`,
/// `scanFirst`, `scanLast` and `scanCount` integer meta values, each under its
/// column name. Column `scan` is not read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MsInspectFile;

impl MsInspectFile {
    /// An msInspect file adapter; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load an msInspect file into a new feature map.
    ///
    /// Empty lines and lines starting with `#` are skipped wherever they
    /// occur; the first remaining line is skipped as the header. Columns past
    /// the eighteenth are ignored.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened (source
    /// `Exception::FileNotFound`); [`Error::Parse`] when a data line has fewer
    /// than 18 columns, or when a numeric column does not convert, with the
    /// source's `ParseError` messages naming the line and the one-based column.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<FeatureMap> {
        let input = TextFile::from_path(path, &ReadOptions::default())?;
        parse(&input)
    }

    /// Load msInspect content from a reader; see [`MsInspectFile::load`].
    ///
    /// # Errors
    ///
    /// As [`MsInspectFile::load`]. Non-finite values in the RT, m/z,
    /// intensity, quality or floating meta-value columns are refused as a
    /// conversion failure of that column; the source stores them. On any error
    /// nothing is returned, whereas the source has already replaced the
    /// caller's map and appended the lines before the failure.
    pub fn load_reader(&self, reader: impl BufRead) -> Result<FeatureMap> {
        let input = TextFile::from_reader(reader, &ReadOptions::default())?;
        parse(&input)
    }

    /// Store a map as an msInspect file: not implemented, as in the source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`], the source's `Exception::NotImplemented`.
    /// The source also prints the file name and spectrum size to standard
    /// error; the port reports them in the error message instead.
    pub fn store(&self, path: impl AsRef<Path>, spectrum: &MSSpectrum) -> Result<()> {
        Err(Error::Unsupported(format!(
            "Store() for MsInspectFile not implemented. Filename was: {}, spec of size {}",
            path.as_ref().display(),
            spectrum.peaks.len()
        )))
    }
}

fn parse(input: &TextFile) -> Result<FeatureMap> {
    let mut map = FeatureMap::new();
    let mut first_line = true;
    for (index, line) in input.lines().iter().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if first_line {
            first_line = false;
            continue;
        }
        let number = index + 1;
        let parts = split(line, '\t');
        if parts.len() < MIN_COLUMN_COUNT {
            return Err(Error::Parse {
                line: number,
                message: format!(
                    "Failed to convert line {number}. Not enough columns (expected 18 or more, got {})",
                    parts.len()
                ),
            });
        }
        map.features.push(feature(&parts, number)?);
    }
    Ok(map)
}

fn feature(parts: &[&str], line: usize) -> Result<Feature> {
    let column_error = |column: usize| Error::Parse {
        line,
        message: format!(
            "Failed to convert value in column {} into a number (line '{line})",
            column + 1
        ),
    };
    let float = |column: usize| {
        to_double(parts[column])
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(|| column_error(column))
    };
    let int = |column: usize| to_int32(parts[column]).map_err(|_| column_error(column));
    let float_meta = |column: usize| MetaValue::try_from(float(column)?);

    let mut feature = Feature::new(float(1)?, float(2)?, 0.0);
    feature.intensity = float(5)? as f32;
    if !feature.intensity.is_finite() {
        return Err(column_error(5));
    }
    feature.charge = int(6)?;
    feature.quality = float(8)? as f32;
    if !feature.quality.is_finite() {
        return Err(column_error(8));
    }
    let meta = &mut feature.metadata;
    meta.insert("accurateMZ".into(), parts[3].into());
    meta.insert("mass".into(), float_meta(4)?);
    meta.insert("chargeStates".into(), int(7)?.into());
    meta.insert("background".into(), float_meta(9)?);
    meta.insert("median".into(), float_meta(10)?);
    meta.insert("peaks".into(), int(11)?.into());
    meta.insert("scanFirst".into(), int(12)?.into());
    meta.insert("scanLast".into(), int(13)?.into());
    meta.insert("scanCount".into(), int(14)?.into());
    meta.insert("totalIntensity".into(), float_meta(15)?);
    meta.insert("sumSquaresDist".into(), float_meta(16)?);
    meta.insert("description".into(), parts[17].into());
    Ok(feature)
}

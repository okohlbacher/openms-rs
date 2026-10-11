// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File adapter for SpecArray (`.pepList`) feature files.
//!
//! Ports `FORMAT/SpecArrayFile.h`/`.cpp`. The first line is a header naming
//! the columns `m/z`, `rt(min)`, `snr`, `charge` and `intensity`; every later
//! line is one feature, entries separated by tabs. See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping.

use crate::concept::string_conversion::{split, to_double, to_int32};
use crate::format::text::{ReadOptions, TextFile};
use crate::kernel::features::{Feature, FeatureMap};
use crate::metadata::MetaValue;
use crate::{Error, MSSpectrum, Result};
use std::io::BufRead;
use std::path::Path;

/// The minimum number of tab-separated columns of a SpecArray data line.
pub const MIN_COLUMN_COUNT: usize = 5;

/// File adapter for SpecArray (`.pepList`) files.
///
/// Each data line becomes a feature with m/z from column 1, RT from column 2
/// converted from minutes to seconds (times 60), the signal-to-noise ratio of
/// column 3 as the floating meta value `s/n`, charge from column 4 and
/// intensity from column 5.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpecArrayFile;

impl SpecArrayFile {
    /// A SpecArray file adapter; it holds no state.
    pub fn new() -> Self {
        Self
    }

    /// Load a SpecArray file into a new feature map.
    ///
    /// The first line is skipped as the header whatever it contains; an empty
    /// file yields an empty map. Columns past the fifth are ignored.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened (source
    /// `Exception::FileNotFound`); [`Error::Parse`] when a data line has fewer
    /// than five columns or a value does not convert to a number, with the
    /// source's `ParseError` messages.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<FeatureMap> {
        let input = TextFile::from_path(path, &ReadOptions::default())?;
        parse(&input)
    }

    /// Load SpecArray content from a reader; see [`SpecArrayFile::load`].
    ///
    /// # Errors
    ///
    /// As [`SpecArrayFile::load`]. Non-finite values are refused as a
    /// conversion failure; the source stores them. On any error nothing is
    /// returned, whereas the source has already replaced the caller's map and
    /// appended the lines before the failure.
    pub fn load_reader(&self, reader: impl BufRead) -> Result<FeatureMap> {
        let input = TextFile::from_reader(reader, &ReadOptions::default())?;
        parse(&input)
    }

    /// Store a map as a SpecArray file: not implemented, as in the source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`], the source's `Exception::NotImplemented`.
    /// The source also prints the file name and spectrum size to standard
    /// error; the port reports them in the error message instead.
    pub fn store(&self, path: impl AsRef<Path>, spectrum: &MSSpectrum) -> Result<()> {
        Err(Error::Unsupported(format!(
            "Store() for SpecArrayFile not implemented. Filename was: {}, spec of size {}",
            path.as_ref().display(),
            spectrum.peaks.len()
        )))
    }
}

fn parse(input: &TextFile) -> Result<FeatureMap> {
    let mut map = FeatureMap::new();
    for (index, line) in input.lines().iter().enumerate().skip(1) {
        let number = index + 1;
        let parts = split(line, '\t');
        if parts.len() < MIN_COLUMN_COUNT {
            // The source's message lacks the spaces around the line number.
            return Err(Error::Parse {
                line: number,
                message: format!(
                    "Failed to convert line{number}not enough columns (expected 5 or more, got {})",
                    parts.len()
                ),
            });
        }
        map.features.push(feature(&parts, number)?);
    }
    Ok(map)
}

fn feature(parts: &[&str], line: usize) -> Result<Feature> {
    let failed = || Error::Parse {
        line,
        message: format!("Failed to convert value into a number (line '{line})"),
    };
    let float = |column: usize| {
        to_double(parts[column])
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(failed)
    };
    let mut feature = Feature::new(0.0, float(0)?, 0.0);
    feature.rt = float(1)? * 60.0;
    if !feature.rt.is_finite() {
        return Err(failed());
    }
    feature
        .metadata
        .insert("s/n".into(), MetaValue::try_from(float(2)?)?);
    feature.charge = to_int32(parts[3]).map_err(|_| failed())?;
    feature.intensity = float(4)? as f32;
    if !feature.intensity.is_finite() {
        return Err(failed());
    }
    Ok(feature)
}

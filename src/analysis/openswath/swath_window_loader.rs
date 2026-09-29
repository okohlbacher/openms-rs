// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reading SWATH-window definition files and annotating SWATH maps with them.
//!
//! Ports `ANALYSIS/OPENSWATH/SwathWindowLoader.h` with `SwathWindowLoader.cpp`,
//! a class of two static members that become free functions. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! # File format
//!
//! Tab- or space-delimited, two columns per row:
//!
//! ```text
//! window_lower window_upper
//! 400 425
//! 425 450
//! ```
//!
//! The first line is a header only when its first token cannot be parsed as a
//! number; otherwise it is the first data row. No particular header text is
//! required, and files without a header are accepted.
//!
//! # Differences from the source
//!
//! The source reads each data row with two `operator>>` extractions and never
//! checks them. A blank (or whitespace-only) row then leaves both values
//! unchanged, so it repeats the previous row, or reads uninitialised values
//! when it is the first row (an empty or header-only file). A row with one
//! column repeats the previous upper bound. The port refuses all of these with
//! [`Error::Parse`] instead of inventing a window. Only the first two columns
//! are read; further columns are ignored, as in the source.

use super::swath_map::SwathMap;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::io::Read;
use std::path::Path;

/// Largest SWATH-window file read, in bytes.
///
/// Native ceiling; the source reads without limit. A window file has one short
/// line per isolation window, so 16 MiB is several orders of magnitude above
/// any real acquisition scheme.
pub const MAX_SWATH_WINDOW_FILE_BYTES: u64 = 16 << 20;

/// Parse one data token as the source's `operator>>` would accept it in full.
fn parse_bound(token: &str, line: usize) -> Result<f64> {
    token.parse::<f64>().map_err(|_| Error::Parse {
        line,
        message: format!("cannot read '{token}' as a SWATH window boundary"),
    })
}

/// Whether the (trimmed, tab-substituted) first line is data, the source's
/// `StringUtils::toDouble(headerSubstrings[0])` succeeding.
fn first_line_is_data(line: &str) -> bool {
    let normalized = line.trim().replace('\t', " ");
    let first = normalized.split(' ').next().unwrap_or("");
    let first = first.strip_prefix('+').unwrap_or(first);
    !first.is_empty() && first.parse::<f64>().is_ok()
}

/// Parse the content of a SWATH-window file into lower and upper bounds.
///
/// The parsing rules of [`read_swath_windows`], on text already in memory;
/// line numbers in errors are one-based.
///
/// # Errors
///
/// As [`read_swath_windows`], without the file errors.
pub fn parse_swath_windows(content: &str) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut lines = content.lines().enumerate();
    let (mut number, mut line) = lines.next().map_or((0, ""), |(i, l)| (i, l));
    if !first_line_is_data(line) {
        // Header: the source reads the next line, which is empty at the end
        // of the file.
        match lines.next() {
            Some((i, l)) => {
                number = i;
                line = l;
            }
            None => {
                number += 1;
                line = "";
            }
        }
    }
    let mut lower_bounds = Vec::new();
    let mut upper_bounds = Vec::new();
    loop {
        let mut tokens = line.split_whitespace();
        let (lower, upper) = match (tokens.next(), tokens.next()) {
            (Some(lower), Some(upper)) => (
                parse_bound(lower, number + 1)?,
                parse_bound(upper, number + 1)?,
            ),
            _ => {
                return Err(Error::Parse {
                    line: number + 1,
                    message: "a SWATH window row needs a lower and an upper boundary".into(),
                });
            }
        };
        if lower.partial_cmp(&upper) != Some(Ordering::Less)
            || !lower.is_finite()
            || !upper.is_finite()
        {
            return Err(Error::InvalidValue(format!(
                "Swath window file contains illegal ranges: {line}"
            )));
        }
        lower_bounds.push(lower);
        upper_bounds.push(upper);
        match lines.next() {
            Some((i, l)) => {
                number = i;
                line = l;
            }
            None => break,
        }
    }
    Ok((lower_bounds, upper_bounds))
}

/// Read a SWATH-window definition file into parallel lower and upper bounds,
/// in file order.
///
/// Source `readSwathWindows`, whose two out-parameters become the returned
/// pair (of equal length). The format is described in the
/// [module documentation](self); numbers use Rust's `f64` grammar, which
/// accepts the source's forms such as `4e2`.
///
/// # Errors
///
/// - [`Error::Io`] when the file is missing or unreadable (the source's
///   `Exception::FileNotFound`), or is not UTF-8.
/// - [`Error::InvalidValue`] when a row has `lower >= upper` or a non-finite
///   bound (the source's `Exception::InvalidValue`, "illegal ranges"). The
///   source has already appended the offending row when it throws; nothing is
///   returned here.
/// - [`Error::Parse`] for a row without two numeric columns, including a
///   blank row and an empty or header-only file, which the source reads as a
///   repeat of the previous row or as uninitialised values.
/// - [`Error::InvalidRange`] when the file exceeds
///   [`MAX_SWATH_WINDOW_FILE_BYTES`].
pub fn read_swath_windows(path: impl AsRef<Path>) -> Result<(Vec<f64>, Vec<f64>)> {
    let path = path.as_ref();
    let file = std::fs::File::open(path).map_err(Error::Io)?;
    let mut bytes = Vec::new();
    file.take(MAX_SWATH_WINDOW_FILE_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(Error::Io)?;
    if bytes.len() as u64 > MAX_SWATH_WINDOW_FILE_BYTES {
        return Err(Error::InvalidRange(format!(
            "SWATH window file {} exceeds {MAX_SWATH_WINDOW_FILE_BYTES} bytes",
            path.display()
        )));
    }
    let content = String::from_utf8(bytes).map_err(|_| {
        Error::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("SWATH window file {} is not UTF-8", path.display()),
        ))
    })?;
    parse_swath_windows(&content)
}

/// Overwrite the precursor boundaries of `swath_maps` with the windows read
/// from `path`.
///
/// Source `annotateSwathMapsFromFile`. The non-MS1 entries of `swath_maps`
/// are matched to the file rows **by position**, so both must be in the same
/// order (usually ascending m/z); with `do_sort` the maps are first sorted by
/// ascending `upper`. Each matched map's `lower` and `upper` are replaced by
/// the row's. A row should narrow the map's window from the raw data
/// (`map.lower <= row.lower` and `row.upper <= map.upper`); a row that
/// extends beyond it is refused, or, with `force`, applied anyway and
/// reported in the returned list (the source logs it at error level).
///
/// Returns the messages of the forced overrides, empty when there were none.
///
/// The sort is stable: maps with equal `upper` keep their input order, where
/// the source's `std::sort` leaves it unspecified.
///
/// # Errors
///
/// All errors leave `swath_maps` unchanged; the source may have sorted and
/// partly annotated them before it throws.
///
/// - [`Error::InvalidValue`] when the number of non-MS1 maps differs from the
///   number of rows (the source's `IllegalArgument`), when a row extends
///   beyond its map without `force` (`IllegalArgument`), or when `do_sort`
///   meets a NaN `upper` (undefined behaviour in the source's sort).
/// - The errors of [`read_swath_windows`].
pub fn annotate_swath_maps_from_file(
    path: impl AsRef<Path>,
    swath_maps: &mut Vec<SwathMap>,
    do_sort: bool,
    force: bool,
) -> Result<Vec<String>> {
    let (lower_bounds, upper_bounds) = read_swath_windows(path)?;
    annotate_swath_maps(&lower_bounds, &upper_bounds, swath_maps, do_sort, force)
}

/// [`annotate_swath_maps_from_file`] with windows already read, for example
/// by [`parse_swath_windows`].
///
/// Native entry point; the source only annotates from a file.
///
/// # Errors
///
/// As [`annotate_swath_maps_from_file`], without the file errors, and
/// [`Error::InvalidValue`] when `lower_bounds` and `upper_bounds` differ in
/// length.
pub fn annotate_swath_maps(
    lower_bounds: &[f64],
    upper_bounds: &[f64],
    swath_maps: &mut Vec<SwathMap>,
    do_sort: bool,
    force: bool,
) -> Result<Vec<String>> {
    if lower_bounds.len() != upper_bounds.len() {
        return Err(Error::InvalidValue(format!(
            "{} lower but {} upper SWATH window boundaries",
            lower_bounds.len(),
            upper_bounds.len()
        )));
    }
    let mut maps = swath_maps.clone();
    if do_sort {
        if maps.iter().any(|m| m.upper.is_nan()) {
            return Err(Error::InvalidValue(
                "cannot sort SWATH maps with a NaN upper boundary".into(),
            ));
        }
        maps.sort_by(|a, b| a.upper.partial_cmp(&b.upper).unwrap_or(Ordering::Equal));
    }
    let mismatch = || {
        Error::InvalidValue(
            "The number of SWATH maps read from the raw data and from the annotation file do not match."
                .into(),
        )
    };
    let mut warnings = Vec::new();
    let mut j = 0usize;
    for map in maps.iter_mut().filter(|m| !m.ms1) {
        if j >= lower_bounds.len() {
            return Err(mismatch());
        }
        let (lower, upper) = (lower_bounds[j], upper_bounds[j]);
        let narrows = map.lower <= lower && upper <= map.upper;
        if !narrows {
            let message = format!(
                "SWATH window #{} from swath_windows_file extends beyond the Swath window of the data. \
                 Did you forget to apply the sort_swath_maps flag? (override with -force)",
                j + 1
            );
            if !force {
                return Err(Error::InvalidValue(message));
            }
            warnings.push(message);
        }
        map.lower = lower;
        map.upper = upper;
        j += 1;
    }
    if j != upper_bounds.len() {
        return Err(mismatch());
    }
    *swath_maps = maps;
    Ok(warnings)
}

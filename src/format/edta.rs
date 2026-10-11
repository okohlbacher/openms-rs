// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Enhanced DTA (`.edta`) feature tables, ported from `FORMAT/EDTAFile.h`.
//!
//! Input is a text file of tab, space or comma separated columns. The separator
//! is chosen from the first line, in that order: a tab anywhere in the first
//! line selects tab, otherwise a space selects space, otherwise a comma selects
//! comma, and a line with none of them is split on space.
//!
//! Three variants of the format are read:
//!
//! - Columns RT, m/z, intensity. A header line is optional.
//! - Columns RT, m/z, intensity, charge, then any number of meta-data columns.
//!   A header line is mandatory, and it names the meta-data columns.
//!
//!   ```text
//!   RT m/z Intensity charge mymeta1 mymeta2
//!   321 405.233 24543534 2 lala  lili
//!   321 406.207 4343344  2 blubb blabb
//!   ```
//!
//! - One or more (RT, m/z, intensity, charge) quadruplets, then meta-data
//!   columns. A header line is mandatory. The first quadruplet is the
//!   consensus; every following one describes a sub-feature. This variant is
//!   told from the previous one by the name of the fifth column, which must be
//!   `RT1` in any letter case. All other sub-feature column names are ignored.
//!
//!   ```text
//!   RT MZ INT CHARGE RT1 MZ1 INT1 CHARGE1 RT2 MZ2 INT2 CHARGE2
//!   321 405 100 2 321 405 100 2 321 406 50 2
//!   323 406 200 2 323 406 200 2 323 407 100 2 323 407 50 2
//!   ```
//!
//! Every variant loads into a [`ConsensusMap`](crate::kernel::ConsensusMap). Writing is available for a
//! [`ConsensusMap`](crate::kernel::ConsensusMap) (consensus plus sub-feature quadruplets) and for a
//! [`FeatureMap`](crate::kernel::FeatureMap) (the four consensus columns only). Neither writer emits meta
//! data, unique IDs, qualities, widths or identifications, exactly as the
//! source; the support document `docs/EDTA_FILE_SUPPORT.md` lists every
//! difference from the source.

use crate::concept::UniqueIdGenerator;
use crate::data_structures::string_utils::{to_double, to_int32};
use crate::format::file_info::text_format::{to_str, to_str_f32};
use crate::format::file_types::{FileType, has_valid_extension};
use crate::format::text::{self, TextFile};
use crate::kernel::{
    ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

/// The value the source's private `checkedToDouble_` / `checkedToInt_` return
/// for a missing column or an `NA` cell.
const MISSING: f64 = -1.0;

/// File adapter for Enhanced DTA files (source `EDTAFile`).
///
/// The source class is stateless apart from its constructor and destructor;
/// so is this one. Reading takes a caller-owned [`UniqueIdGenerator`], which
/// replaces the process-wide singleton the source's `setUniqueId()` draws from:
/// seed it to obtain a reproducible sequence of feature IDs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EDTAFile;

/// The variant of an input file, decided from its first line (source anonymous
/// `enum { TYPE_UNDEFINED, TYPE_OLD_NOCHARGE, TYPE_OLD_CHARGE, TYPE_CONSENSUS }`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum InputType {
    /// Three columns: RT, m/z, intensity. The charge is never set.
    OldNoCharge,
    /// RT, m/z, intensity, charge and meta-data columns.
    OldCharge,
    /// Quadruplets for the consensus and each sub-feature, then meta data.
    Consensus,
}

impl EDTAFile {
    /// A stateless adapter (source default constructor).
    pub fn new() -> Self {
        Self
    }

    /// Loads an EDTA file into a consensus map.
    ///
    /// Equivalent to [`Self::load_with_limits`] with the default
    /// [`text::Limits`].
    ///
    /// # Errors
    ///
    /// See [`Self::load_with_limits`].
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        self.load_with_limits(path, text::Limits::default(), ids)
    }

    /// Loads an EDTA file into a consensus map, reading at most what `limits`
    /// admits.
    ///
    /// The file is read as plain text (the source's `TextFile`), so a
    /// compressed file is not decompressed. The column headers recorded in the
    /// map carry the path as given, converted lossily to UTF-8.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be opened, where the source
    /// throws `Exception::FileNotFound`; [`Error::InvalidValue`] when the file
    /// exceeds `limits` or is not UTF-8; and the parse errors of
    /// [`Self::read_with_limits`].
    pub fn load_with_limits(
        &self,
        path: impl AsRef<Path>,
        limits: text::Limits,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        let path = path.as_ref();
        let options = text::ReadOptions {
            limits,
            ..text::ReadOptions::default()
        };
        let input = TextFile::from_path(path, &options)?;
        parse(&input, &path.to_string_lossy(), ids)
    }

    /// Reads EDTA text from `reader`; `filename` is recorded in the map's
    /// column headers, as the source records the path it loaded.
    ///
    /// # Errors
    ///
    /// See [`Self::read_with_limits`].
    pub fn read(
        &self,
        reader: impl BufRead,
        filename: &str,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        self.read_with_limits(reader, filename, text::Limits::default(), ids)
    }

    /// Reads EDTA text from `reader` under explicit `limits`.
    ///
    /// Unique IDs are drawn from `ids` in the source's order: one for the map,
    /// then for each data line one for the consensus feature followed by one
    /// for every sub-feature quadruplet, including a quadruplet that is
    /// skipped because its RT or m/z is missing.
    ///
    /// A line that is empty after trimming space, tab, CR and LF is skipped.
    /// Numbers are converted as the source's `StringUtils::toDouble` and
    /// `StringUtils::toInt32`, and a cell that is exactly `NA` (or a column
    /// that is absent) stands for -1. A sub-feature is kept only when both its
    /// RT and its m/z differ from -1. Meta-data cells that are empty after
    /// trimming are skipped; the others are stored as string meta values under
    /// their trimmed column name, a later column overwriting an earlier one of
    /// the same name.
    ///
    /// The map has one column header per sub-feature column (at least one),
    /// labelled `EDTA_Map <index>` and each recording `filename` and the number
    /// of consensus features.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] with the source's messages and one-based line
    /// numbers when:
    ///
    /// - the first line has fewer than three columns (also for an empty input,
    ///   which the source dereferences without checking);
    /// - the first line holds numbers but more than three columns, since only
    ///   the three-column variant may omit its header;
    /// - a data line has fewer than three columns;
    /// - one of the consensus cells, including a fourth cell of a three-column
    ///   file, which the source converts even though it never uses it, is not
    ///   a number;
    /// - one of a sub-feature's four cells is not a number;
    /// - a non-empty meta-data cell has no (or an empty) column name.
    ///
    /// As a native check the source lacks, it also returns [`Error::Parse`]
    /// when a converted value is NaN or infinite, or an intensity overflows
    /// `f32`: the source stores such values and every later consumer of the
    /// map rejects them. [`Error::InvalidValue`] is returned when the input
    /// exceeds `limits` or is not UTF-8.
    pub fn read_with_limits(
        &self,
        reader: impl BufRead,
        filename: &str,
        limits: text::Limits,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ConsensusMap> {
        let options = text::ReadOptions {
            limits,
            ..text::ReadOptions::default()
        };
        let input = TextFile::from_reader(reader, &options)?;
        parse(&input, filename, ids)
    }

    /// Stores a consensus map as an enhanced DTA file.
    ///
    /// The header is `RT`, `m/z`, `intensity`, `charge`, followed by `RT<i>`,
    /// `m/z<i>`, `intensity<i>`, `charge<i>` for `i` from 1 to the largest
    /// number of sub-features of any consensus feature, all tab separated.
    /// Each line holds the consensus position, intensity and charge, then its
    /// sub-feature handles in handle order (map index, then unique ID), then
    /// `NA` quadruplets up to the column count. Numbers are written in the
    /// source's full-precision `StringUtils::toStr` form (intensities as
    /// `float`). Lines end in LF on every platform, where the source's
    /// text-mode stream writes CRLF on Windows.
    ///
    /// The source header documents this overload as "NOT IMPLEMENTED"; the
    /// source implementation nevertheless writes the file described here, and
    /// this port follows the implementation.
    ///
    /// # Information the format cannot hold
    ///
    /// Handles are written in consecutive columns, so a consensus feature whose
    /// handles come from maps 0 and 2 reads back with sub-features from maps 0
    /// and 1: the map index is not preserved. Unique IDs, meta values, widths,
    /// qualities, column headers and identifications are not written. This is
    /// the source's behaviour.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `path` is not UTF-8 or names a
    /// known file type other than EDTA, where the source throws
    /// `Exception::UnableToCreateFile` ("invalid file extension"); an unknown
    /// extension is accepted, as in the source. Returns [`Error::Io`] when the
    /// file cannot be created or written (the source's `UnableToCreateFile`).
    pub fn store_consensus(&self, path: impl AsRef<Path>, map: &ConsensusMap) -> Result<()> {
        let path = path.as_ref();
        let name = path.to_str().ok_or_else(|| {
            Error::InvalidValue(format!("EDTA output path is not UTF-8: {}", path.display()))
        })?;
        if !has_valid_extension(name, FileType::Edta) {
            return Err(Error::InvalidValue(format!(
                "cannot create file '{name}': invalid file extension, expected '{}'",
                FileType::Edta.name()
            )));
        }
        store_with(path, |writer| write_consensus_lines(writer, map))
    }

    /// Writes a consensus map in the layout of [`Self::store_consensus`] to
    /// `writer`, without the file-extension check.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when writing fails.
    pub fn write_consensus(&self, writer: impl Write, map: &ConsensusMap) -> Result<()> {
        let mut writer = writer;
        write_consensus_lines(&mut writer, map)?;
        writer.flush()?;
        Ok(())
    }

    /// Stores a feature map as an enhanced DTA file with the tab-separated
    /// columns `RT`, `m/z`, `intensity` and `charge`.
    ///
    /// Unlike [`Self::store_consensus`], the source does not check the file
    /// extension here, and neither does this port.
    ///
    /// # Information the format cannot hold
    ///
    /// Only the four columns are written; subordinates, convex hulls,
    /// qualities, meta values, identifications and unique IDs are dropped, as
    /// in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be created or written (the
    /// source's `Exception::UnableToCreateFile`).
    pub fn store_features(&self, path: impl AsRef<Path>, map: &FeatureMap) -> Result<()> {
        store_with(path.as_ref(), |writer| write_feature_lines(writer, map))
    }

    /// Writes a feature map in the layout of [`Self::store_features`] to
    /// `writer`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when writing fails.
    pub fn write_features(&self, writer: impl Write, map: &FeatureMap) -> Result<()> {
        let mut writer = writer;
        write_feature_lines(&mut writer, map)?;
        writer.flush()?;
        Ok(())
    }
}

/// Source `StringUtils::split(s, char, out)`: an empty string yields no
/// fields; otherwise every separator ends a field, so consecutive separators
/// yield empty fields.
fn split(line: &str, separator: char) -> Vec<&str> {
    if line.is_empty() {
        Vec::new()
    } else {
        line.split(separator).collect()
    }
}

/// Source private `checkedToDouble_`: the column's value, or -1 when the column
/// is absent or exactly `NA`. `Ok(None)` is a conversion failure.
fn checked_double(parts: &[&str], index: usize) -> Option<f64> {
    match parts.get(index) {
        Some(&cell) if cell != "NA" => to_double(cell).ok(),
        _ => Some(MISSING),
    }
}

/// Source private `checkedToInt_`: the column's value, or -1 when the column is
/// absent or exactly `NA`.
fn checked_int(parts: &[&str], index: usize) -> Option<i32> {
    match parts.get(index) {
        Some(&cell) if cell != "NA" => to_int32(cell).ok(),
        _ => Some(-1),
    }
}

/// One parsed quadruplet. `None` is the source's conversion exception.
fn quadruplet(parts: &[&str], first: usize) -> Option<(f64, f64, f64, i32)> {
    let rt = checked_double(parts, first)?;
    let mz = checked_double(parts, first + 1)?;
    let intensity = checked_double(parts, first + 2)?;
    let charge = checked_int(parts, first + 3)?;
    Some((rt, mz, intensity, charge))
}

/// Native check: the three coordinates finite and the intensity within `f32`.
fn finite(rt: f64, mz: f64, intensity: f64, line: usize, what: &str) -> Result<f32> {
    let narrowed = intensity as f32;
    if rt.is_finite() && mz.is_finite() && narrowed.is_finite() {
        Ok(narrowed)
    } else {
        Err(Error::Parse {
            line,
            message: format!(
                "{what} holds a value that is NaN, infinite or outside the f32 intensity range"
            ),
        })
    }
}

fn parse(input: &TextFile, filename: &str, ids: &mut UniqueIdGenerator) -> Result<ConsensusMap> {
    let lines = input.lines();
    let first = lines.first().map(String::as_str).unwrap_or("");
    let mut map = ConsensusMap::new();
    map.unique_id = ids.get_unique_id();

    let separator = if first.contains('\t') {
        '\t'
    } else if first.contains(' ') {
        ' '
    } else if first.contains(',') {
        ','
    } else {
        ' '
    };
    let headers: Vec<&str> = split(first, separator)
        .into_iter()
        .map(text::trim)
        .collect();
    let header_trimmed = text::trim(first);

    let mut input_type = match headers.len() {
        0..=2 => {
            return Err(Error::Parse {
                line: 1,
                message: format!(
                    "Failed parsing in line 1: not enough columns! Expected at least 3 columns!\nOffending line: '{header_trimmed}'  (line 1)\n"
                ),
            });
        }
        3 => InputType::OldNoCharge,
        // The source leaves the type undefined for five or more columns until
        // the fifth column is inspected below.
        _ => InputType::OldCharge,
    };
    // The source converts the first three cells; any failure marks a header.
    let has_header = headers[..3].iter().any(|cell| to_double(cell).is_err());
    let offset = usize::from(has_header);
    if headers.len() >= 5 {
        input_type = if headers[4].eq_ignore_ascii_case("RT1") {
            InputType::Consensus
        } else {
            InputType::OldCharge
        };
    }
    let input_features = if input_type == InputType::Consensus {
        headers.len() / 4
    } else {
        1
    };
    if !has_header && input_type != InputType::OldNoCharge {
        return Err(Error::Parse {
            line: 1,
            message: format!(
                "Failed parsing in line 1: No HEADER provided. This is only allowed for three columns. You have more!\nOffending line: '{header_trimmed}'  (line 1)\n"
            ),
        });
    }

    let mut features = Vec::new();
    features
        .try_reserve(lines.len())
        .map_err(|_| Error::InvalidValue("EDTA feature allocation failed".into()))?;
    for (index, raw) in lines.iter().enumerate().skip(offset) {
        let line = index + 1;
        let line_trimmed = text::trim(raw);
        if line_trimmed.is_empty() {
            // The source logs "Notice: Empty line ignored" for every empty
            // line but the last; this port emits no log output.
            continue;
        }
        let parts = split(raw, separator);
        if parts.len() < 3 {
            return Err(Error::Parse {
                line,
                message: format!(
                    "Failed parsing in line {line}: At least three columns are needed! (got  {})\nOffending line: '{line_trimmed}'  (line {line})\n",
                    parts.len()
                ),
            });
        }

        let mut consensus = ConsensusFeature::new();
        consensus.base.unique_id = ids.get_unique_id();
        let (rt, mz, intensity, charge) = quadruplet(&parts, 0).ok_or_else(|| Error::Parse {
            line,
            message: format!(
                "Failed parsing in line {line}: Could not convert the first three columns to a number!\nOffending line: '{line_trimmed}'  (line {line})\n"
            ),
        })?;
        consensus.base.intensity = finite(rt, mz, intensity, line, "the consensus quadruplet")?;
        consensus.base.rt = rt;
        consensus.base.mz = mz;
        if input_type != InputType::OldNoCharge {
            consensus.base.charge = charge;
        }

        for sub in 1..input_features {
            let unique_id = ids.get_unique_id();
            let (rt, mz, intensity, charge) =
                quadruplet(&parts, sub * 4).ok_or_else(|| Error::Parse {
                    line,
                    message: format!(
                        "Failed parsing in line {line}: Could not convert one of the four sub-feature columns (starting at column {}) to a number! Is the correct separator specified?\nOffending line: '{line_trimmed}'  (line {line})\n",
                        sub * 4 + 1
                    ),
                })?;
            // Only features with at least RT and m/z set are accepted.
            if rt != MISSING && mz != MISSING {
                let mut feature = Feature::new(rt, mz, 0.0);
                feature.base.intensity =
                    finite(rt, mz, intensity, line, "a sub-feature quadruplet")?;
                feature.base.charge = charge;
                feature.base.unique_id = unique_id;
                consensus.insert(FeatureHandle::new((sub - 1) as u64, &feature.base))?;
            }
        }

        for (column, cell) in parts.iter().enumerate().skip(input_features * 4) {
            let value = text::trim(cell);
            if value.is_empty() {
                continue;
            }
            let name = headers.get(column).copied().unwrap_or("");
            if name.is_empty() {
                return Err(Error::Parse {
                    line,
                    message: format!(
                        "Error: Missing meta data header for column {}!Offending header line: '{header_trimmed}'  (line 1)",
                        column + 1
                    ),
                });
            }
            consensus
                .base
                .metadata
                .insert(name.to_owned(), MetaValue::from(value));
        }
        features.push(consensus);
    }

    let count = features.len();
    map.features = features;
    let maps = input_features.saturating_sub(1).max(1);
    for index in 0..maps {
        map.column_headers.insert(
            index as u64,
            ColumnHeader {
                filename: filename.to_owned(),
                label: format!("EDTA_Map {index}"),
                size: count,
                ..ColumnHeader::default()
            },
        );
    }
    Ok(map)
}

fn write_consensus_lines(writer: &mut impl Write, map: &ConsensusMap) -> Result<()> {
    let max_sub = map
        .features
        .iter()
        .map(ConsensusFeature::len)
        .max()
        .unwrap_or(0);
    let mut header = String::from("RT\tm/z\tintensity\tcharge");
    for i in 1..=max_sub {
        header.push_str(&format!("\tRT{i}\tm/z{i}\tintensity{i}\tcharge{i}"));
    }
    writeln!(writer, "{header}")?;
    for feature in &map.features {
        let base = &feature.base;
        write!(
            writer,
            "{}\t{}\t{}\t{}",
            to_str(base.rt),
            to_str(base.mz),
            to_str_f32(base.intensity),
            base.charge
        )?;
        for handle in feature.handles() {
            write!(
                writer,
                "\t{}\t{}\t{}\t{}",
                to_str(handle.rt),
                to_str(handle.mz),
                to_str_f32(handle.intensity),
                handle.charge
            )?;
        }
        for _ in feature.len()..max_sub {
            writer.write_all(b"\tNA\tNA\tNA\tNA")?;
        }
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn write_feature_lines(writer: &mut impl Write, map: &FeatureMap) -> Result<()> {
    writer.write_all(b"RT\tm/z\tintensity\tcharge\n")?;
    for feature in &map.features {
        let base = &feature.base;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            to_str(base.rt),
            to_str(base.mz),
            to_str_f32(base.intensity),
            base.charge
        )?;
    }
    Ok(())
}

/// Creates `path` and streams `body` into it, as the source's `TextFile::store`
/// truncates and writes the destination.
fn store_with(path: &Path, body: impl FnOnce(&mut BufWriter<File>) -> Result<()>) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    body(&mut writer)?;
    writer.flush()?;
    Ok(())
}

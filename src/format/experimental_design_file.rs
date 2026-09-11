// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The tab-separated experimental design reader, in both source table layouts.

pub use super::text::Limits;
use super::text::{self, TextFile};
use crate::metadata::{ExperimentalDesign, MSFileSectionEntry, SampleSection};
use crate::system::file;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

const FRACTION_GROUP: &str = "Fraction_Group";
const FRACTION: &str = "Fraction";
const SPECTRA_FILEPATH: &str = "Spectra_Filepath";
const LABEL: &str = "Label";
const SAMPLE: &str = "Sample";
/// Columns of the MS file section, excluded from the sample metadata factors.
const FILE_COLUMNS: [&str; 4] = [FRACTION_GROUP, FRACTION, SPECTRA_FILEPATH, LABEL];

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadOptions {
    /// Fail when a resolved `Spectra_Filepath` does not exist.
    pub require_spectra_files: bool,
    pub limits: Limits,
}

fn parse_error(line: usize, message: impl Into<String>) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}

/// Read a design from a tab-separated file, discarding advisory warnings.
pub fn load(path: impl AsRef<Path>, options: &ReadOptions) -> Result<ExperimentalDesign> {
    Ok(load_with_warnings(path, options)?.0)
}

/// Read a design and return the source's advisory warnings. The source writes
/// these to its log stream; this port returns them so a caller owns them.
pub fn load_with_warnings(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<(ExperimentalDesign, Vec<String>)> {
    let path = path.as_ref();
    let text = TextFile::from_path(
        path,
        &text::ReadOptions {
            trim_lines: true,
            limits: options.limits,
            ..Default::default()
        },
    )?;
    let design_path = path.to_string_lossy().into_owned();
    load_text(&text, &design_path, options)
}

/// Parse an already-loaded line buffer. `design_path` locates relative spectra
/// paths and is the source's `filename` diagnostic argument; it is not read.
pub fn load_text(
    text: &TextFile,
    design_path: &str,
    options: &ReadOptions,
) -> Result<(ExperimentalDesign, Vec<String>)> {
    if is_one_table(text) {
        parse_one_table(text, design_path, options)
    } else {
        parse_two_table(text, design_path, options)
    }
}

/// The source format detector: it scans every line rather than only headers,
/// does not trim cells and does not skip comment lines. A file is two-table as
/// soon as one line has no `Fraction_Group` cell and exactly one `Sample` cell.
fn is_one_table(text: &TextFile) -> bool {
    for line in text.iter() {
        if line.trim().is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.trim().split('\t').collect();
        if !cells.contains(&FRACTION_GROUP)
            && cells.iter().filter(|cell| **cell == SAMPLE).count() == 1
        {
            return false;
        }
    }
    true
}

/// Header columns to their indices, with the source's duplicate, missing and
/// (optionally) unexpected column checks.
fn parse_header(
    cells: &[String],
    line: usize,
    required: &[&str],
    optional: &[&str],
    allow_other: bool,
) -> Result<BTreeMap<String, usize>> {
    let mut columns = BTreeMap::new();
    for (index, name) in cells.iter().enumerate() {
        if columns.insert(name.clone(), index).is_some() {
            return Err(parse_error(
                line,
                "some column headers of the table appear multiple times",
            ));
        }
    }
    for name in required {
        if !columns.contains_key(*name) {
            return Err(parse_error(line, format!("missing column header: {name}")));
        }
    }
    if !allow_other {
        for name in columns.keys() {
            if !required.contains(&name.as_str()) && !optional.contains(&name.as_str()) {
                return Err(parse_error(
                    line,
                    format!("header not allowed in this section of the design: {name}"),
                ));
            }
        }
    }
    Ok(columns)
}

fn cell<'a>(
    cells: &'a [String],
    columns: &BTreeMap<String, usize>,
    name: &str,
    line: usize,
) -> Result<&'a str> {
    let index = columns
        .get(name)
        .ok_or_else(|| parse_error(line, format!("missing column {name}")))?;
    cells
        .get(*index)
        .map(String::as_str)
        .ok_or_else(|| parse_error(line, "wrong number of records in line"))
}

/// Source `toInt32` followed by an unsigned assignment. Negative values are
/// rejected instead of wrapping (see `OpenMS_CPP_ISSUES.md`, CPP-060).
fn index_cell(
    cells: &[String],
    columns: &BTreeMap<String, usize>,
    name: &str,
    line: usize,
) -> Result<u32> {
    let text = cell(cells, columns, name, line)?;
    text.parse::<i32>()
        .ok()
        .and_then(|value| u32::try_from(value).ok())
        .ok_or_else(|| parse_error(line, format!("{name} must be a nonnegative integer")))
}

fn split_trimmed(line: &str) -> Vec<String> {
    line.split('\t')
        .map(|cell| cell.trim().to_owned())
        .collect()
}

/// Resolve a spectra path: an absolute path is kept; a relative one is tried
/// against the design file's directory, then the working directory, and is
/// otherwise kept as written.
fn find_spectra_file(
    spec_file: &str,
    design_path: &str,
    require: bool,
    line: usize,
) -> Result<String> {
    let mut result = spec_file.to_owned();
    if Path::new(spec_file).is_relative() {
        let relative = format!("{}/{}", file::path(design_path), spec_file);
        if file::exists(&relative) {
            result = relative;
        } else if let Ok(absolute) = file::absolute_path(spec_file)
            && file::exists(&absolute)
        {
            result = absolute.to_string_lossy().into_owned();
        }
    }
    if require && !file::exists(&result) {
        return Err(parse_error(
            line,
            format!("spectra file does not exist: '{result}'"),
        ));
    }
    Ok(result)
}

fn parse_one_table(
    text: &TextFile,
    design_path: &str,
    options: &ReadOptions,
) -> Result<(ExperimentalDesign, Vec<String>)> {
    let mut warnings = Vec::new();
    let mut msfile_section = Vec::new();
    let mut content: Vec<Vec<String>> = Vec::new();
    let mut sample_columns: BTreeMap<String, usize> = BTreeMap::new();
    let mut columns: BTreeMap<String, usize> = BTreeMap::new();
    let mut samplename_to_index: BTreeMap<String, usize> = BTreeMap::new();
    let (mut has_sample, mut has_label) = (false, false);
    let mut header_seen = false;
    let mut n_col = 0;

    for (index, raw) in text.iter().enumerate() {
        let number = index + 1;
        let line = raw.trim();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut cells = split_trimmed(line);
        if !header_seen {
            header_seen = true;
            columns = parse_header(
                &cells,
                number,
                &[FRACTION_GROUP, FRACTION, SPECTRA_FILEPATH],
                &[LABEL, SAMPLE],
                true,
            )?;
            has_label = columns.contains_key(LABEL);
            has_sample = columns.contains_key(SAMPLE);
            // Absent optional columns are appended to the header, as in source.
            if !has_label {
                columns.insert(LABEL.into(), columns.len());
                cells.push(LABEL.into());
            }
            if !has_sample {
                columns.insert(SAMPLE.into(), columns.len());
                cells.push(SAMPLE.into());
            }
            n_col = columns.len();
            // Every other column is sample metadata; Sample itself is a factor.
            for (position, name) in cells.iter().enumerate() {
                if !FILE_COLUMNS.contains(&name.as_str()) {
                    sample_columns.insert(name.clone(), position);
                }
            }
            continue;
        }

        if !has_label {
            cells.push("1".into());
        }
        let label = index_cell(&cells, &columns, LABEL, number)?;
        let fraction = index_cell(&cells, &columns, FRACTION, number)?;
        let fraction_group = index_cell(&cells, &columns, FRACTION_GROUP, number)?;
        if !has_sample {
            if label > 1 {
                return Err(parse_error(
                    number,
                    "column 'Sample' is required for multiplexed one-table designs (Label > 1)",
                ));
            }
            // Without a Sample column each fraction group becomes its own sample.
            cells.push(fraction_group.to_string());
        }
        let sample_name = cell(&cells, &columns, SAMPLE, number)?.to_owned();
        if n_col != cells.len() {
            return Err(parse_error(number, "wrong number of records in line"));
        }
        let next = samplename_to_index.len();
        let (sample, inserted) = match samplename_to_index.get(&sample_name) {
            Some(existing) => (*existing, false),
            None => {
                samplename_to_index.insert(sample_name.clone(), next);
                (next, true)
            }
        };
        let mut sample_cells = Vec::with_capacity(sample_columns.len());
        for position in sample_columns.values() {
            sample_cells.push(
                cells
                    .get(*position)
                    .cloned()
                    .ok_or_else(|| parse_error(number, "wrong number of records in line"))?,
            );
        }
        if inserted {
            content.push(sample_cells);
        } else if content[sample] != sample_cells {
            warnings.push(format!(
                "factors for sample '{sample_name}' do not match those of its first row"
            ));
        }

        msfile_section.push(MSFileSectionEntry {
            fraction_group,
            fraction,
            label,
            sample: u32::try_from(sample)
                .map_err(|_| parse_error(number, "too many samples in design"))?,
            sample_name,
            path: find_spectra_file(
                cell(&cells, &columns, SPECTRA_FILEPATH, number)?,
                design_path,
                options.require_spectra_files,
                number,
            )?,
        });
    }

    // The stored sample columns are renumbered to their position in the
    // alphabetically ordered sample section, matching the content rows above.
    for (position, index) in sample_columns.values_mut().enumerate() {
        *index = position;
    }
    let section = SampleSection::from_table(content, samplename_to_index, sample_columns)?;
    Ok((
        ExperimentalDesign::from_sections(msfile_section, section)?,
        warnings,
    ))
}

fn parse_two_table(
    text: &TextFile,
    design_path: &str,
    options: &ReadOptions,
) -> Result<(ExperimentalDesign, Vec<String>)> {
    #[derive(PartialEq)]
    enum State {
        RunHeader,
        RunContent,
        SampleHeader,
        SampleContent,
    }

    let mut msfile_section = Vec::new();
    let mut content: Vec<Vec<String>> = Vec::new();
    let mut sample_to_rowindex: BTreeMap<String, usize> = BTreeMap::new();
    let mut sample_columns: BTreeMap<String, usize> = BTreeMap::new();
    let mut columns: BTreeMap<String, usize> = BTreeMap::new();
    let (mut has_sample, mut has_label) = (false, false);
    let mut state = State::RunHeader;
    let mut n_col = 0;

    for (index, raw) in text.iter().enumerate() {
        let number = index + 1;
        let line = raw.trim();
        // Empty lines separate the two sections, so they are only skipped
        // outside the file section.
        if line.starts_with('#') || (line.is_empty() && state != State::RunContent) {
            continue;
        }
        let cells = split_trimmed(line);
        match state {
            State::RunHeader => {
                state = State::RunContent;
                columns = parse_header(
                    &cells,
                    number,
                    &[FRACTION_GROUP, FRACTION, SPECTRA_FILEPATH],
                    &[LABEL, SAMPLE],
                    false,
                )?;
                has_label = columns.contains_key(LABEL);
                has_sample = columns.contains_key(SAMPLE);
                n_col = columns.len();
            }
            State::RunContent if line.is_empty() => state = State::SampleHeader,
            State::RunContent => {
                if n_col != cells.len() {
                    return Err(parse_error(number, "wrong number of records in line"));
                }
                let fraction_group = index_cell(&cells, &columns, FRACTION_GROUP, number)?;
                msfile_section.push(MSFileSectionEntry {
                    fraction_group,
                    fraction: index_cell(&cells, &columns, FRACTION, number)?,
                    label: if has_label {
                        index_cell(&cells, &columns, LABEL, number)?
                    } else {
                        1
                    },
                    sample: 0,
                    sample_name: if has_sample {
                        cell(&cells, &columns, SAMPLE, number)?.to_owned()
                    } else {
                        format!("Fraction group {fraction_group}")
                    },
                    path: find_spectra_file(
                        cell(&cells, &columns, SPECTRA_FILEPATH, number)?,
                        design_path,
                        options.require_spectra_files,
                        number,
                    )?,
                });
            }
            State::SampleHeader => {
                state = State::SampleContent;
                sample_columns = parse_header(&cells, number, &[SAMPLE], &[], true)?;
                n_col = sample_columns.len();
            }
            State::SampleContent => {
                if n_col != cells.len() {
                    return Err(parse_error(number, "wrong number of records in line"));
                }
                let sample = cell(&cells, &sample_columns, SAMPLE, number)?.to_owned();
                if sample_to_rowindex.contains_key(&sample) {
                    return Err(parse_error(
                        number,
                        format!("sample '{sample}' appears multiple times in the sample table"),
                    ));
                }
                sample_to_rowindex.insert(sample, content.len());
                content.push(cells);
            }
        }
    }

    for row in &mut msfile_section {
        let index = sample_to_rowindex.get(&row.sample_name).ok_or_else(|| {
            Error::InvalidValue(format!(
                "the sample section has no sample named '{}'",
                row.sample_name
            ))
        })?;
        row.sample =
            u32::try_from(*index).map_err(|_| Error::InvalidValue("too many samples".into()))?;
    }
    let section = SampleSection::from_table(content, sample_to_rowindex, sample_columns)?;
    Ok((
        ExperimentalDesign::from_sections(msfile_section, section)?,
        Vec::new(),
    ))
}

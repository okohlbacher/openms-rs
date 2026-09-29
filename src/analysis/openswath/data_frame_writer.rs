// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Row-by-row writers for named numeric tables (score matrices).
//!
//! Ports `OPENSWATHALGO/DATAACCESS/DataFrameWriter.h` with
//! `DataFrameWriter.cpp`: the interface `IDataFrameWriter` becomes the trait
//! [`DataFrameWriter`](crate::analysis::openswath::data_frame_writer::DataFrameWriter), `DataMatrix` an in-memory table, and `CSVWriter` a
//! tab-separated file writer. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_ALGO_SUPPORT.md).

use crate::{Error, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// A sink for a table with named columns and named rows. Source `IDataFrameWriter`.
pub trait DataFrameWriter {
    /// Set (or write) the column names.
    ///
    /// # Errors
    ///
    /// Implementation-specific; a file writer reports I/O errors.
    fn colnames(&mut self, colnames: &[String]) -> Result<()>;

    /// Store one row of values under `rowname`.
    ///
    /// # Errors
    ///
    /// Implementation-specific; a file writer reports I/O errors.
    fn store(&mut self, rowname: &str, values: &[f64]) -> Result<()>;
}

/// An in-memory table. Source `DataMatrix`.
///
/// The source keeps the column names, row names and rows private and offers
/// no way to read them back; the port adds read accessors, since a table
/// nobody can read is only useful as a discarding sink.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DataMatrix {
    colnames: Vec<String>,
    rownames: Vec<String>,
    store: Vec<Vec<f64>>,
}

impl DataMatrix {
    /// An empty table.
    pub fn new() -> Self {
        Self::default()
    }

    /// The column names last set.
    pub fn column_names(&self) -> &[String] {
        &self.colnames
    }

    /// The row names, in storage order.
    pub fn row_names(&self) -> &[String] {
        &self.rownames
    }

    /// The rows, parallel to [`Self::row_names`].
    pub fn rows(&self) -> &[Vec<f64>] {
        &self.store
    }
}

impl DataFrameWriter for DataMatrix {
    /// Replace the column names; never fails.
    fn colnames(&mut self, colnames: &[String]) -> Result<()> {
        self.colnames = colnames.to_vec();
        Ok(())
    }

    /// Append the row; never fails. Row lengths are not checked against the
    /// column count, as in the source.
    fn store(&mut self, rowname: &str, values: &[f64]) -> Result<()> {
        self.rownames.push(rowname.to_string());
        self.store.push(values.to_vec());
        Ok(())
    }
}

/// A tab-separated text writer. Source `CSVWriter` (tab-separated despite the
/// name).
///
/// [`DataFrameWriter::colnames`] writes the names joined by tabs;
/// [`DataFrameWriter::store`] writes the row name, a tab, and the values joined
/// by tabs, each formatted as a C++ stream does after `std::setprecision(5)`
/// (`%g` with five significant digits: `1234.5678` → `1234.6`, `0.000012345`
/// → `1.2345e-05`, `3.0` → `3`). Every line ends with `\n`.
///
/// # Differences from the source
///
/// - [`Self::create`] fails when the file cannot be created; the source opens
///   the stream in its constructor without checking and then writes nothing.
/// - Write errors are reported; the source ignores them.
/// - [`Self::finish`] flushes and reports a flush error. The source destructor
///   flushes, closes, and prints "have flushed and closed the file stream" to
///   standard output; the port prints nothing. Dropping the writer without
///   `finish` still flushes, ignoring errors.
#[derive(Debug)]
pub struct CsvWriter {
    out: BufWriter<File>,
    sep: &'static str,
    eol: &'static str,
}

impl CsvWriter {
    /// Create (truncate) `filename` for writing.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created.
    pub fn create(filename: impl AsRef<Path>) -> Result<Self> {
        Ok(Self {
            out: BufWriter::new(File::create(filename)?),
            sep: "\t",
            eol: "\n",
        })
    }

    /// Flush and close the file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<()> {
        self.out.flush().map_err(Error::from)
    }
}

impl DataFrameWriter for CsvWriter {
    fn colnames(&mut self, colnames: &[String]) -> Result<()> {
        self.out.write_all(colnames.join(self.sep).as_bytes())?;
        self.out.write_all(self.eol.as_bytes())?;
        Ok(())
    }

    fn store(&mut self, rowname: &str, values: &[f64]) -> Result<()> {
        let mut line = String::from(rowname);
        line.push_str(self.sep);
        let formatted: Vec<String> = values.iter().map(|&v| ostream_g5(v)).collect();
        line.push_str(&formatted.join(self.sep));
        line.push_str(self.eol);
        self.out.write_all(line.as_bytes())?;
        Ok(())
    }
}

/// `ostream << value` after `std::setprecision(5)` in the default float
/// format: `%.5g`, trailing zeros removed, exponent at least two digits, and
/// libstdc++'s spellings of the non-finite values.
fn ostream_g5(value: f64) -> String {
    const DIGITS: i32 = 5;
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{:.*e}", (DIGITS - 1) as usize, value);
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        return scientific;
    };
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if !(-4..DIGITS).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", trim_zeros(mantissa), exponent.abs());
    }
    let decimals = usize::try_from(DIGITS - 1 - exponent).unwrap_or(0);
    trim_zeros(&format!("{value:.decimals$}")).to_string()
}

fn trim_zeros(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

#[cfg(test)]
mod tests {
    use super::ostream_g5;

    #[test]
    fn g5_matches_the_c_stream_rule() {
        // Independently derived from the C `%.5g` rule.
        assert_eq!(ostream_g5(1234.5678), "1234.6");
        assert_eq!(ostream_g5(123456.0), "1.2346e+05");
        assert_eq!(ostream_g5(0.000012345), "1.2345e-05");
        assert_eq!(ostream_g5(0.0001), "0.0001");
        assert_eq!(ostream_g5(3.0), "3");
        assert_eq!(ostream_g5(-2.5), "-2.5");
        assert_eq!(ostream_g5(0.0), "0");
        assert_eq!(ostream_g5(f64::INFINITY), "inf");
    }
}

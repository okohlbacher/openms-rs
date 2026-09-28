// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writer for the meta-value table consumed by GNPS Feature-Based Molecular
//! Networking, ported from `FORMAT/GNPSMetaValueFile.h`.
//!
//! GNPS FBMN can pick up per-sample metadata so the web platform can group or
//! colour network nodes by sample attribute. This writer produces the minimal
//! seed of such a table: one row per mzML file used to build the input
//! [`ConsensusMap`](crate::kernel::ConsensusMap), holding that file's basename and a sequential `MAP<i>`
//! identifier. Users typically extend the file with further attribute columns
//! before uploading it (see
//! <https://ccms-ucsd.github.io/GNPSDocumentation/metadata/>). The file is
//! optional for FBMN; the `GNPSExport` TOPP tool produces it. See
//! `docs/GNPS_TABLES_SUPPORT.md`.

use crate::Result;
use crate::format::sv_out_stream::{QuotingMethod, SVOutStream};
use crate::kernel::ConsensusMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Writer for the meta-value (`.tsv`) table of GNPS FBMN (source
/// `GNPSMetaValueFile`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GNPSMetaValueFile;

impl GNPSMetaValueFile {
    /// Write the FBMN meta-value TSV seed of `consensus_map` to `output_file`,
    /// which is overwritten if it exists.
    ///
    /// The layout is that of [`Self::write`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::Io`] when the file cannot be created or written.
    /// The source opens the file without checking it, so a destination that
    /// cannot be opened is silently not written; this port reports it. Also
    /// the errors of [`Self::write`].
    pub fn store(consensus_map: &ConsensusMap, output_file: impl AsRef<Path>) -> Result<()> {
        let mut writer = BufWriter::new(File::create(output_file.as_ref())?);
        Self::write(consensus_map, &mut writer)?;
        writer.flush()?;
        Ok(())
    }

    /// Write the FBMN meta-value TSV seed to `writer`.
    ///
    /// The primary MS-run paths of `consensus_map` (the filename of every
    /// column header, in key order; `ConsensusMap::getPrimaryMSRunPath`)
    /// provide the file list. The file starts with the header row
    /// `<TAB>filename<TAB>ATTRIBUTE_MAPID`, whose first column name is empty;
    /// each data row is `<index><TAB><basename><TAB>MAP<index>` with `index`
    /// counting from 0. The basename is the text after the last `/` or `\`.
    /// Fields are written through [`SVOutStream`] without quoting, so a tab
    /// inside a basename becomes `_`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when a basename contains a line
    /// break, which the source refuses with `Exception::IllegalArgument`, and
    /// [`crate::Error::Io`] when writing fails.
    pub fn write(consensus_map: &ConsensusMap, writer: impl Write) -> Result<()> {
        let mut out = SVOutStream::with_options(writer, "\t", "_", QuotingMethod::None)?;
        out.write_field("")?;
        out.write_field("filename")?;
        out.write_field("ATTRIBUTE_MAPID")?;
        out.end_line()?;
        for (index, path) in consensus_map.primary_ms_run_path().iter().enumerate() {
            let basename = path
                .rfind(['/', '\\'])
                .map_or(path.as_str(), |slash| &path[slash + 1..]);
            out.write_field(&index.to_string())?;
            out.write_field(basename)?;
            out.write_field(&format!("MAP{index}"))?;
            out.end_line()?;
        }
        out.finish()?;
        Ok(())
    }
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writer for the feature-quantification table of GNPS Feature-Based Molecular
//! Networking, ported from `FORMAT/GNPSQuantificationFile.h`.
//!
//! The format matches the schema at
//! <https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/#feature-quantification-table>.
//! The file is **required** for FBMN, in contrast to the optional meta-value
//! table of [`crate::format::gnps_meta_value`]. The `GNPSExport` TOPP tool
//! produces it. See `docs/GNPS_TABLES_SUPPORT.md`.

use crate::Result;
use crate::concept::constants::user_param::{
    IIMN_ADDUCT_PARTNERS, IIMN_ANNOTATION_NETWORK_NUMBER, IIMN_BEST_ION, IIMN_ROW_ID,
};
use crate::data_structures::list::ListFormat;
use crate::format::sv_out_stream::{QuotingMethod, SVOutStream};
use crate::kernel::{ConsensusMap, FeatureHandle};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// The ion-identity molecular networking meta values, in column order.
const IIMN_META_VALUES: [&str; 4] = [
    IIMN_ROW_ID,
    IIMN_BEST_ION,
    IIMN_ADDUCT_PARTNERS,
    IIMN_ANNOTATION_NETWORK_NUMBER,
];

/// Writer for the feature-quantification (`.txt`) table required by GNPS
/// Feature-Based Molecular Networking (source `GNPSQuantificationFile`).
///
/// For each consensus feature this writes one row giving the consensus-level
/// `rt_cf`, `mz_cf`, `intensity_cf`, `charge_cf`, `width_cf` and `quality_cf`,
/// optionally followed by the IIMN columns when the first consensus feature
/// carries the `row ID` meta value, followed by the per-map sub-feature
/// columns. The file starts with the `#MAP` and `#CONSENSUS` header lines and
/// one `MAP` row per column header (id, filename, label, size).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GNPSQuantificationFile;

impl GNPSQuantificationFile {
    /// Write the FBMN feature-quantification table of `consensus_map` to
    /// `output_file`, which is overwritten if it exists.
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

    /// Write the FBMN feature-quantification table to `writer`.
    ///
    /// Output structure, tab separated through [`SVOutStream`] without
    /// quoting (a tab inside a text field becomes `_`):
    ///
    /// - `#MAP`, `id`, `filename`, `label`, `size`.
    /// - `#CONSENSUS`, `rt_cf`, `mz_cf`, `intensity_cf`, `charge_cf`,
    ///   `width_cf`, `quality_cf`, then `row ID`, `best ion`, `partners`,
    ///   `annotation network number` if and only if the first consensus
    ///   feature carries `row ID`, then `rt_<i>`, `mz_<i>`, `intensity_<i>`,
    ///   `charge_<i>`, `width_<i>` for `i` from 0 to the number of column
    ///   headers minus one.
    /// - One `MAP` row per column header, in key order: the key, filename,
    ///   label and size.
    /// - One `CONSENSUS` row per consensus feature: its RT, m/z, intensity,
    ///   charge, width and quality; the four IIMN meta values when those
    ///   columns exist (empty when a feature lacks one); then for each `i` the
    ///   RT, m/z, intensity, charge and width of the handle whose map index is
    ///   `i`, or five empty fields when there is none.
    ///
    /// Numbers use the source's `StringUtils::toStr`; a meta value is written
    /// as the source streams a `DataValue` (strings unquoted, lists as
    /// `[a, b]`).
    ///
    /// # Source behaviour kept
    ///
    /// The sub-feature columns are numbered `0..n` for `n` column headers,
    /// whatever the header keys are: when the keys are not exactly `0..n`, a
    /// handle whose map index is `n` or larger is not written. When two handles
    /// share a map index, the later one in handle order (the larger unique ID)
    /// is written.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when a text field (a filename,
    /// label or string meta value) contains a line break, which the source
    /// refuses with `Exception::IllegalArgument`, and [`crate::Error::Io`]
    /// when writing fails.
    pub fn write(consensus_map: &ConsensusMap, writer: impl Write) -> Result<()> {
        let iimn = consensus_map
            .features
            .first()
            .is_some_and(|feature| feature.base.metadata.contains_key(IIMN_ROW_ID));
        let maps = consensus_map.column_headers.len();
        let mut out = SVOutStream::with_options(writer, "\t", "_", QuotingMethod::None)?;

        for field in ["#MAP", "id", "filename", "label", "size"] {
            out.write_field(field)?;
        }
        out.end_line()?;
        for field in [
            "#CONSENSUS",
            "rt_cf",
            "mz_cf",
            "intensity_cf",
            "charge_cf",
            "width_cf",
            "quality_cf",
        ] {
            out.write_field(field)?;
        }
        if iimn {
            for field in IIMN_META_VALUES {
                out.write_field(field)?;
            }
        }
        for i in 0..maps {
            for prefix in ["rt_", "mz_", "intensity_", "charge_", "width_"] {
                out.write_field(&format!("{prefix}{i}"))?;
            }
        }
        out.end_line()?;

        for (id, header) in &consensus_map.column_headers {
            out.write_field("MAP")?;
            out.write_number(*id)?;
            out.write_field(&header.filename)?;
            out.write_field(&header.label)?;
            out.write_number(header.size)?;
            out.end_line()?;
        }

        for feature in &consensus_map.features {
            let base = &feature.base;
            out.write_field("CONSENSUS")?;
            out.write_number(base.rt)?;
            out.write_number(base.mz)?;
            out.write_number(base.intensity)?;
            out.write_number(base.charge)?;
            out.write_number(base.width)?;
            out.write_number(base.quality)?;
            if iimn {
                for key in IIMN_META_VALUES {
                    match base.metadata.get(key) {
                        Some(value) => out.write_display(value.to_list_text()?)?,
                        None => out.write_display("")?,
                    }
                }
            }
            let mut by_map: BTreeMap<u64, &FeatureHandle> = BTreeMap::new();
            for handle in feature.handles() {
                by_map.insert(handle.map_index, handle);
            }
            for i in 0..maps {
                match by_map.get(&(i as u64)) {
                    Some(handle) => {
                        out.write_number(handle.rt)?;
                        out.write_number(handle.mz)?;
                        out.write_number(handle.intensity)?;
                        out.write_number(handle.charge)?;
                        out.write_number(handle.width)?;
                    }
                    None => {
                        for _ in 0..5 {
                            out.write_field("")?;
                        }
                    }
                }
            }
            out.end_line()?;
        }
        out.finish()?;
        Ok(())
    }
}

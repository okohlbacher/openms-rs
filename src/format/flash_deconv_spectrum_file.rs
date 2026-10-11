// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FLASHDeconv spectrum-level output: the `*.tsv` table and the TopFD
//! `*.msalign` file for TopPIC, ported from
//! `FORMAT/FLASHDeconvSpectrumFile.h`.
//!
//! Only the header writers are ported. `writeDeconvolvedMasses`,
//! `writeTopFD`, `writeIsobaricQuantification` and `writeMzML` take
//! `DeconvolvedSpectrum`, `PeakGroup` and
//! `FLASHHelperClasses::PrecalculatedAveragine` records from
//! `ANALYSIS/TOPDOWN`, which have no Rust counterpart yet. See
//! `docs/FLASHDECONV_FILES_SUPPORT.md`.

use crate::Result;
use crate::param::Param;
use std::io::Write;

/// Writer of FLASHDeconv spectrum-level files (source
/// `FLASHDeconvSpectrumFile`, a class of static functions).
///
/// The source's private `topFD_min_peak_count_` (3) and
/// `topFD_max_peak_count_` (500) bound the peaks per mass in `writeTopFD`,
/// which is not ported.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FLASHDeconvSpectrumFile;

/// The columns every spectrum-level header starts with.
const LEAD: &[u8] = b"Index\tFileName\tScanNum\tFeatureIndex\t";
/// The mass columns after the optional decoy column.
const MASS: &[u8] = b"RetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t\
                      SumIntensity\tMinCharge\tMaxCharge\tPeakCount\t";
/// The per-peak detail columns.
const DETAIL: &[u8] = b"PeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\
                        \tPeakPPMErrors\tNoisePeakMZs\tNoisePeakIntensities\tNoisePeakCharges\
                        \tNoisePeakMasses\tNoisePeakIsotopeIndices\tNoisePeakPPMErrors\t";
/// The precursor columns of MS levels above 1.
const PRECURSOR: &[u8] = b"PrecursorScanNum\tPrecursorFeatureIndex\tPrecursorMz\
                           \tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\
                           \tPrecursorMonoisotopicMass\tPrecursorQscore\tPrecursorQscore2D\t";
/// The score columns.
const SCORES: &[u8] = b"IsotopeCosine\tChargeCosine\tChargeScore\tMassSNR\tChargeSNR\
                        \tAveragePPMError\tRepresentativeCharge\tRepresentativeMzStart\
                        \tRepresentativeMzEnd\tQscore\tQscore2D\t";

impl FLASHDeconvSpectrumFile {
    /// Write the header line of the spectrum-level `*.tsv` output (source
    /// `writeDeconvolvedMassesHeader`, `FLASHDeconvSpectrumFile.cpp:236-324`).
    ///
    /// The tab-separated columns are `Index`, `FileName`, `ScanNum`,
    /// `FeatureIndex`, `TargetDecoyType` (with `report_decoy`), the mass
    /// columns `RetentionTime` through `PeakCount`, the peak lists `PeakMZs`
    /// through `NoisePeakPPMErrors` (with `detail`), the precursor columns
    /// `PrecursorScanNum` through `PrecursorQscore2D` and `PrecursorQvalue`
    /// (with `report_decoy`) for an `ms_level` other than 1, the scores
    /// `IsotopeCosine` through `Qscore2D`, `Qvalue` (with `report_decoy`) and,
    /// with `detail`, `PerChargeIntensity` and `PerIsotopeIntensity`.
    ///
    /// The source's spelling is kept exactly, including its trailing tabs:
    /// without `detail`, the line ends in `Qscore2D\t\n` for MS1 without
    /// decoys, in `Qvalue\n` for MS1 with them, and in `Qscore2D\t\n` or
    /// `Qvalue\t\n` above MS1.
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when writing fails; the source writes to an
    /// unchecked stream.
    pub fn write_deconvolved_masses_header(
        mut writer: impl Write,
        ms_level: u32,
        detail: bool,
        report_decoy: bool,
    ) -> Result<()> {
        writer.write_all(LEAD)?;
        if report_decoy {
            writer.write_all(b"TargetDecoyType\t")?;
        }
        writer.write_all(MASS)?;
        if detail {
            writer.write_all(DETAIL)?;
        }
        if ms_level != 1 {
            writer.write_all(PRECURSOR)?;
            if report_decoy {
                writer.write_all(b"PrecursorQvalue\t")?;
            }
        }
        writer.write_all(SCORES)?;
        if report_decoy {
            // The MS1 summary header alone writes no tab after `Qvalue`.
            if ms_level == 1 && !detail {
                writer.write_all(b"Qvalue")?;
            } else {
                writer.write_all(b"Qvalue\t")?;
            }
        }
        if detail {
            writer.write_all(b"PerChargeIntensity\tPerIsotopeIntensity")?;
        }
        writer.write_all(b"\n")?;
        Ok(())
    }

    /// Write the parameter header of a TopFD `*.msalign` file (source
    /// `writeTopFDHeader`, `FLASHDeconvSpectrumFile.cpp:516-525`).
    ///
    /// The header is `#FLASHDeconv generated msalign file`, a
    /// `####################### Parameters ######################` rule, one
    /// `#<name>: <value>` line per parameter in the parameter tree's
    /// traversal order, and the rule again. `<name>` is the entry's own name
    /// without its section path, as the source's `ParamIterator` yields it,
    /// and `<value>` is the source's `ParamValue` stream output
    /// ([`crate::param::ParamValue::to_stream_text`]).
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when writing fails, and a parameter-tree resource
    /// failure.
    pub fn write_top_fd_header(mut writer: impl Write, param: &Param) -> Result<()> {
        const RULE: &[u8] = b"####################### Parameters ######################\n";
        writer.write_all(b"#FLASHDeconv generated msalign file\n")?;
        writer.write_all(RULE)?;
        for item in param.iter()? {
            writeln!(
                writer,
                "#{}: {}",
                item.entry.name,
                item.entry.value.to_stream_text()?
            )?;
        }
        writer.write_all(RULE)?;
        Ok(())
    }
}

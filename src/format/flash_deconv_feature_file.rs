// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FLASHDeconv feature-level output: the `*.tsv` table and the TopFD
//! `*.feature` file for TopPIC, ported from `FORMAT/FLASHDeconvFeatureFile.h`.
//!
//! Only the header writers are ported. The row writers `writeFeatures` and
//! `writeTopFDFeatures` take `FLASHHelperClasses::MassFeature` and
//! `DeconvolvedSpectrum` records from `ANALYSIS/TOPDOWN`, which have no Rust
//! counterpart yet. See `docs/FLASHDECONV_FILES_SUPPORT.md`.

use crate::Result;
use std::io::Write;

/// Writer of FLASHDeconv feature-level files (source
/// `FLASHDeconvFeatureFile`, a class of static functions).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FLASHDeconvFeatureFile;

impl FLASHDeconvFeatureFile {
    /// Write the header line of the regular feature `*.tsv` output (source
    /// `writeHeader`, `FLASHDeconvFeatureFile.cpp:18-27`).
    ///
    /// The tab-separated columns are `FeatureIndex`, `FileName`, `MSLevel`,
    /// then `IsDecoy` when `report_decoy` is set (the source default is
    /// `false`), then `MonoisotopicMass` through `PerIsotopeIntensity`,
    /// followed by a newline.
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when writing fails; the source writes to an
    /// unchecked stream.
    pub fn write_header(mut writer: impl Write, report_decoy: bool) -> Result<()> {
        writer.write_all(b"FeatureIndex\tFileName\tMSLevel")?;
        if report_decoy {
            writer.write_all(b"\tIsDecoy")?;
        }
        writer.write_all(
            b"\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime\
              \tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime\
              \tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\
              \tChargeCount\tIsotopeCosineScore\tQscore2D\tPerChargeIntensity\
              \tPerIsotopeIntensity\n",
        )?;
        Ok(())
    }

    /// Write the header line of the TopFD feature file (source
    /// `writeTopFDFeatureHeader`, `FLASHDeconvFeatureFile.cpp:29-44`).
    ///
    /// MS level 1 writes the MS1 feature columns (`File_name` through
    /// `EC_score`); every other level writes the MS2 columns (`File_name`,
    /// `Fraction_ID`, `Spectrum_ID` through `Precursor_intensity`). The
    /// writer is flushed afterwards, as the source flushes its stream.
    ///
    /// # Errors
    ///
    /// [`crate::Error::Io`] when writing or flushing fails.
    pub fn write_top_fd_feature_header(mut writer: impl Write, ms_level: u32) -> Result<()> {
        if ms_level == 1 {
            writer.write_all(
                b"File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\
                  \tMin_scan\tMax_scan\tMin_charge\tMax_charge\tApex_time\tApex_scan\
                  \tApex_intensity\tRep_charge\tRep_average_mz\tEnvelope_num\tEC_score\n",
            )?;
        } else {
            writer.write_all(
                b"File_name\tFraction_ID\tSpectrum_ID\tScans\tMS_one_ID\tMS_one_scans\
                  \tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_score\
                  \tFraction_feature_min_time\tFraction_feature_max_time\
                  \tFraction_feature_apex_time\tPrecursor_monoisotopic_mz\
                  \tPrecursor_average_mz\tPrecursor_charge\tPrecursor_intensity\n",
            )?;
        }
        writer.flush()?;
        Ok(())
    }
}

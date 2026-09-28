// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Extraction and look-up of spectrum meta data, and the four repair helpers
//! that annotate peptide identifications from raw data.
//!
//! Ports `src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h` and
//! `src/openms/source/METADATA/SpectrumMetaDataLookup.cpp` (core SDK
//! `bc9cc12`). See `docs/SPECTRUM_META_DATA_LOOKUP_SUPPORT.md`.
//!
//! A common use is importing identification results from a search engine's
//! own format, where some meta data - above all retention times - has to be
//! looked up in the raw data; [`SpectrumMetaDataLookup::add_missing_rts_to_peptide_ids`]
//! is one example.
//!
//! Meta data of a spectrum is a [`SpectrumMetaData`] record. Which items to
//! extract or look up is chosen with [`MetaDataFlags`] bits. Meta data comes
//! from spectra or from spectrum reference strings, whose format is a regular
//! expression with named groups:
//!
//! | [`SpectrumMetaData`] field | flag | group | comment (* undefined for MS1) |
//! |---|---|---|---|
//! | `rt` | [`MDF_RT`] | `RT` | retention time of the spectrum |
//! | `precursor_rt` | [`MDF_PRECURSORRT`] | `PRECRT` | retention time of the precursor spectrum* |
//! | `precursor_mz` | [`MDF_PRECURSORMZ`] | `MZ` | m/z of the precursor ion* |
//! | `precursor_charge` | [`MDF_PRECURSORCHARGE`] | `CHARGE` | charge of the precursor ion* |
//! | `ms_level` | [`MDF_MSLEVEL`] | `LEVEL` | MS level |
//! | `scan_number` | [`MDF_SCANNUMBER`] | `SCAN` | scan number, from the native ID |
//! | `native_id` | [`MDF_NATIVEID`] | `ID` | native ID of the spectrum |
//! | | [`MDF_ALL`] | | every flag set |
//! | | | `INDEX0` | look-up only: position counting from 0 |
//! | | | `INDEX1` | look-up only: position counting from 1 |
//!
//! # Example
//!
//! ```
//! use openms::kernel::{MSSpectrum, Precursor};
//! use openms::metadata::spectrum_lookup::DEFAULT_SCAN_REGEXP;
//! use openms::metadata::spectrum_meta_data_lookup::{
//!     MDF_PRECURSORMZ, MDF_RT, SpectrumMetaData, SpectrumMetaDataLookup,
//! };
//!
//! let survey = MSSpectrum { native_id: "spectrum=0".into(), rt: 1.0, ms_level: 1, ..MSSpectrum::default() };
//! let fragment = MSSpectrum {
//!     native_id: "spectrum=1".into(),
//!     rt: 2.0,
//!     ms_level: 2,
//!     precursors: vec![Precursor::new(1000.0, 2)],
//!     ..MSSpectrum::default()
//! };
//! let mut lookup = SpectrumMetaDataLookup::new();
//! lookup.read_spectra(&[survey, fragment], DEFAULT_SCAN_REGEXP, true)?;
//! assert_eq!(lookup.spectrum_meta_data(1)?.precursor_rt, 1.0);
//!
//! lookup.add_reference_format(r"rt=(?<RT>\d+(\.\d+)?),mz=(?<MZ>\d+(\.\d+)?)")?;
//! let mut meta = SpectrumMetaData::default();
//! lookup.spectrum_meta_data_from_reference("rt=5.0,mz=1000.0", &mut meta, MDF_RT | MDF_PRECURSORMZ)?;
//! assert_eq!((meta.rt, meta.precursor_mz), (5.0, 1000.0)); // no look-up needed
//! # Ok::<(), openms::Error>(())
//! ```

use std::collections::BTreeMap;
use std::ops::Deref;
use std::path::Path;

use crate::concept::boost_regex::BoostRegex;
use crate::concept::constants::user_param::{FAIMS_CV, IM};
use crate::concept::string_conversion::{to_double, to_int32};
use crate::format::{FileHandler, FileType, PeakFileOptions};
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::faims_helper::FaimsHelper;
use crate::kernel::{MSExperiment, MSSpectrum};
use crate::metadata::spectrum_lookup::SpectrumLookup;
use crate::metadata::spectrum_native_id_parser;
use crate::metadata::{DriftTimeUnit, ImTypes, IonMobilityFormat, MetaValue};
use crate::{Error, Result};

/// Bit set choosing which meta data to extract or look up
/// (`SpectrumMetaDataLookup::MetaDataFlags`, an `unsigned char`).
pub type MetaDataFlags = u8;

/// Retention time (`MDF_RT`).
pub const MDF_RT: MetaDataFlags = 1;
/// Precursor retention time (`MDF_PRECURSORRT`).
pub const MDF_PRECURSORRT: MetaDataFlags = 2;
/// Precursor m/z (`MDF_PRECURSORMZ`).
pub const MDF_PRECURSORMZ: MetaDataFlags = 4;
/// Precursor charge (`MDF_PRECURSORCHARGE`).
pub const MDF_PRECURSORCHARGE: MetaDataFlags = 8;
/// MS level (`MDF_MSLEVEL`).
pub const MDF_MSLEVEL: MetaDataFlags = 16;
/// Scan number (`MDF_SCANNUMBER`).
pub const MDF_SCANNUMBER: MetaDataFlags = 32;
/// Native ID (`MDF_NATIVEID`).
pub const MDF_NATIVEID: MetaDataFlags = 64;
/// Every flag (`MDF_ALL`).
pub const MDF_ALL: MetaDataFlags = 127;

/// Meta data of one spectrum (`SpectrumMetaDataLookup::SpectrumMetaData`).
///
/// [`Default`] is the source's constructor: every floating-point item quiet
/// NaN (not known), charge and MS level 0, scan number -1, an empty native ID.
/// The source's sentinels are kept because the reference extraction fills
/// items one by one into a record the caller owns.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumMetaData {
    /// Retention time.
    pub rt: f64,
    /// Retention time of the precursor spectrum.
    pub precursor_rt: f64,
    /// Precursor mass-to-charge ratio.
    pub precursor_mz: f64,
    /// Precursor charge.
    pub precursor_charge: i32,
    /// MS level.
    pub ms_level: usize,
    /// Scan number, -1 when not known.
    pub scan_number: i32,
    /// Native ID.
    pub native_id: String,
}

impl Default for SpectrumMetaData {
    fn default() -> Self {
        Self {
            rt: f64::NAN,
            precursor_rt: f64::NAN,
            precursor_mz: f64::NAN,
            precursor_charge: 0,
            ms_level: 0,
            scan_number: -1,
            native_id: String::new(),
        }
    }
}

/// What a repair helper did: the source's `bool` result and the lines it
/// logs, in order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Annotation {
    /// The source's return value.
    pub success: bool,
    /// The messages the source writes to its log streams, each with the
    /// source's own prefix where it has one (`Error: `).
    pub log: Vec<String>,
}

/// Helper class for looking up spectrum meta data (`SpectrumMetaDataLookup`).
///
/// The source derives from `SpectrumLookup`; here the lookup is composed and
/// reachable through [`Deref`], so every `find_by_*` method, `is_empty`,
/// `reference_formats` and `rt_tolerance` read as on the base. Mutation goes
/// through [`add_reference_format`](Self::add_reference_format) and
/// [`set_rt_tolerance`](Self::set_rt_tolerance) only: the base's own
/// `read_spectra` would re-index the spectra without the meta data this type
/// keeps beside the index, which the source's name hiding merely discourages.
#[derive(Clone, Debug, Default)]
pub struct SpectrumMetaDataLookup {
    lookup: SpectrumLookup,
    metadata: Vec<SpectrumMetaData>,
    spectra_data_ref: String,
}

impl Deref for SpectrumMetaDataLookup {
    type Target = SpectrumLookup;

    fn deref(&self) -> &SpectrumLookup {
        &self.lookup
    }
}

impl SpectrumMetaDataLookup {
    /// An empty lookup, the source's constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a spectrum reference format, as
    /// [`SpectrumLookup::add_reference_format`].
    ///
    /// # Errors
    ///
    /// As [`SpectrumLookup::add_reference_format`].
    pub fn add_reference_format(&mut self, regexp: &str) -> Result<()> {
        self.lookup.add_reference_format(regexp)
    }

    /// Set the retention-time tolerance, the source's public `rt_tolerance`.
    pub fn set_rt_tolerance(&mut self, rt_tolerance: f64) {
        self.lookup.rt_tolerance = rt_tolerance;
    }

    /// The spectra-data reference recorded by
    /// [`set_spectra_data_ref`](Self::set_spectra_data_ref).
    pub fn spectra_data_ref(&self) -> &str {
        &self.spectra_data_ref
    }

    /// Record the file the spectra came from (`setSpectraDataRef`).
    pub fn set_spectra_data_ref(&mut self, spectra_data: impl Into<String>) {
        self.spectra_data_ref = spectra_data.into();
    }

    /// Read `spectra`, extract their meta data and index them for look-up
    /// (`readSpectra`); returns the errors the source logs.
    ///
    /// Scan numbers are extracted with `scan_regexp`
    /// ([`DEFAULT_SCAN_REGEXP`](crate::metadata::spectrum_lookup::DEFAULT_SCAN_REGEXP)
    /// is the source's default argument); an empty expression disables them.
    /// With `get_precursor_rt` (the source's default is `false`), the
    /// precursor retention time of each spectrum with a precursor is the
    /// retention time of the last preceding spectrum one MS level lower.
    ///
    /// Unlike the source, a second call replaces what the first read: the
    /// source's `readSpectra` clears neither its meta data nor the base's
    /// indexes, so after re-reading, index *i* still returns the first read's
    /// meta data and stale retention times, IDs and scan numbers keep
    /// resolving (recorded as a source defect candidate). For a lookup read
    /// once, the result is the source's.
    ///
    /// # Errors
    ///
    /// As [`SpectrumLookup::read_spectra`], which builds the index; nothing is
    /// changed on error.
    pub fn read_spectra(
        &mut self,
        spectra: &[MSSpectrum],
        scan_regexp: &str,
        get_precursor_rt: bool,
    ) -> Result<Vec<String>> {
        let mut lookup = self.lookup.clone();
        // The index's own warnings are the base class's wording; this class
        // logs its own error per spectrum instead, as the source does.
        lookup.read_spectra(spectra, scan_regexp)?;
        let regex = if scan_regexp.is_empty() {
            None
        } else {
            Some(BoostRegex::new(scan_regexp)?)
        };
        let mut log = Vec::new();
        let mut metadata = Vec::with_capacity(spectra.len());
        let mut precursor_rts = BTreeMap::new();
        for spectrum in spectra {
            let mut meta = SpectrumMetaData::default();
            log.extend(Self::spectrum_meta_data_of(
                spectrum,
                &mut meta,
                regex.as_ref(),
                &precursor_rts,
            )?);
            if get_precursor_rt {
                precursor_rts.insert(meta.ms_level, meta.rt);
            }
            metadata.push(meta);
        }
        self.lookup = lookup;
        self.metadata = metadata;
        Ok(log)
    }

    /// The meta data of the spectrum at `index`
    /// (`getSpectrumMetaData(Size, SpectrumMetaData&)`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `index` is not below the number of spectra
    /// read, where the source throws `Exception::IndexOverflow`.
    pub fn spectrum_meta_data(&self, index: usize) -> Result<SpectrumMetaData> {
        self.metadata.get(index).cloned().ok_or_else(|| {
            Error::InvalidRange(format!(
                "index {index} is not below the number of spectra ({})",
                self.metadata.len()
            ))
        })
    }

    /// Extract meta data from a spectrum into `meta`
    /// (static `getSpectrumMetaData(const MSSpectrum&, ...)`); returns the
    /// errors the source logs.
    ///
    /// Native ID, retention time and MS level are always written. The scan
    /// number is written only when `scan_regexp` is given (the source's
    /// empty `boost::regex` default is `None`): the extracted number, or -1
    /// with the source's error line when extraction fails. For a spectrum
    /// with a precursor, the first precursor's m/z and charge are written,
    /// and, when `precursor_rts` is not empty, the entry for the MS level one
    /// below the spectrum's (with an error line when there is none). Other
    /// fields of `meta` are left as they were.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the regular-expression search exceeds its
    /// work budget; the source has no budget.
    pub fn spectrum_meta_data_of(
        spectrum: &MSSpectrum,
        meta: &mut SpectrumMetaData,
        scan_regexp: Option<&BoostRegex>,
        precursor_rts: &BTreeMap<usize, f64>,
    ) -> Result<Vec<String>> {
        let mut log = Vec::new();
        meta.native_id.clone_from(&spectrum.native_id);
        meta.rt = spectrum.rt;
        meta.ms_level = spectrum.ms_level as usize;
        if let Some(regex) = scan_regexp {
            let scan =
                spectrum_native_id_parser::try_extract_scan_number_with(&meta.native_id, regex)?
                    .unwrap_or(-1);
            meta.scan_number = scan;
            if scan < 0 {
                log.push(format!(
                    "Error: Could not extract scan number from spectrum native ID '{}' using regular expression '{}'.",
                    meta.native_id,
                    regex.as_str()
                ));
            }
        }
        if let Some(precursor) = spectrum.precursors.first() {
            meta.precursor_mz = precursor.mz;
            meta.precursor_charge = precursor.charge;
            if !precursor_rts.is_empty() {
                // The source's `ms_level - 1` wraps for level 0 to a key no
                // map holds; checked subtraction finds nothing there either.
                match meta
                    .ms_level
                    .checked_sub(1)
                    .and_then(|level| precursor_rts.get(&level))
                {
                    Some(rt) => meta.precursor_rt = *rt,
                    None => log.push(format!(
                        "Error: Could not set precursor RT for spectrum with native ID '{}' - precursor spectrum not found.",
                        meta.native_id
                    )),
                }
            }
        }
        Ok(log)
    }

    /// Meta data for a spectrum reference, extracted from the reference where
    /// possible and looked up otherwise
    /// (`getSpectrumMetaData(const std::string&, SpectrumMetaData&, MetaDataFlags)`).
    ///
    /// The reference formats are searched for in `spectrum_ref` in order and
    /// the first that matches decides; when none matches, `meta` is left
    /// unchanged and this returns `Ok`, as the source does. From the match,
    /// each item requested in `flags` whose group captured a non-empty value
    /// is written to `meta` and its flag cleared (a native ID group is written
    /// even when empty, and then clears nothing). If any requested flag is
    /// left, the spectrum is looked up as in
    /// [`SpectrumLookup::find_by_reference`] and `meta` is **replaced** by its
    /// stored meta data, discarding what was extracted, as in the source.
    ///
    /// Numbers convert as `StringUtils::toDouble` (retention times, m/z and
    /// the charge, which the source then truncates to `Int`) and `toInt32`
    /// (MS level, scan number).
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when a captured number does not convert (the
    ///   source's `ConversionError`), when a charge lies outside `i32` (the
    ///   source's conversion is undefined there), or when an MS level is
    ///   negative (the source wraps it to a huge `Size`);
    /// - the look-up's errors, as [`SpectrumLookup::find_by_reference`]
    ///   (`ElementNotFound` becomes [`Error::InvalidRange`]);
    /// - [`Error::InvalidValue`] when a search exceeds its work budget.
    pub fn spectrum_meta_data_from_reference(
        &self,
        spectrum_ref: &str,
        meta: &mut SpectrumMetaData,
        flags: MetaDataFlags,
    ) -> Result<()> {
        let mut flags = flags;
        for regex in &self.lookup.reference_formats {
            let Some(captures) = regex.search(spectrum_ref.as_bytes())? else {
                continue;
            };
            let group = |name: &str| {
                captures
                    .name(name)
                    .and_then(|range| spectrum_ref.get(range))
            };
            let wanted = |flags: MetaDataFlags, flag: MetaDataFlags| flags & flag == flag;
            let double = |value: &str| to_double(value).map_err(Error::InvalidValue);
            let int = |value: &str| to_int32(value).map_err(Error::InvalidValue);
            if wanted(flags, MDF_RT) {
                if let Some(value) = group("RT").filter(|v| !v.is_empty()) {
                    meta.rt = double(value)?;
                    flags &= !MDF_RT;
                }
            }
            if wanted(flags, MDF_PRECURSORRT) {
                if let Some(value) = group("PRECRT").filter(|v| !v.is_empty()) {
                    meta.precursor_rt = double(value)?;
                    flags &= !MDF_PRECURSORRT;
                }
            }
            if wanted(flags, MDF_PRECURSORMZ) {
                if let Some(value) = group("MZ").filter(|v| !v.is_empty()) {
                    meta.precursor_mz = double(value)?;
                    flags &= !MDF_PRECURSORMZ;
                }
            }
            if wanted(flags, MDF_PRECURSORCHARGE) {
                if let Some(value) = group("CHARGE").filter(|v| !v.is_empty()) {
                    meta.precursor_charge = truncate_charge(double(value)?, value)?;
                    flags &= !MDF_PRECURSORCHARGE;
                }
            }
            if wanted(flags, MDF_MSLEVEL) {
                if let Some(value) = group("LEVEL").filter(|v| !v.is_empty()) {
                    let level = int(value)?;
                    meta.ms_level = usize::try_from(level).map_err(|_| {
                        Error::InvalidValue(format!("MS level '{value}' is negative"))
                    })?;
                    flags &= !MDF_MSLEVEL;
                }
            }
            if wanted(flags, MDF_SCANNUMBER) {
                if let Some(value) = group("SCAN").filter(|v| !v.is_empty()) {
                    meta.scan_number = int(value)?;
                    flags &= !MDF_SCANNUMBER;
                }
            }
            if wanted(flags, MDF_NATIVEID) {
                if let Some(value) = group("ID") {
                    meta.native_id = value.to_owned();
                    if !meta.native_id.is_empty() {
                        flags &= !MDF_NATIVEID;
                    }
                }
            }
            if flags != 0 {
                let index =
                    self.lookup
                        .find_by_regexp_match(spectrum_ref, regex.as_str(), &captures)?;
                *meta = self.spectrum_meta_data(index)?;
            }
            return Ok(());
        }
        Ok(())
    }

    /// Fill in missing retention times of peptide identifications from the
    /// spectra of `experiment`, matched by native ID
    /// (`addMissingRTsToPeptideIDs`).
    ///
    /// A retention time is missing when it is `None` or NaN (the source's
    /// quiet-NaN default). `success` is `false` for an experiment without
    /// spectra (nothing is changed, with the source's info line) and when a
    /// spectrum reference names no spectrum (that identification keeps its
    /// missing time, with the source's error line).
    ///
    /// # Errors
    ///
    /// As [`SpectrumLookup::read_spectra`] with the default scan expression.
    pub fn add_missing_rts_to_peptide_ids(
        peptides: &mut [PeptideIdentification],
        experiment: &MSExperiment,
    ) -> Result<Annotation> {
        let mut log = Vec::new();
        if experiment.spectra.is_empty() {
            log.push("No spectra found in the experiment. Skipping RT annotation.".to_owned());
            return Ok(Annotation {
                success: false,
                log,
            });
        }
        let lookup = default_lookup(&experiment.spectra)?;
        let mut success = true;
        for peptide in peptides.iter_mut() {
            if !peptide.rt.is_none_or(f64::is_nan) {
                continue;
            }
            let native_id = peptide.spectrum_reference();
            match lookup.find_by_native_id(&native_id) {
                Ok(index) => peptide.rt = Some(experiment.spectra[index].rt),
                Err(_) => {
                    log.push(format!(
                        "Error: Failed to look up retention time for peptide identification with spectrum reference '{native_id}' - no spectrum with corresponding native ID found."
                    ));
                    success = false;
                }
            }
        }
        Ok(Annotation { success, log })
    }

    /// Annotate every peptide identification with the drift time of its
    /// spectrum as the `IM` meta value (`addMissingIMToPeptideIDs`).
    ///
    /// Only spectra whose ion-mobility format is one value per spectrum
    /// ([`IonMobilityFormat::PerSpectrum`], the source's `IM_SPECTRUM`) give
    /// a value; `success` is `false` when any identification's spectrum is
    /// of another format, and for an experiment without spectra (with the
    /// source's info line). The format is determined from the data
    /// ([`ImTypes::determine_im_format`]); the source first consults a
    /// stored format, which this crate's spectrum does not carry.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when a spectrum reference names no spectrum:
    /// the source does not catch that `ElementNotFound` here, so it leaves the
    /// function with the identifications before it already annotated, as
    /// this does. As [`SpectrumLookup::read_spectra`] otherwise, and
    /// [`Error::InvalidValue`] for a non-finite drift time.
    pub fn add_missing_im_to_peptide_ids(
        peptides: &mut [PeptideIdentification],
        experiment: &MSExperiment,
    ) -> Result<Annotation> {
        let mut log = Vec::new();
        if experiment.spectra.is_empty() {
            log.push("No spectra found in the experiment. Skipping IM annotation.".to_owned());
            return Ok(Annotation {
                success: false,
                log,
            });
        }
        let lookup = default_lookup(&experiment.spectra)?;
        let mut all_ids_have_im = true;
        for peptide in peptides.iter_mut() {
            let index = lookup.find_by_native_id(&peptide.spectrum_reference())?;
            let spectrum = &experiment.spectra[index];
            if ImTypes::determine_im_format(spectrum) == IonMobilityFormat::PerSpectrum {
                peptide
                    .metadata
                    .insert(IM.into(), MetaValue::try_from(spectrum.drift_time)?);
            } else {
                all_ids_have_im = false;
            }
        }
        Ok(Annotation {
            success: all_ids_have_im,
            log,
        })
    }

    /// Annotate peptide identifications with the FAIMS compensation voltage
    /// of their spectrum as the `FAIMS_CV` meta value
    /// (`addMissingFAIMSToPeptideIDs`).
    ///
    /// A spectrum whose drift-time unit is the FAIMS compensation voltage
    /// gives its drift time; an MS2+ spectrum without one inherits the last
    /// voltage seen before it. Identifications that already carry `FAIMS_CV`
    /// are counted and left alone; a reference naming no spectrum gets the
    /// source's warning. `success` is whether any identification is
    /// annotated, and `false` without spectra or without FAIMS data. The log
    /// ends with the source's summary line.
    ///
    /// # Errors
    ///
    /// As [`FaimsHelper::get_compensation_voltages`] (a NaN voltage) and
    /// [`SpectrumLookup::read_spectra`]; [`Error::InvalidValue`] for a
    /// non-finite voltage stored as a meta value.
    pub fn add_missing_faims_to_peptide_ids(
        peptides: &mut [PeptideIdentification],
        experiment: &MSExperiment,
    ) -> Result<Annotation> {
        let mut log = Vec::new();
        if experiment.spectra.is_empty() {
            log.push("No spectra found in the experiment. Skipping FAIMS annotation.".to_owned());
            return Ok(Annotation {
                success: false,
                log,
            });
        }
        if FaimsHelper::get_compensation_voltages(experiment)?
            .voltages
            .is_empty()
        {
            log.push(
                "No FAIMS data found in the experiment. Skipping FAIMS annotation.".to_owned(),
            );
            return Ok(Annotation {
                success: false,
                log,
            });
        }
        let lookup = default_lookup(&experiment.spectra)?;
        let mut index_to_faims_cv = BTreeMap::new();
        let mut last_faims_cv = f64::NAN;
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            if spectrum.drift_time_unit == DriftTimeUnit::FaimsCompensationVoltage {
                last_faims_cv = spectrum.drift_time;
                index_to_faims_cv.insert(index, last_faims_cv);
            } else if !last_faims_cv.is_nan() && spectrum.ms_level > 1 {
                index_to_faims_cv.insert(index, last_faims_cv);
            }
        }
        let mut annotated_count = 0_usize;
        for peptide in peptides.iter_mut() {
            if peptide.metadata.contains_key(FAIMS_CV) {
                annotated_count += 1;
                continue;
            }
            let native_id = peptide.spectrum_reference();
            match lookup.find_by_native_id(&native_id) {
                Ok(index) => {
                    if let Some(cv) = index_to_faims_cv.get(&index) {
                        peptide
                            .metadata
                            .insert(FAIMS_CV.into(), MetaValue::try_from(*cv)?);
                        annotated_count += 1;
                    }
                }
                Err(_) => log.push(format!(
                    "Could not find spectrum with native ID '{native_id}' for FAIMS annotation."
                )),
            }
        }
        log.push(format!(
            "Annotated {annotated_count} of {} peptide IDs with FAIMS CV values.",
            peptides.len()
        ));
        Ok(Annotation {
            success: annotated_count > 0,
            log,
        })
    }

    /// Set the spectrum reference of peptide identifications to the native ID
    /// of the spectrum at their retention time in the raw-data file
    /// `filename` (`addMissingSpectrumReferences`).
    ///
    /// The file (mzXML, mzML, mzData or MGF) is loaded without peak data, as
    /// the source does. Identifications that already have a
    /// `spectrum_reference` are skipped unless `override_spectra_references`.
    /// A retention time matching no spectrum within the default tolerance
    /// logs the source's error line and makes `success` `false`; with
    /// `stop_on_error`, processing stops there. With `override_spectra_data`,
    /// every protein identification's primary MS run path (its
    /// `spectra_data`) becomes `["file://" + filename]`.
    ///
    /// Divergence: the source takes its protein identifications **by value**
    /// and sets `spectra_data` on that copy, so the caller never sees it
    /// (recorded as a source defect candidate); this takes them by mutable
    /// reference so that `override_spectra_data` has the documented effect.
    /// Pass an empty slice for the source's default argument.
    ///
    /// # Errors
    ///
    /// The loader's error when the file is missing or unreadable, where the
    /// source throws `Exception::FileNotFound` before anything is changed, as
    /// here; its errors for a type outside the four; and as
    /// [`SpectrumMetaDataLookup::read_spectra`].
    pub fn add_missing_spectrum_references(
        peptides: &mut [PeptideIdentification],
        filename: &str,
        stop_on_error: bool,
        override_spectra_data: bool,
        override_spectra_references: bool,
        proteins: &mut [ProteinIdentification],
    ) -> Result<Annotation> {
        let mut options = PeakFileOptions::default();
        options.fill_data = false;
        options.skip_xml_checks = true;
        let experiment = FileHandler::load_experiment_with_options(
            Path::new(filename),
            &[
                FileType::MzXml,
                FileType::MzMl,
                FileType::MzData,
                FileType::Mgf,
            ],
            &options,
        )?;
        let mut lookup = Self::new();
        lookup.read_spectra(
            &experiment.spectra,
            crate::metadata::spectrum_lookup::DEFAULT_SCAN_REGEXP,
            false,
        )?;
        lookup.set_spectra_data_ref(filename);
        if override_spectra_data {
            let spectra_data = vec![format!("file://{}", lookup.spectra_data_ref)];
            for protein in proteins.iter_mut() {
                protein.primary_ms_run_paths.clone_from(&spectra_data);
            }
        }
        let mut log = Vec::new();
        let mut success = true;
        for peptide in peptides.iter_mut() {
            if !override_spectra_references && peptide.metadata.contains_key("spectrum_reference") {
                continue;
            }
            let rt = peptide.rt.unwrap_or(f64::NAN);
            match lookup.find_by_rt(rt) {
                Ok(index) => {
                    let meta = lookup.spectrum_meta_data(index)?;
                    peptide.set_spectrum_reference(meta.native_id);
                }
                Err(_) => {
                    log.push(format!(
                        "Error: Failed to look up spectrum native ID for peptide identification with retention time '{}'.",
                        crate::param::value::format_float(rt, true)
                    ));
                    success = false;
                    if stop_on_error {
                        break;
                    }
                }
            }
        }
        Ok(Annotation { success, log })
    }
}

/// A `SpectrumLookup` over `spectra` with the default scan expression, as the
/// repair helpers build one.
fn default_lookup(spectra: &[MSSpectrum]) -> Result<SpectrumLookup> {
    let mut lookup = SpectrumLookup::new();
    lookup.read_spectra(
        spectra,
        crate::metadata::spectrum_lookup::DEFAULT_SCAN_REGEXP,
    )?;
    Ok(lookup)
}

/// The source's implicit `double` to `Int` conversion of a charge:
/// truncation toward zero, refused where it is undefined.
fn truncate_charge(value: f64, text: &str) -> Result<i32> {
    let truncated = value.trunc();
    if truncated.is_finite() && truncated >= f64::from(i32::MIN) && truncated <= f64::from(i32::MAX)
    {
        // In range and integral, so the cast is exact.
        Ok(truncated as i32)
    } else {
        Err(Error::InvalidValue(format!(
            "precursor charge '{text}' is outside the range of a 32-bit integer"
        )))
    }
}

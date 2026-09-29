// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Quality control for OpenSWATH: the MS1 charge-state distribution.
//!
//! Ports `ANALYSIS/OPENSWATH/SwathQC.h` with `SwathQC.cpp` (namespace
//! `OpenSwath`). A [`SwathQC`] collects information while spectra stream past,
//! as a plugin of a transforming data consumer; when the data are already in
//! memory, [`SwathQC::charge_distribution_of_maps`] does the same in one call.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! The source hands out two lambdas that capture `this`
//! (`getExpSettingsFunc`, `getSpectraProcessingFunc`) for a consumer's
//! callbacks. Their bodies are the methods
//! [`SwathQC::apply_experimental_settings`] and [`SwathQC::process_spectrum`];
//! [`SwathQC::exp_settings_func`] and [`SwathQC::spectra_processing_func`]
//! return closures over them, which borrow the collector mutably for as long
//! as they live.

use super::data_access_helper::convert_to_openms_spectrum;
use super::swath_map::SwathMap;
use crate::comparison::Tolerance;
use crate::kernel::{MSSpectrum, SpectrumType};
use crate::metadata::{ExperimentalSettings, MetaValueData};
use crate::processing::deisotoping::Deisotoper;
use crate::processing::peak_picking::PeakPickerHiRes;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Charge state to number of isotopic envelopes, source `ChargeDistribution`
/// (`std::map<int,int>`, ordered by charge).
pub type ChargeDistribution = BTreeMap<i32, i32>;

/// Collector of OpenSWATH quality-control data.
///
/// Source `OpenSwath::SwathQC`; its deleted default constructor has no
/// counterpart, so [`SwathQC::new`] is the only way to build one.
#[derive(Clone, Debug, PartialEq)]
pub struct SwathQC {
    /// The charge distribution collected so far.
    cd: ChargeDistribution,
    /// Number of MS1 spectra expected (0 = unknown, sample everything).
    nr_ms1_spectra: usize,
    /// Number of MS1 spectra to inspect for the charge distribution.
    cd_spectra: usize,
    /// m/z tolerance (Th) for isotope deconvolution.
    decon_ms1_mz_tol: f64,
    /// Number of spectra taken by [`SwathQC::process_spectrum`] so far.
    ms1_spectra_seen: usize,
}

impl SwathQC {
    /// A collector that inspects `cd_spectra` MS1 spectra for the charge
    /// distribution, deconvolving isotopes within `decon_ms1_mz_tol` Th.
    ///
    /// Source constructor `SwathQC(cd_spectra, decon_ms1_mz_tol)`. The number
    /// of expected MS1 spectra starts at 0 (unknown).
    pub fn new(cd_spectra: usize, decon_ms1_mz_tol: f64) -> Self {
        Self {
            cd: ChargeDistribution::new(),
            nr_ms1_spectra: 0,
            cd_spectra,
            decon_ms1_mz_tol,
            ms1_spectra_seen: 0,
        }
    }

    /// Take the number of expected MS1 spectra from the meta value
    /// `nr_ms1_spectra` of `es`, the body of the source's
    /// `getExpSettingsFunc` lambda.
    ///
    /// When a number is already set this does nothing (some parsers call the
    /// callback during parsing, when the information may still be missing).
    /// Without the meta value the number stays 0.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `nr_ms1_spectra` is not an integer (the
    /// source's `ConversionError`), or is negative or beyond `i32` (the
    /// source truncates to `int` and wraps to `size_t`); the collector is then
    /// unchanged.
    pub fn apply_experimental_settings(&mut self, es: &ExperimentalSettings) -> Result<()> {
        if self.nr_ms1_spectra > 0 {
            return Ok(());
        }
        match es.metadata.get("nr_ms1_spectra") {
            None => self.nr_ms1_spectra = 0,
            Some(value) => match value.data() {
                MetaValueData::Integer(n) => {
                    let n = i32::try_from(*n)
                        .ok()
                        .and_then(|n| usize::try_from(n).ok())
                        .ok_or_else(|| {
                            Error::InvalidValue(format!(
                                "meta value 'nr_ms1_spectra' = {n} is not a count"
                            ))
                        })?;
                    self.nr_ms1_spectra = n;
                }
                _ => {
                    return Err(Error::InvalidValue(
                        "meta value 'nr_ms1_spectra' is not an integer".into(),
                    ));
                }
            },
        }
        Ok(())
    }

    /// A callback over [`SwathQC::apply_experimental_settings`], the source's
    /// `getExpSettingsFunc`.
    pub fn exp_settings_func(&mut self) -> impl FnMut(&ExperimentalSettings) -> Result<()> + '_ {
        move |es| self.apply_experimental_settings(es)
    }

    /// Count the charge states of one MS1 spectrum into the distribution, the
    /// body of the source's `getSpectraProcessingFunc` lambda.
    ///
    /// Only MS1 spectra are looked at, and of those only the ones picked by
    /// [`SwathQC::is_subsampled_spectrum`] from the expected MS1 count (set by
    /// [`SwathQC::set_nr_ms1_spectra`] or the experimental settings; with 0
    /// every spectrum is taken). A profile spectrum is first centroided by
    /// [`PeakPickerHiRes`] with default parameters; a centroided one is used
    /// as is; a spectrum of unknown type (after estimating it from the data)
    /// is skipped. The isotopic envelopes are then found by
    /// [`Deisotoper`] (tolerance `decon_ms1_mz_tol` Th, charges 1 to 10, 3 to
    /// 10 isotope peaks, keeping only deisotoped peaks, annotating the charge,
    /// sharing isotopes between envelopes as the source does), and each kept
    /// peak's charge is counted. This picks up non-peptide signals too.
    ///
    /// # Errors
    ///
    /// The errors of type estimation, peak picking and deisotoping, and
    /// [`Error::MissingInformation`] when the deisotoped spectrum carries no
    /// trailing `charge` integer array (the source's
    /// `Exception::Postcondition`). The distribution is unchanged on error,
    /// but the spectrum still counts as seen, as in the source.
    pub fn process_spectrum(&mut self, spec: &MSSpectrum) -> Result<()> {
        if spec.ms_level != 1 {
            return Ok(());
        }
        if !Self::is_subsampled_spectrum(
            self.nr_ms1_spectra,
            self.cd_spectra,
            self.ms1_spectra_seen,
        ) {
            return Ok(());
        }
        self.ms1_spectra_seen += 1;

        let tmp = match spec.get_type(true)? {
            SpectrumType::Profile => PeakPickerHiRes::default().pick_spectrum(spec)?.spectrum,
            SpectrumType::Centroid => spec.clone(),
            // Unknown: too dangerous to analyse.
            _ => return Ok(()),
        };
        if tmp.is_empty() {
            return Ok(());
        }
        let deisotoper = Deisotoper {
            tolerance: Tolerance::Absolute(self.decon_ms1_mz_tol),
            min_charge: 1,
            max_charge: 10,
            keep_only_deisotoped: true,
            min_isotope_peaks: 3,
            max_isotope_peaks: 10,
            make_single_charged: false,
            annotate_charge: true,
            allow_shared_isotopes: true,
            ..Deisotoper::default()
        };
        let result = deisotoper.deisotope(&tmp)?.spectrum;
        let ida = result.integer_data_arrays.last().ok_or_else(|| {
            Error::MissingInformation("IntegerDataArray must not be empty!".into())
        })?;
        if ida.name != "charge" {
            return Err(Error::MissingInformation(
                "IntegerDataArray.back().getName() != \"charge\"".into(),
            ));
        }
        let mut cd = self.cd.clone();
        for &q in &ida.data {
            let count = cd.entry(q).or_insert(0);
            *count = count.checked_add(1).ok_or_else(|| {
                Error::InvalidRange("charge count exceeds the source's int range".into())
            })?;
        }
        self.cd = cd;
        Ok(())
    }

    /// A callback over [`SwathQC::process_spectrum`], the source's
    /// `getSpectraProcessingFunc`.
    pub fn spectra_processing_func(&mut self) -> impl FnMut(&MSSpectrum) -> Result<()> + '_ {
        move |spec| self.process_spectrum(spec)
    }

    /// Sample the spectra of every MS1 SWATH map and count their charge
    /// states.
    ///
    /// Source static `getChargeDistribution(swath_maps, nr_samples, mz_tol)`.
    /// From each map with `ms1` set, `nr_samples` spectra are taken evenly
    /// (by [`SwathQC::is_subsampled_spectrum`]; `usize::MAX`, the source's
    /// `-1`, takes all), converted to kernel spectra (peaks only, MS level 1,
    /// type estimated from the data) and passed through
    /// [`SwathQC::process_spectrum`] of a collector that samples everything.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when an MS1 map has no spectrum access
    /// (the source dereferences a null pointer), and the errors of spectrum
    /// access, conversion and [`SwathQC::process_spectrum`].
    pub fn charge_distribution_of_maps(
        swath_maps: &[SwathMap],
        nr_samples: usize,
        mz_tol: f64,
    ) -> Result<ChargeDistribution> {
        let mut sq = SwathQC::new(nr_samples, mz_tol);
        sq.set_nr_ms1_spectra(0);
        for map in swath_maps.iter().filter(|m| m.ms1) {
            let sptr = map.sptr.as_ref().ok_or_else(|| {
                Error::MissingInformation("an MS1 SWATH map has no spectrum access".into())
            })?;
            let nr_spec = sptr.nr_spectra();
            // One spectrum object is reused, as in the source: only its peaks
            // are replaced, so it keeps MS level 1 and an unknown type.
            let mut s = MSSpectrum::default();
            for i in 0..nr_spec {
                if !Self::is_subsampled_spectrum(nr_spec, nr_samples, i) {
                    continue;
                }
                convert_to_openms_spectrum(&sptr.spectrum_by_id(i)?, &mut s)?;
                sq.process_spectrum(&s)?;
            }
        }
        Ok(sq.cd)
    }

    /// The collected distribution as the source's JSON document: an object
    /// with one key `ChargeDistributionMS1` whose value is an array of
    /// `[charge, count]` pairs in ascending charge order, pretty-printed with
    /// two-space indentation and a trailing newline (nlohmann's
    /// `std::setw(2) << json << std::endl`).
    pub fn to_json(&self) -> String {
        let mut out = String::from("{\n  \"ChargeDistributionMS1\": [");
        if self.cd.is_empty() {
            out.push_str("]\n}\n");
            return out;
        }
        let pairs: Vec<String> = self
            .cd
            .iter()
            .map(|(charge, count)| format!("\n    [\n      {charge},\n      {count}\n    ]"))
            .collect();
        out.push_str(&pairs.join(","));
        out.push_str("\n  ]\n}\n");
        out
    }

    /// Save the collected data as JSON (see [`SwathQC::to_json`]) to `path`.
    ///
    /// Source `storeJSON`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be written (the source's
    /// `Exception::FileNotWritable`, checked after writing so that a full disk
    /// is caught too).
    pub fn store_json(&self, path: impl AsRef<Path>) -> Result<()> {
        std::fs::write(path, self.to_json()).map_err(Error::Io)
    }

    /// The charge distribution collected by [`SwathQC::process_spectrum`].
    ///
    /// Source member `getChargeDistribution()`.
    pub fn charge_distribution(&self) -> &ChargeDistribution {
        &self.cd
    }

    /// Set the number of expected MS1 spectra, which
    /// [`SwathQC::process_spectrum`] needs to sample uniformly; 0 inspects
    /// every spectrum passed in. Source `setNrMS1Spectra`.
    pub fn set_nr_ms1_spectra(&mut self, nr: usize) {
        self.nr_ms1_spectra = nr;
    }

    /// Whether spectrum `idx` is one of `subsample_count` spectra sampled
    /// uniformly from `total_spec_count`.
    ///
    /// Source protected static `isSubsampledSpectrum_`, public here because
    /// it is a pure function the class test calls. With `total_spec_count`
    /// 0 (unknown) every index is sampled; an index at or beyond the total, or
    /// a `subsample_count` of 0, is never sampled. Otherwise, with step
    /// `total / min(total, subsample)` (at least 1), `idx` is sampled when it
    /// equals the rounded (half away from zero) `floor(idx / step) * step` or
    /// `ceil(idx / step) * step`. For 4 of 10 spectra this samples indices 0,
    /// 3, 5 and 8.
    pub fn is_subsampled_spectrum(
        total_spec_count: usize,
        subsample_count: usize,
        idx: usize,
    ) -> bool {
        if total_spec_count == 0 {
            return true;
        }
        if idx >= total_spec_count || subsample_count == 0 {
            return false;
        }
        let spec_count = total_spec_count as f64;
        let step_size = spec_count / spec_count.min(subsample_count as f64);
        let steps = idx as f64 / step_size;
        let steps_low = steps.floor() * step_size;
        let steps_high = steps.ceil() * step_size;
        let target = idx as f64;
        steps_low.round() == target || steps_high.round() == target
    }
}

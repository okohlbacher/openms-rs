// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Predicates for range operations on spectra and peaks.
//!
//! Source: `OpenMS/KERNEL/RangeUtils.h`, a header-only group of fifteen
//! functor templates. Support document: `docs/RANGE_UTILS_SUPPORT.md`.
//!
//! Each predicate answers the source `operator()` through
//! [`SpectrumPredicate::keep`](crate::kernel::range_utils::SpectrumPredicate::keep)
//! or [`PeakPredicate::keep`](crate::kernel::range_utils::PeakPredicate::keep),
//! and carries the source `reverse` flag that inverts the answer. The source
//! pairs the functors with `std::remove_if`/`std::erase_if`, which *removes*
//! matches;
//! [`MSExperiment::retain_spectra`](crate::kernel::MSExperiment::retain_spectra)
//! and
//! [`MSSpectrum::retain_peaks_where`](crate::kernel::MSSpectrum::retain_peaks_where)
//! *keep* matches, like `Vec::retain`. A source `erase_if(spectra, P(args))`
//! therefore becomes `retain_spectra(&P::new(args, /* reverse */ true))`, or
//! equivalently the negated closure.
//!
//! Two source predicates have no struct here because the crate already covers
//! them: `InRTRange` is
//! [`MSExperiment::spectra_in_rt_range`](crate::kernel::MSExperiment::spectra_in_rt_range)
//! for a sorted borrowing query and a closure on `rt` for removal;
//! `InMSLevelRange` is
//! [`MSExperiment::ms_levels`](crate::kernel::MSExperiment::ms_levels) /
//! [`MSExperiment::contains_scan_of_level`](crate::kernel::MSExperiment::contains_scan_of_level)
//! plus a closure on `ms_level`. Closures implement
//! [`SpectrumPredicate`](crate::kernel::range_utils::SpectrumPredicate) directly.
//!
//! The source example removing spectra in a retention-time range, rewritten:
//!
//! ```rust
//! use openms::kernel::{MSExperiment, MSSpectrum};
//!
//! let mut experiment = MSExperiment::new();
//! for rt in [10.0, 30.0, 50.0] {
//!     let mut spectrum = MSSpectrum::new();
//!     spectrum.rt = rt;
//!     experiment.spectra.push(spectrum);
//! }
//! // Source: `std::erase_if(spectra, InRTRange<MSSpectrum>(0.0, 36.0))`.
//! // Removal keeps what lies *outside* the closed range 0.0..=36.0 s.
//! let removed = experiment.retain_spectra(&|s: &MSSpectrum| !(0.0 <= s.rt && s.rt <= 36.0))?;
//! assert_eq!(removed, 2);
//! assert_eq!(experiment.spectra[0].rt, 50.0);
//! # Ok::<(), openms::Error>(())
//! ```
//!
//! The source example removing peaks in an intensity range, rewritten:
//!
//! ```rust
//! use openms::kernel::range_utils::InIntensityRange;
//! use openms::kernel::{MSSpectrum, Peak1D};
//!
//! let mut spectrum = MSSpectrum::from_peaks(vec![
//!     Peak1D::new(100.0, 10.0),
//!     Peak1D::new(200.0, 6000.0),
//! ]);
//! // Source: `std::erase_if(spectrum, InIntensityRange<Peak1D>(0.0, 5000.0))`.
//! let removed = spectrum.retain_peaks_where(&InIntensityRange::new(0.0, 5000.0, true)?)?;
//! assert_eq!(removed, 1);
//! assert_eq!(spectrum.peaks[0].intensity, 6000.0);
//! # Ok::<(), openms::Error>(())
//! ```

use super::{MSExperiment, MSSpectrum, Peak1D, Precursor};
use crate::error::{Error, Result};
use crate::metadata::{ActivationMethod, MetaInfo, Polarity, ScanMode};
use std::collections::BTreeSet;
use std::str::FromStr;

/// Metadata key under which the source mzML handler stores a precursor's
/// collision energy (`MS:1000045`), read by [`IsInCollisionEnergyRange`].
pub const COLLISION_ENERGY_KEY: &str = "collision energy";

/// PSI-MS accession of the collision-energy term, consulted when the source
/// metadata key is absent but the term was retained on the precursor.
pub const COLLISION_ENERGY_ACCESSION: &str = "MS:1000045";

/// A yes/no decision about one spectrum; the source functor `operator()`.
///
/// `true` means the spectrum *matches* the predicate. Removal in the source
/// (`remove_if`) drops matches; [`MSExperiment::retain_spectra`] keeps them.
/// Any `Fn(&MSSpectrum) -> bool` closure is a predicate.
pub trait SpectrumPredicate {
    /// Whether `spectrum` satisfies the predicate, after the `reverse` flag.
    fn keep(&self, spectrum: &MSSpectrum) -> bool;
}

impl<F: Fn(&MSSpectrum) -> bool> SpectrumPredicate for F {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self(spectrum)
    }
}

/// A yes/no decision about one peak; the source functor `operator()` for the
/// `PeakType` templates. Any `Fn(&Peak1D) -> bool` closure is a predicate.
pub trait PeakPredicate {
    /// Whether `peak` satisfies the predicate, after the `reverse` flag.
    fn keep(&self, peak: &Peak1D) -> bool;
}

impl<F: Fn(&Peak1D) -> bool> PeakPredicate for F {
    fn keep(&self, peak: &Peak1D) -> bool {
        self(peak)
    }
}

/// Ceiling on the number of spectra or peaks one filtering call may visit.
///
/// Native only: the source functors carry no bound. The count is checked
/// before any mutation, so a rejected call leaves its container unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeFilterLimits {
    /// Maximum spectra (for [`MSExperiment::retain_spectra_with_limits`]) or
    /// peaks (for [`MSSpectrum::retain_peaks_where_with_limits`]) per call.
    pub max_items: usize,
}

impl Default for RangeFilterLimits {
    fn default() -> Self {
        Self {
            max_items: 50_000_000,
        }
    }
}

fn preflight(count: usize, what: &str, limits: RangeFilterLimits) -> Result<()> {
    if count > limits.max_items {
        return Err(Error::InvalidValue(format!(
            "{what} count {count} exceeds the filter limit {}",
            limits.max_items
        )));
    }
    Ok(())
}

/// Reject nonfinite or inverted closed-range bounds before a predicate exists.
/// The source stores any pair silently; an inverted pair then never matches.
fn closed_range(min: f64, max: f64, what: &str) -> Result<()> {
    if !min.is_finite() || !max.is_finite() {
        return Err(Error::InvalidValue(format!("{what} bounds must be finite")));
    }
    if min > max {
        return Err(Error::InvalidRange(format!(
            "{what} minimum {min} exceeds maximum {max}"
        )));
    }
    Ok(())
}

impl MSExperiment {
    /// Keep the spectra for which `predicate` is `true`, returning how many
    /// were removed. Chromatograms and settings are untouched.
    ///
    /// The source applies its functors through `std::erase_if`, which removes
    /// matches; pass the predicate with `reverse = true` (or negate a closure)
    /// to reproduce that. Predicates are infallible, so removal cannot stop
    /// halfway.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the spectrum count exceeds
    /// [`RangeFilterLimits::default`], before anything is removed.
    pub fn retain_spectra<P: SpectrumPredicate + ?Sized>(
        &mut self,
        predicate: &P,
    ) -> Result<usize> {
        self.retain_spectra_with_limits(predicate, RangeFilterLimits::default())
    }

    /// As [`Self::retain_spectra`], with an explicit spectrum-count ceiling.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `spectra.len()` exceeds
    /// `limits.max_items`; the experiment is then unchanged.
    pub fn retain_spectra_with_limits<P: SpectrumPredicate + ?Sized>(
        &mut self,
        predicate: &P,
        limits: RangeFilterLimits,
    ) -> Result<usize> {
        preflight(self.spectra.len(), "spectrum", limits)?;
        let before = self.spectra.len();
        self.spectra.retain(|spectrum| predicate.keep(spectrum));
        Ok(before - self.spectra.len())
    }
}

impl MSSpectrum {
    /// Keep the peaks for which `predicate` is `true`, moving every aligned
    /// annotation array with them, and return how many were removed.
    ///
    /// This is [`MSSpectrum::retain_peaks`] driven by a [`PeakPredicate`];
    /// the source pairs `InMzRange`/`InIntensityRange` with `std::erase_if`,
    /// which removes matches, so pass `reverse = true` to reproduce that.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the peak count exceeds
    /// [`RangeFilterLimits::default`] or when a nonempty data array does not
    /// have one entry per peak. Either error leaves the spectrum unchanged.
    pub fn retain_peaks_where<P: PeakPredicate + ?Sized>(
        &mut self,
        predicate: &P,
    ) -> Result<usize> {
        self.retain_peaks_where_with_limits(predicate, RangeFilterLimits::default())
    }

    /// As [`Self::retain_peaks_where`], with an explicit peak-count ceiling.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `peaks.len()` exceeds
    /// `limits.max_items` or an aligned array is inconsistent; the spectrum is
    /// then unchanged.
    pub fn retain_peaks_where_with_limits<P: PeakPredicate + ?Sized>(
        &mut self,
        predicate: &P,
        limits: RangeFilterLimits,
    ) -> Result<usize> {
        preflight(self.peaks.len(), "peak", limits)?;
        let before = self.peaks.len();
        self.retain_peaks(|peak| predicate.keep(peak))?;
        Ok(before - self.peaks.len())
    }
}

/// Predicate that determines if a record has a certain metavalue.
///
/// The source template accepts any `MetaInfoInterface`; here the spectrum
/// form is the [`SpectrumPredicate`] and [`HasMetaValue::evaluate`] serves
/// every other owner of a [`MetaInfo`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasMetaValue {
    key: String,
    reverse: bool,
}

impl HasMetaValue {
    /// `key` is the metavalue that needs to be present. If `reverse` is true,
    /// the predicate is true when the metavalue does not exist.
    pub fn new(key: impl Into<String>, reverse: bool) -> Self {
        Self {
            key: key.into(),
            reverse,
        }
    }

    /// The decision for any metadata map, not only a spectrum's.
    pub fn evaluate(&self, metadata: &MetaInfo) -> bool {
        self.reverse ^ metadata.contains_key(&self.key)
    }
}

impl SpectrumPredicate for HasMetaValue {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self.evaluate(&spectrum.metadata)
    }
}

/// Predicate that determines if a spectrum has a certain scan mode.
///
/// The source compares an integer against the cast enumerator; the port
/// compares [`ScanMode`] values directly, so no out-of-range integer exists.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasScanMode {
    mode: ScanMode,
    reverse: bool,
}

impl HasScanMode {
    /// `mode` is the scan mode to look for. If `reverse` is true, the
    /// predicate is true when the spectrum has a different scan mode.
    pub fn new(mode: ScanMode, reverse: bool) -> Self {
        Self { mode, reverse }
    }
}

impl SpectrumPredicate for HasScanMode {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self.reverse ^ (spectrum.instrument_settings.scan_mode == self.mode)
    }
}

/// Predicate that determines if a spectrum has a certain scan polarity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasScanPolarity {
    polarity: Polarity,
    reverse: bool,
}

impl HasScanPolarity {
    /// `polarity` is the scan polarity to look for. If `reverse` is true, the
    /// predicate is true when the spectrum has a different scan polarity.
    pub fn new(polarity: Polarity, reverse: bool) -> Self {
        Self { polarity, reverse }
    }
}

impl SpectrumPredicate for HasScanPolarity {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self.reverse ^ (spectrum.instrument_settings.polarity == self.polarity)
    }
}

/// Predicate that determines if a spectrum is empty (has no peaks).
///
/// Only the peak list counts, as the source `empty()`; metadata, precursors
/// and data arrays do not make a spectrum nonempty.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IsEmptySpectrum {
    reverse: bool,
}

impl IsEmptySpectrum {
    /// If `reverse` is true, the predicate is true when the spectrum is not
    /// empty.
    pub fn new(reverse: bool) -> Self {
        Self { reverse }
    }
}

impl SpectrumPredicate for IsEmptySpectrum {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self.reverse ^ spectrum.peaks.is_empty()
    }
}

/// Predicate that determines if a spectrum is a zoom (enhanced resolution)
/// spectrum, read from its instrument settings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IsZoomSpectrum {
    reverse: bool,
}

impl IsZoomSpectrum {
    /// If `reverse` is true, the predicate is true when the spectrum is not a
    /// zoom spectrum.
    pub fn new(reverse: bool) -> Self {
        Self { reverse }
    }
}

impl SpectrumPredicate for IsZoomSpectrum {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        self.reverse ^ spectrum.instrument_settings.zoom_scan
    }
}

/// Predicate that determines if a spectrum was generated using any activation
/// method given in the constructor list.
///
/// Any precursor carrying any listed method matches. The source compares the
/// long names in `Precursor::NamesOfActivationMethod` as strings; the port
/// holds typed [`ActivationMethod`] values, and [`HasActivationMethod::from_names`]
/// resolves names the way the source list did.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasActivationMethod {
    methods: BTreeSet<ActivationMethod>,
    reverse: bool,
}

impl HasActivationMethod {
    /// `methods` are compared against the precursor activation methods. If
    /// `reverse` is true, the predicate is true when the spectrum is not using
    /// one of the specified activation methods.
    pub fn new(methods: impl IntoIterator<Item = ActivationMethod>, reverse: bool) -> Self {
        Self {
            methods: methods.into_iter().collect(),
            reverse,
        }
    }

    /// Build from method names: the source long names
    /// (`"Collision-induced dissociation"`) or the short names (`"CID"`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a name that is neither. The source
    /// keeps unknown strings and they simply never match; here a misspelt
    /// method is reported instead of silently filtering nothing.
    pub fn from_names(
        names: impl IntoIterator<Item = impl AsRef<str>>,
        reverse: bool,
    ) -> Result<Self> {
        let methods = names
            .into_iter()
            .map(|name| ActivationMethod::from_str(name.as_ref()))
            .collect::<Result<BTreeSet<_>>>()?;
        Ok(Self { methods, reverse })
    }
}

impl SpectrumPredicate for HasActivationMethod {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        let found = spectrum.precursors.iter().any(|precursor| {
            precursor
                .activation_methods
                .iter()
                .any(|method| self.methods.contains(method))
        });
        self.reverse ^ found
    }
}

/// Predicate that determines if a spectrum's precursor is within a certain
/// m/z range.
///
/// If multiple precursors are present, all must fulfill the closed-interval
/// criterion; a spectrum without precursors fulfills it vacuously.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InPrecursorMZRange {
    mz_left: f64,
    mz_right: f64,
    reverse: bool,
}

impl InPrecursorMZRange {
    /// `mz_left` and `mz_right` are the closed-interval boundaries. If
    /// `reverse` is true, the predicate is true when a precursor's m/z is
    /// outside the interval.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite bound and
    /// [`Error::InvalidRange`] when `mz_left > mz_right`; the source accepts
    /// both silently.
    pub fn new(mz_left: f64, mz_right: f64, reverse: bool) -> Result<Self> {
        closed_range(mz_left, mz_right, "precursor m/z")?;
        Ok(Self {
            mz_left,
            mz_right,
            reverse,
        })
    }
}

impl SpectrumPredicate for InPrecursorMZRange {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        let all_inside = spectrum
            .precursors
            .iter()
            .all(|precursor| self.mz_left <= precursor.mz && precursor.mz <= self.mz_right);
        self.reverse ^ all_inside
    }
}

/// Predicate that determines if a spectrum has a certain precursor charge as
/// given in the constructor list. Any precursor with a listed charge matches.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasPrecursorCharge {
    charges: BTreeSet<i32>,
    reverse: bool,
}

impl HasPrecursorCharge {
    /// `charges` are compared against each precursor charge. If `reverse` is
    /// true, the predicate is true when the spectrum has none of the specified
    /// precursor charges.
    pub fn new(charges: impl IntoIterator<Item = i32>, reverse: bool) -> Self {
        Self {
            charges: charges.into_iter().collect(),
            reverse,
        }
    }
}

impl SpectrumPredicate for HasPrecursorCharge {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        let found = spectrum
            .precursors
            .iter()
            .any(|precursor| self.charges.contains(&precursor.charge));
        self.reverse ^ found
    }
}

/// Predicate that determines if a peak lies inside/outside a specific closed
/// m/z range.
///
/// The source note assumes the m/z dimension is dimension 0 of the position;
/// [`Peak1D::mz`] is that dimension, so nothing is assumed here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InMzRange {
    min: f64,
    max: f64,
    reverse: bool,
}

impl InMzRange {
    /// `min` and `max` are the closed boundaries. If `reverse` is true, the
    /// predicate is true when the peak lies outside the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite bound and
    /// [`Error::InvalidRange`] when `min > max`; the source accepts both.
    pub fn new(min: f64, max: f64, reverse: bool) -> Result<Self> {
        closed_range(min, max, "m/z")?;
        Ok(Self { min, max, reverse })
    }

    /// The decision for a bare m/z, for peak types other than [`Peak1D`].
    pub fn contains(&self, mz: f64) -> bool {
        self.reverse ^ (self.min <= mz && mz <= self.max)
    }
}

impl PeakPredicate for InMzRange {
    fn keep(&self, peak: &Peak1D) -> bool {
        self.contains(peak.mz)
    }
}

/// Predicate that determines if a peak lies inside/outside a specific closed
/// intensity range.
///
/// The `f32` intensity is widened to `f64` before comparison, as the source
/// `double tmp = p.getIntensity()`; `4.9f32` is therefore below `5.0` and
/// `10.1f32` above `10.0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InIntensityRange {
    min: f64,
    max: f64,
    reverse: bool,
}

impl InIntensityRange {
    /// `min` and `max` are the closed boundaries. If `reverse` is true, the
    /// predicate is true when the peak lies outside the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite bound and
    /// [`Error::InvalidRange`] when `min > max`; the source accepts both.
    pub fn new(min: f64, max: f64, reverse: bool) -> Result<Self> {
        closed_range(min, max, "intensity")?;
        Ok(Self { min, max, reverse })
    }

    /// The decision for a bare intensity, for peak types other than [`Peak1D`].
    pub fn contains(&self, intensity: f32) -> bool {
        let value = f64::from(intensity);
        self.reverse ^ (self.min <= value && value <= self.max)
    }
}

impl PeakPredicate for InIntensityRange {
    fn keep(&self, peak: &Peak1D) -> bool {
        self.contains(peak.intensity)
    }
}

/// A precursor's collision energy as the source predicate reads it.
///
/// The source reads the metavalue `"collision energy"`, which its mzML handler
/// fills from `MS:1000045`. The port looks up [`COLLISION_ENERGY_KEY`] in
/// [`Precursor::cv_terms`]' metadata first and then the retained
/// [`COLLISION_ENERGY_ACCESSION`] term's value. Non-numeric values count as
/// absent, where the source `DataValue` conversion would misbehave.
pub fn collision_energy(precursor: &Precursor) -> Option<f64> {
    if let Some(value) = precursor.cv_terms.metadata.get(COLLISION_ENERGY_KEY) {
        return value.as_f64().ok();
    }
    precursor
        .cv_terms
        .get(COLLISION_ENERGY_ACCESSION)
        .and_then(|terms| terms.iter().find_map(|term| term.value.as_f64().ok()))
}

/// Predicate that determines if an MSn spectrum was generated with a collision
/// energy in the given closed range.
///
/// This applies only to spectra whose precursors record a collision energy
/// (CID and HCD in practice). The source notes say the predicate "returns
/// true" for MS1 spectra and for spectra without a collision energy; its code
/// returns **false** for both, whatever `reverse` is, and the port preserves
/// the code. Because the source drives this predicate with `remove_if`, that
/// `false` is what keeps those spectra, which is what the notes mean.
///
/// Any precursor in range makes the spectrum match.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsInCollisionEnergyRange {
    min_energy: f64,
    max_energy: f64,
    reverse: bool,
}

impl IsInCollisionEnergyRange {
    /// `min` and `max` are the minimum and maximum collision energy included
    /// in the range. If `reverse` is true, the predicate is true when the
    /// collision energy lies outside the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite bound and
    /// [`Error::InvalidRange`] when `min > max`; the source accepts both.
    pub fn new(min: f64, max: f64, reverse: bool) -> Result<Self> {
        closed_range(min, max, "collision energy")?;
        Ok(Self {
            min_energy: min,
            max_energy: max,
            reverse,
        })
    }
}

impl SpectrumPredicate for IsInCollisionEnergyRange {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        // Leave non-fragmentation spectra untouched; the source tests level 1 only.
        if spectrum.ms_level == 1 {
            return false;
        }
        let mut is_in = false;
        let mut has_collision_energy = false;
        for precursor in &spectrum.precursors {
            if let Some(energy) = collision_energy(precursor) {
                has_collision_energy = true;
                is_in |= !(energy > self.max_energy || energy < self.min_energy);
            }
        }
        // The source accepts every spectrum without a collision energy value.
        if !has_collision_energy {
            return false;
        }
        self.reverse ^ is_in
    }
}

/// Predicate that determines if the width of the isolation window of an MSn
/// spectrum is in the given closed range.
///
/// The width is the sum of the lower and upper isolation offsets. The source
/// note says the predicate "returns true" for MS1 spectra; its code returns
/// **false** regardless of `reverse`, which is what keeps them under
/// `remove_if`, and the port preserves the code. Unlike
/// [`IsInCollisionEnergyRange`], a spectrum without precursors is not special:
/// nothing is in range, so the result is `reverse`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsInIsolationWindowSizeRange {
    min_size: f64,
    max_size: f64,
    reverse: bool,
}

impl IsInIsolationWindowSizeRange {
    /// `min_size` and `max_size` are the minimum and maximum isolation window
    /// widths. If `reverse` is true, the predicate is true when the width lies
    /// outside the range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite bound and
    /// [`Error::InvalidRange`] when `min_size > max_size`; the source accepts
    /// both.
    pub fn new(min_size: f64, max_size: f64, reverse: bool) -> Result<Self> {
        closed_range(min_size, max_size, "isolation window width")?;
        Ok(Self {
            min_size,
            max_size,
            reverse,
        })
    }
}

impl SpectrumPredicate for IsInIsolationWindowSizeRange {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        if spectrum.ms_level == 1 {
            return false;
        }
        let mut is_in = false;
        for precursor in &spectrum.precursors {
            let width =
                precursor.isolation_window_upper_offset + precursor.isolation_window_lower_offset;
            is_in |= !(width > self.max_size || width < self.min_size);
        }
        self.reverse ^ is_in
    }
}

/// Predicate that determines if the isolation window covers ANY of the given
/// m/z values.
///
/// For each precursor the window is `[mz - lower_offset, mz + upper_offset]`;
/// the smallest listed value at or above the lower edge must lie at or below
/// the upper edge. The source note says the predicate "returns true" for MS1
/// spectra; its code returns **false** regardless of `reverse`, which keeps
/// them under `remove_if`, and the port preserves the code.
///
/// The source logs a warning for every precursor whose lower or upper offset
/// is zero ("Filtering will probably be too strict (unless you hit the exact
/// precursor m/z)!") and then evaluates normally. The port evaluates the same
/// way but emits no warning; the module has no log sink.
#[derive(Clone, Debug, PartialEq)]
pub struct IsInIsolationWindow {
    mz_values: Vec<f64>,
    reverse: bool,
}

impl IsInIsolationWindow {
    /// `mz_values` are the m/z values of which at least one needs to be
    /// covered; they are sorted on construction, as in the source. If
    /// `reverse` is true, the predicate is true when the isolation window
    /// misses ALL m/z values.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite value; the source would
    /// sort NaN into an unspecified position and search past it.
    pub fn new(mz_values: impl IntoIterator<Item = f64>, reverse: bool) -> Result<Self> {
        let mut mz_values: Vec<f64> = mz_values.into_iter().collect();
        if mz_values.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidValue(
                "isolation window m/z values must be finite".into(),
            ));
        }
        mz_values.sort_by(f64::total_cmp);
        Ok(Self { mz_values, reverse })
    }

    /// The m/z values in the sorted order the search uses.
    pub fn mz_values(&self) -> &[f64] {
        &self.mz_values
    }
}

impl SpectrumPredicate for IsInIsolationWindow {
    fn keep(&self, spectrum: &MSSpectrum) -> bool {
        if spectrum.ms_level == 1 {
            return false;
        }
        let mut is_in = false;
        for precursor in &spectrum.precursors {
            let lower_mz = precursor.mz - precursor.isolation_window_lower_offset;
            // Source: std::lower_bound, the first value not below the lower edge.
            let index = self.mz_values.partition_point(|&value| value < lower_mz);
            if let Some(&candidate) = self.mz_values.get(index) {
                let upper_mz = precursor.mz + precursor.isolation_window_upper_offset;
                is_in |= candidate <= upper_mz;
            }
        }
        self.reverse ^ is_in
    }
}

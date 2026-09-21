// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{CVTermList, MetaInfo, MetaMergePolicy, merge_meta, validate_meta};
use crate::kernel::{NumericRange, Precursor, SpectrumType};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::{fmt, str::FromStr};

macro_rules! named_enum {
    ($(#[$enum_attr:meta])* $name:ident { $($(#[$variant_attr:meta])* $variant:ident => $label:literal),+ $(,)? }) => {
        $(#[$enum_attr])*
        /// Enumeration with one OpenMS string label per variant, in source order.
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
        pub enum $name { $($(#[$variant_attr])* $variant),+ }
        impl $name {
            /// Every variant, in the source's declaration order.
            pub const ALL: &'static [Self] = &[$(Self::$variant),+];
            /// The variant's OpenMS label, which its `Display` and `FromStr` also use.
            pub const fn name(self) -> &'static str { match self { $(Self::$variant => $label),+ } }
        }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str(self.name()) }
        }
        impl FromStr for $name {
            type Err = Error;
            fn from_str(value: &str) -> Result<Self> {
                Self::ALL.iter().copied().find(|item| item.name() == value)
                    .ok_or_else(|| invalid(concat!("unknown ",stringify!($name))))
            }
        }
    };
}

named_enum!(#[derive(Default)] ScanMode {
    #[default] Unknown => "Unknown", MassSpectrum => "MassSpectrum", Ms1Spectrum => "MS1Spectrum", MsnSpectrum => "MSnSpectrum",
    SelectedIonMonitoring => "SelectedIonMonitoring", SelectedReactionMonitoring => "SelectedReactionMonitoring",
    ConsecutiveReactionMonitoring => "ConsecutiveReactionMonitoring", ConstantNeutralGain => "ConstantNeutralGain",
    ConstantNeutralLoss => "ConstantNeutralLoss", Precursor => "Precursor", EnhancedMultiplyCharged => "EnhancedMultiplyCharged",
    TimeDelayedFragmentation => "TimeDelayedFragmentation", ElectromagneticRadiation => "ElectromagneticRadiation",
    Emission => "Emission", Absorption => "Absorption"
});
named_enum!(#[derive(Default)] Polarity { #[default] Unknown => "unknown", Positive => "positive", Negative => "negative" });
named_enum!(#[derive(Default)] DriftTimeUnit { #[default] None => "<NONE>", Millisecond => "ms", InverseReducedMobility => "1/K0", FaimsCompensationVoltage => "FAIMS_CV", CollisionCrossSection => "CCS" });
named_enum!(
    #[doc = "How a spectrum represents ion mobility, the source `IMFormat` (`IONMOBILITY/IMTypes.h:43-53`). `Unknown`, the default, means not yet determined: the file handler or the ion-mobility peak picker is expected to set a known value, and `ImTypes::determine_im_format` derives one from the data. See `docs/IM_TYPES_SUPPORT.md`."]
    #[derive(Default)] IonMobilityFormat { None => "none", PerPeak => "im_peak", PerSpectrum => "im_spectrum", #[default] Unknown => "unknown" });
named_enum!(
    #[doc = "Processing state of a spectrum's ion-mobility dimension, the source `IMPeakType` (`IONMOBILITY/IMTypes.h:64-72`); the analogue of `SpectrumType` for the m/z dimension. `Profile` is raw data such as a full TIMS frame before ion-mobility centroiding, `Centroid` is data centroided in the ion-mobility dimension, and `Unknown`, the default, means not yet determined."]
    #[derive(Default)] IonMobilityPeakType { Profile => "im_profile", Centroid => "im_centroided", #[default] Unknown => "unknown" });
named_enum!(#[derive(Default)] ChecksumType { #[default] Unknown => "Unknown", Sha1 => "SHA-1", Md5 => "MD5" });
named_enum!(#[derive(Default)] ChromatogramType {
    #[default] Mass => "mass chromatogram", TotalIonCurrent => "total ion current chromatogram", SelectedIonCurrent => "selected ion current chromatogram",
    BasePeak => "base peak chromatogram", SelectedIonMonitoring => "selected ion monitoring chromatogram", SelectedReactionMonitoring => "selected reaction monitoring chromatogram",
    ElectromagneticRadiation => "electromagnetic radiation chromatogram", Absorption => "absorption chromatogram", Emission => "emission chromatogram", Unknown => "unknown chromatogram"
});
named_enum!(ProcessingAction {
    DataProcessing => "Data processing action", ChargeDeconvolution => "Charge deconvolution", Deisotoping => "Deisotoping", Smoothing => "Smoothing",
    ChargeCalculation => "Charge calculation", PrecursorRecalculation => "Precursor recalculation", BaselineReduction => "Baseline reduction",
    PeakPicking => "Peak picking", RetentionTimeAlignment => "Retention time alignment", MzCalibration => "Calibration of m/z positions",
    IntensityNormalization => "Intensity normalization", DataFiltering => "Data filtering", Quantitation => "Quantitation", FeatureGrouping => "Feature grouping",
    IdentificationMapping => "Identification mapping", FormatConversion => "File format conversion", ConversionMzData => "Conversion to mzData format",
    ConversionMzML => "Conversion to mzML format", ConversionMzXML => "Conversion to mzXML format", ConversionDta => "Conversion to DTA format",
    Identification => "Identification", IonMobilityBinning => "Ion mobility binning"
});

/// The 19 activation methods in pinned OpenMS order. No invalid size sentinel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ActivationMethod {
    Cid,
    Psd,
    Pd,
    Sid,
    Bird,
    Ecd,
    Imd,
    Sori,
    Hcid,
    Lcid,
    Phd,
    Etd,
    Etcid,
    Ethcd,
    Pqd,
    Trap,
    Hcd,
    InSource,
    Lift,
}
impl ActivationMethod {
    /// Every activation method, in the pinned OpenMS order.
    pub const ALL: &'static [Self] = &[
        Self::Cid,
        Self::Psd,
        Self::Pd,
        Self::Sid,
        Self::Bird,
        Self::Ecd,
        Self::Imd,
        Self::Sori,
        Self::Hcid,
        Self::Lcid,
        Self::Phd,
        Self::Etd,
        Self::Etcid,
        Self::Ethcd,
        Self::Pqd,
        Self::Trap,
        Self::Hcd,
        Self::InSource,
        Self::Lift,
    ];
    /// The method's abbreviation, such as `CID`; `FromStr` accepts it as well
    /// as the full [`ActivationMethod::name`].
    pub const fn short_name(self) -> &'static str {
        match self {
            Self::Cid => "CID",
            Self::Psd => "PSD",
            Self::Pd => "PD",
            Self::Sid => "SID",
            Self::Bird => "BIRD",
            Self::Ecd => "ECD",
            Self::Imd => "IMD",
            Self::Sori => "SORI",
            Self::Hcid => "HCID",
            Self::Lcid => "LCID",
            Self::Phd => "PHD",
            Self::Etd => "ETD",
            Self::Etcid => "ETciD",
            Self::Ethcd => "EThcD",
            Self::Pqd => "PQD",
            Self::Trap => "TRAP",
            Self::Hcd => "HCD",
            Self::InSource => "INSOURCE",
            Self::Lift => "LIFT",
        }
    }
    /// The method's full name, as written by the mzML transport.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Cid => "Collision-induced dissociation",
            Self::Psd => "Post-source decay",
            Self::Pd => "Plasma desorption",
            Self::Sid => "Surface-induced dissociation",
            Self::Bird => "Blackbody infrared radiative dissociation",
            Self::Ecd => "Electron capture dissociation",
            Self::Imd => "Infrared multiphoton dissociation",
            Self::Sori => "Sustained off-resonance irradiation",
            Self::Hcid => "High-energy collision-induced dissociation",
            Self::Lcid => "Low-energy collision-induced dissociation",
            Self::Phd => "Photodissociation",
            Self::Etd => "Electron transfer dissociation",
            Self::Etcid => "Electron transfer and collision-induced dissociation",
            Self::Ethcd => "Electron transfer and higher-energy collision dissociation",
            Self::Pqd => "Pulsed q dissociation",
            Self::Trap => "trap-type collision-induced dissociation",
            Self::Hcd => "beam-type collision-induced dissociation",
            Self::InSource => "in-source collision-induced dissociation",
            Self::Lift => "Bruker proprietary method",
        }
    }
}
impl fmt::Display for ActivationMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}
impl FromStr for ActivationMethod {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        Self::ALL
            .iter()
            .copied()
            .find(|method| method.name() == value || method.short_name() == value)
            .ok_or_else(|| invalid("unknown activation method"))
    }
}

/// Compatibility wrapper for the precursor record now stored directly in spectra.
/// Dereferencing accesses the same acquisition fields; there is no duplicate state.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PrecursorInfo {
    pub peak: Precursor,
}
impl From<Precursor> for PrecursorInfo {
    fn from(peak: Precursor) -> Self {
        Self { peak }
    }
}
impl From<PrecursorInfo> for Precursor {
    fn from(info: PrecursorInfo) -> Self {
        info.peak
    }
}
impl std::ops::Deref for PrecursorInfo {
    type Target = Precursor;
    fn deref(&self) -> &Precursor {
        &self.peak
    }
}
impl std::ops::DerefMut for PrecursorInfo {
    fn deref_mut(&mut self) -> &mut Precursor {
        &mut self.peak
    }
}
impl Precursor {
    /// Check the acquisition fields and the attached CV terms.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite m/z, intensity or drift
    /// value, a negative isolation target, activation energy or window offset, a
    /// parent spectrum reference that is empty or holds a control character, or
    /// an invalid CV term.
    pub fn validate(&self) -> Result<()> {
        finite(self.mz, "precursor m/z")?;
        finite(f64::from(self.intensity), "precursor intensity")?;
        if let Some(target) = self.isolation_target_mz {
            nonnegative(target, "isolation target m/z")?;
        }
        nonnegative(self.activation_energy, "activation energy")?;
        for value in [
            self.isolation_window_lower_offset,
            self.isolation_window_upper_offset,
            self.drift_window_lower_offset,
            self.drift_window_upper_offset,
        ] {
            nonnegative(value, "isolation offset")?;
        }
        if let Some(drift) = self.drift_time {
            finite(drift, "drift time or voltage")?;
        }
        if self
            .spectrum_reference
            .as_ref()
            .is_some_and(|s| s.is_empty() || s.chars().any(char::is_control))
        {
            return Err(invalid(
                "parent spectrum reference must be nonempty and control-free",
            ));
        }
        self.cv_terms.validate()
    }
    /// The isolation window as an absolute m/z range, the target m/z widened by
    /// the two stored offsets.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when [`Precursor::validate`] fails or
    /// either bound is not finite.
    pub fn isolation_window(&self) -> Result<NumericRange> {
        self.validate()?;
        checked_window(
            self.mz,
            self.isolation_window_lower_offset,
            self.isolation_window_upper_offset,
        )
    }
    /// Preserves getUnchargedMass: unknown charge assumes 2, and negative charge
    /// is used with its sign. Use feature decharging for an absolute-charge mass.
    pub fn uncharged_mass(&self) -> Result<f64> {
        self.validate()?;
        let charge = f64::from(if self.charge == 0 { 2 } else { self.charge });
        let mass = self.mz * charge - charge * crate::constants::PROTON_MASS_U;
        finite(mass, "uncharged mass")?;
        Ok(mass)
    }
}

/// Isolation description of a product ion: a target m/z with its window.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Product {
    pub mz: f64,
    pub isolation_window_lower_offset: f64,
    pub isolation_window_upper_offset: f64,
    pub cv_terms: CVTermList,
}
impl Hash for Product {
    fn hash<H: Hasher>(&self, state: &mut H) {
        super::value::hash_float(self.mz, state);
        super::value::hash_float(self.isolation_window_lower_offset, state);
        super::value::hash_float(self.isolation_window_upper_offset, state);
        self.cv_terms.hash(state);
    }
}
impl Product {
    /// Check the m/z, both window offsets and the attached CV terms.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite m/z, a negative or
    /// nonfinite window offset, or an invalid CV term.
    pub fn validate(&self) -> Result<()> {
        finite(self.mz, "product m/z")?;
        nonnegative(self.isolation_window_lower_offset, "lower isolation offset")?;
        nonnegative(self.isolation_window_upper_offset, "upper isolation offset")?;
        self.cv_terms.validate()
    }
    /// The isolation window as an absolute m/z range, the target m/z widened by
    /// the two stored offsets.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when [`Product::validate`] fails or
    /// either bound is not finite.
    pub fn isolation_window(&self) -> Result<NumericRange> {
        self.validate()?;
        checked_window(
            self.mz,
            self.isolation_window_lower_offset,
            self.isolation_window_upper_offset,
        )
    }
}

/// One inclusive m/z acquisition window of a scan.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScanWindow {
    pub begin: f64,
    pub end: f64,
    pub metadata: MetaInfo,
}
impl ScanWindow {
    /// A window from `begin` to `end`, with no metadata.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either bound is not finite or
    /// `begin` exceeds `end`.
    pub fn new(begin: f64, end: f64) -> Result<Self> {
        let result = Self {
            begin,
            end,
            ..Self::default()
        };
        result.validate()?;
        Ok(result)
    }
    /// Check the bounds and the metadata.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either bound is not finite, `begin`
    /// exceeds `end`, or the metadata is invalid.
    pub fn validate(&self) -> Result<()> {
        finite(self.begin, "scan window begin")?;
        finite(self.end, "scan window end")?;
        if self.begin > self.end {
            return Err(invalid("scan window begin exceeds end"));
        }
        validate_meta(&self.metadata)
    }
    /// Whether `mz` lies in the window, both bounds included.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when [`ScanWindow::validate`] fails or
    /// `mz` is not finite.
    pub fn contains(&self, mz: f64) -> Result<bool> {
        self.validate()?;
        finite(mz, "query m/z")?;
        Ok(self.begin <= mz && mz <= self.end)
    }
}

/// Instrument state for one scan: scan mode, polarity and acquisition windows.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InstrumentSettings {
    pub scan_mode: ScanMode,
    pub zoom_scan: bool,
    pub polarity: Polarity,
    pub scan_windows: Vec<ScanWindow>,
    pub metadata: MetaInfo,
}
impl InstrumentSettings {
    /// Check the metadata and every scan window.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the metadata is invalid or any
    /// window fails [`ScanWindow::validate`].
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        for window in &self.scan_windows {
            window.validate()?;
        }
        Ok(())
    }
}

/// Origin of the data: file name and path, size, type and checksum.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SourceFile {
    pub name: String,
    pub path: String,
    /// Megabytes, following the source public float accessor.
    pub size_mb: f32,
    pub file_type: String,
    pub checksum: String,
    pub checksum_type: ChecksumType,
    pub native_id_type: String,
    pub native_id_type_accession: String,
    pub cv_terms: CVTermList,
}
impl SourceFile {
    /// Check the size, the checksum against its declared algorithm, and the CV
    /// terms.
    ///
    /// A nonempty checksum must have the hexadecimal length its
    /// [`ChecksumType`] implies — 40 for SHA-1, 32 for MD5 — and
    /// [`ChecksumType::Unknown`] imposes no shape at all.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative or nonfinite size, a
    /// checksum that does not match its algorithm, or an invalid CV term.
    pub fn validate(&self) -> Result<()> {
        nonnegative(f64::from(self.size_mb), "source file size")?;
        if !self.checksum.is_empty() {
            let length = match self.checksum_type {
                ChecksumType::Sha1 => Some(40),
                ChecksumType::Md5 => Some(32),
                ChecksumType::Unknown => None,
            };
            if length.is_some_and(|length| {
                self.checksum.len() != length
                    || !self.checksum.bytes().all(|b| b.is_ascii_hexdigit())
            }) {
                return Err(invalid("checksum does not match its declared algorithm"));
            }
        }
        self.cv_terms.validate()
    }
}

/// One acquisition that contributed to a combined spectrum.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Acquisition {
    pub identifier: String,
    pub metadata: MetaInfo,
}
impl Acquisition {
    /// Check the metadata.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the metadata is invalid. The
    /// identifier is free text and is not examined.
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)
    }
}

/// The acquisitions behind one spectrum and how they were combined.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AcquisitionInfo {
    pub acquisitions: Vec<Acquisition>,
    pub method_of_combination: String,
    pub metadata: MetaInfo,
}
impl AcquisitionInfo {
    /// Check the metadata and every acquisition.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the metadata is invalid or any
    /// acquisition fails [`Acquisition::validate`].
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        for acquisition in &self.acquisitions {
            acquisition.validate()?;
        }
        Ok(())
    }
}

/// Name, version and CV terms of a piece of software.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Software {
    pub name: String,
    pub version: String,
    pub cv_terms: CVTermList,
}
impl Software {
    /// Check the attached CV terms.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid CV term. Name and version
    /// are free text and are not examined.
    pub fn validate(&self) -> Result<()> {
        self.cv_terms.validate()
    }
}

/// Validated Gregorian wall-clock timestamp, without an inferred time zone.
/// Parses the OpenMS display form YYYY-MM-DD HH:MM:SS (also accepts T separator).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CompletionTime {
    year: u16,
    month: u8,
    day: u8,
    hour: u8,
    minute: u8,
    second: u8,
}
impl FromStr for CompletionTime {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self> {
        if text.len() != 19 || !text.is_ascii() {
            return Err(invalid("invalid completion timestamp"));
        }
        let b = text.as_bytes();
        if b[4] != b'-'
            || b[7] != b'-'
            || !matches!(b[10], b' ' | b'T')
            || b[13] != b':'
            || b[16] != b':'
        {
            return Err(invalid("invalid completion timestamp"));
        }
        let parse = |range: std::ops::Range<usize>| -> Result<u16> {
            if !b[range.clone()].iter().all(u8::is_ascii_digit) {
                return Err(invalid("invalid timestamp digits"));
            }
            text[range]
                .parse()
                .map_err(|_| invalid("invalid timestamp number"))
        };
        let year = parse(0..4)?;
        let month = parse(5..7)? as u8;
        let day = parse(8..10)? as u8;
        let hour = parse(11..13)? as u8;
        let minute = parse(14..16)? as u8;
        let second = parse(17..19)? as u8;
        // chrono decides whether the day exists (proleptic Gregorian). The
        // fixed-width parse above stays, because `NaiveDateTime::parse_from_str`
        // also accepts a sign, unpadded fields and second 60.
        let calendar_day =
            chrono::NaiveDate::from_ymd_opt(year.into(), month.into(), day.into()).is_some();
        if year == 0 || !calendar_day || hour > 23 || minute > 59 || second > 59 {
            return Err(invalid("completion timestamp is outside calendar bounds"));
        }
        Ok(Self {
            year,
            month,
            day,
            hour,
            minute,
            second,
        })
    }
}
impl fmt::Display for CompletionTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            self.year, self.month, self.day, self.hour, self.minute, self.second
        )
    }
}

/// One processing step: the software that ran, what it did and when.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DataProcessing {
    pub software: Software,
    pub actions: BTreeSet<ProcessingAction>,
    /// Complete source date-time, including stored milliseconds. None replaces
    /// the source unset sentinel; formats reject invalid/partial Some values.
    pub completion_time: Option<crate::data_structures::DateTime>,
    pub metadata: MetaInfo,
}
impl DataProcessing {
    /// Check the software record and the metadata.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when [`Software::validate`] fails or the
    /// metadata is invalid. The action set and the completion time are already
    /// typed, so neither needs checking here.
    pub fn validate(&self) -> Result<()> {
        self.software.validate()?;
        validate_meta(&self.metadata)
    }
}

/// Standalone acquisition settings; existing kernel scalar fields are not
/// implicitly synchronized when this object is used separately.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumSettings {
    pub spectrum_type: SpectrumType,
    /// Stored ion-mobility format, the source `im_type_`. Defaults to
    /// [`IonMobilityFormat::Unknown`], as in the source. Read and written
    /// through [`SpectrumSettings::im_format`] and
    /// [`SpectrumSettings::set_im_format`].
    pub ion_mobility_format: IonMobilityFormat,
    /// Stored ion-mobility peak type, the source `im_peak_type_`. Defaults to
    /// [`IonMobilityPeakType::Unknown`], as in the source. Read and written
    /// through [`SpectrumSettings::im_peak_type`] and
    /// [`SpectrumSettings::set_im_peak_type`].
    pub ion_mobility_peak_type: IonMobilityPeakType,
    pub native_id: String,
    pub comment: String,
    pub instrument_settings: InstrumentSettings,
    pub acquisition_info: AcquisitionInfo,
    pub source_file: SourceFile,
    pub precursors: Vec<PrecursorInfo>,
    pub products: Vec<Product>,
    pub data_processing: Vec<DataProcessing>,
    /// Native settings attachment; C++ stores these on MSSpectrum itself.
    pub peptide_identifications: Vec<crate::identification::PeptideIdentification>,
    pub metadata: MetaInfo,
}
impl SpectrumSettings {
    /// Names of the spectrum types, the source
    /// `SpectrumSettings::NamesOfSpectrumType` (`SpectrumSettings.cpp:21`), in
    /// enum order.
    ///
    /// The source array is sized by `SIZE_OF_SPECTRUMTYPE` and so has one entry
    /// per real type; the sentinel itself has no name.
    pub const NAMES_OF_SPECTRUM_TYPE: [&'static str; 3] = ["Unknown", "Centroid", "Profile"];

    /// All spectrum type names known to OpenMS, the source
    /// `getAllNamesOfSpectrumType`, which copies the array into a `StringList`.
    ///
    /// [`SpectrumSettings::NAMES_OF_SPECTRUM_TYPE`] is the same list without the
    /// allocation.
    pub fn all_names_of_spectrum_type() -> Vec<String> {
        Self::NAMES_OF_SPECTRUM_TYPE
            .iter()
            .map(|name| (*name).to_string())
            .collect()
    }

    /// The name of a spectrum type, the source `spectrumTypeToString`.
    ///
    /// Infallible: the source's `@throws Exception::InvalidValue` fires only for
    /// `SIZE_OF_SPECTRUMTYPE`, which [`SpectrumType`] does not have.
    pub const fn spectrum_type_to_string(spectrum_type: SpectrumType) -> &'static str {
        match spectrum_type {
            SpectrumType::Unknown => Self::NAMES_OF_SPECTRUM_TYPE[0],
            SpectrumType::Centroid => Self::NAMES_OF_SPECTRUM_TYPE[1],
            SpectrumType::Profile => Self::NAMES_OF_SPECTRUM_TYPE[2],
        }
    }

    /// Convert an entry of [`SpectrumSettings::NAMES_OF_SPECTRUM_TYPE`] to its
    /// spectrum type, the source `toSpectrumType`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `name` is not one of the names, as
    /// the source throws `Exception::InvalidValue` with the offending string.
    /// Matching is exact and case-sensitive, as the source `std::find`.
    pub fn to_spectrum_type(name: &str) -> Result<SpectrumType> {
        match name {
            "Unknown" => Ok(SpectrumType::Unknown),
            "Centroid" => Ok(SpectrumType::Centroid),
            "Profile" => Ok(SpectrumType::Profile),
            other => Err(invalid(&format!("unknown spectrum type '{other}'"))),
        }
    }

    /// Check every nested record and the metadata.
    ///
    /// The source validates nothing; this is the native guard that keeps a
    /// settings object usable by the transports, and `unify` runs it on both
    /// sides before merging anything.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a nonfinite or out-of-domain value in
    /// the instrument settings, the acquisition info, the source file, any
    /// precursor, product or processing record, any attached peptide
    /// identification, or the metadata.
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        self.instrument_settings.validate()?;
        self.acquisition_info.validate()?;
        self.source_file.validate()?;
        for precursor in &self.precursors {
            precursor.validate()?;
        }
        for product in &self.products {
            product.validate()?;
        }
        for processing in &self.data_processing {
            processing.validate()?;
        }
        for identification in &self.peptide_identifications {
            identification.validate()?;
        }
        Ok(())
    }
    /// The stored ion-mobility format, or
    /// [`IonMobilityFormat::Unknown`] (the default) when none was set.
    ///
    /// Ports `SpectrumSettings::getIMFormat` (`METADATA/SpectrumSettings.h:103`).
    /// The source `@note` still applies: when this is `Unknown`, derive the
    /// format from the data with
    /// [`ImTypes::determine_im_format`](crate::metadata::ImTypes::determine_im_format),
    /// or with
    /// [`ImTypes::determine_im_format_with_stored`](crate::metadata::ImTypes::determine_im_format_with_stored)
    /// to reproduce the source's stored-value short circuit exactly.
    ///
    /// The source `operator==` does not compare `im_type_` or `im_peak_type_`,
    /// so two C++ settings objects differing only in their ion-mobility
    /// annotation compare equal. The derived [`PartialEq`] here compares both
    /// fields; `tests/metadata.rs` asserts that divergence rather than hiding
    /// it, and `unify` leaves both fields local, as the source does.
    pub const fn im_format(&self) -> IonMobilityFormat {
        self.ion_mobility_format
    }
    /// Set the stored ion-mobility format.
    ///
    /// Ports `SpectrumSettings::setIMFormat` (`METADATA/SpectrumSettings.h:97`),
    /// which stores the value unconditionally; there is nothing to validate, so
    /// this cannot fail. Setting it does not touch the spectrum's data:
    /// annotating a format the peaks do not carry is the caller's error, exactly
    /// as in the source.
    pub const fn set_im_format(&mut self, im_format: IonMobilityFormat) {
        self.ion_mobility_format = im_format;
    }
    /// The stored ion-mobility peak type, or
    /// [`IonMobilityPeakType::Unknown`] (the default) when none was set.
    ///
    /// Ports `SpectrumSettings::getIMPeakType`
    /// (`METADATA/SpectrumSettings.h:111`). The source readers use `Unknown` as
    /// "not annotated" and substitute
    /// [`IonMobilityPeakType::Profile`] for ion-mobility data that arrived
    /// without a peak type; the ion-mobility peak picker writes
    /// [`IonMobilityPeakType::Centroid`] together with
    /// [`IonMobilityFormat::PerPeak`] on its output
    /// (`PROCESSING/CENTROIDING/PeakPickerIM.cpp:975-976`).
    pub const fn im_peak_type(&self) -> IonMobilityPeakType {
        self.ion_mobility_peak_type
    }
    /// Set the stored ion-mobility peak type.
    ///
    /// Ports `SpectrumSettings::setIMPeakType`
    /// (`METADATA/SpectrumSettings.h:107`), an unchecked store like
    /// [`SpectrumSettings::set_im_format`].
    pub const fn set_im_peak_type(&mut self, im_peak_type: IonMobilityPeakType) {
        self.ion_mobility_peak_type = im_peak_type;
    }
    /// Source unify: incoming metadata overwrites, comments concatenate without
    /// a delimiter, lists append, conflicting spectrum types become Unknown.
    /// Native ID, instrument/source/acquisition and ion-mobility settings stay local.
    pub fn unify(&mut self, other: &Self) -> Result<()> {
        self.validate()?;
        other.validate()?;
        let mut merged = self.clone();
        merge_meta(
            &mut merged.metadata,
            &other.metadata,
            MetaMergePolicy::Overwrite,
        )?;
        if merged.spectrum_type != other.spectrum_type {
            merged.spectrum_type = SpectrumType::Unknown;
        }
        merged.comment.push_str(&other.comment);
        merged.precursors.extend(other.precursors.iter().cloned());
        merged.products.extend(other.products.iter().cloned());
        merged
            .data_processing
            .extend(other.data_processing.iter().cloned());
        merged
            .peptide_identifications
            .extend(other.peptide_identifications.iter().cloned());
        *self = merged;
        Ok(())
    }
}

/// Acquisition settings of one chromatogram, with a single precursor and
/// product rather than the spectrum's lists.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChromatogramSettings {
    pub chromatogram_type: ChromatogramType,
    pub native_id: String,
    pub comment: String,
    pub instrument_settings: InstrumentSettings,
    pub acquisition_info: AcquisitionInfo,
    pub source_file: SourceFile,
    pub precursor: PrecursorInfo,
    pub product: Product,
    pub data_processing: Vec<DataProcessing>,
    pub metadata: MetaInfo,
}
impl ChromatogramSettings {
    /// Check every nested record and the metadata.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid value in the instrument
    /// settings, the acquisition info, the source file, the precursor, the
    /// product, any processing record or the metadata.
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        self.instrument_settings.validate()?;
        self.acquisition_info.validate()?;
        self.source_file.validate()?;
        self.precursor.validate()?;
        self.product.validate()?;
        for processing in &self.data_processing {
            processing.validate()?;
        }
        Ok(())
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!("{name} must be finite")))
    }
}
fn nonnegative(value: f64, name: &str) -> Result<()> {
    finite(value, name)?;
    if value < 0.0 {
        Err(Error::InvalidValue(format!("{name} must be nonnegative")))
    } else {
        Ok(())
    }
}
fn checked_window(center: f64, lower: f64, upper: f64) -> Result<NumericRange> {
    let min = center - lower;
    let max = center + upper;
    finite(min, "window lower bound")?;
    finite(max, "window upper bound")?;
    Ok(NumericRange { min, max })
}

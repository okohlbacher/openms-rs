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
        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
        pub enum $name { $($(#[$variant_attr])* $variant),+ }
        impl $name {
            pub const ALL: &'static [Self] = &[$(Self::$variant),+];
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
named_enum!(#[derive(Default)] IonMobilityFormat { None => "none", PerPeak => "im_peak", PerSpectrum => "im_spectrum", #[default] Unknown => "unknown" });
named_enum!(#[derive(Default)] IonMobilityPeakType { Profile => "im_profile", Centroid => "im_centroided", #[default] Unknown => "unknown" });
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
        let mass = self.mz * charge - charge * crate::chemistry::PROTON_MASS_U;
        finite(mass, "uncharged mass")?;
        Ok(mass)
    }
}

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
    pub fn validate(&self) -> Result<()> {
        finite(self.mz, "product m/z")?;
        nonnegative(self.isolation_window_lower_offset, "lower isolation offset")?;
        nonnegative(self.isolation_window_upper_offset, "upper isolation offset")?;
        self.cv_terms.validate()
    }
    pub fn isolation_window(&self) -> Result<NumericRange> {
        self.validate()?;
        checked_window(
            self.mz,
            self.isolation_window_lower_offset,
            self.isolation_window_upper_offset,
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScanWindow {
    pub begin: f64,
    pub end: f64,
    pub metadata: MetaInfo,
}
impl ScanWindow {
    pub fn new(begin: f64, end: f64) -> Result<Self> {
        let result = Self {
            begin,
            end,
            ..Self::default()
        };
        result.validate()?;
        Ok(result)
    }
    pub fn validate(&self) -> Result<()> {
        finite(self.begin, "scan window begin")?;
        finite(self.end, "scan window end")?;
        if self.begin > self.end {
            return Err(invalid("scan window begin exceeds end"));
        }
        validate_meta(&self.metadata)
    }
    pub fn contains(&self, mz: f64) -> Result<bool> {
        self.validate()?;
        finite(mz, "query m/z")?;
        Ok(self.begin <= mz && mz <= self.end)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct InstrumentSettings {
    pub scan_mode: ScanMode,
    pub zoom_scan: bool,
    pub polarity: Polarity,
    pub scan_windows: Vec<ScanWindow>,
    pub metadata: MetaInfo,
}
impl InstrumentSettings {
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        for window in &self.scan_windows {
            window.validate()?;
        }
        Ok(())
    }
}

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

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Acquisition {
    pub identifier: String,
    pub metadata: MetaInfo,
}
impl Acquisition {
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AcquisitionInfo {
    pub acquisitions: Vec<Acquisition>,
    pub method_of_combination: String,
    pub metadata: MetaInfo,
}
impl AcquisitionInfo {
    pub fn validate(&self) -> Result<()> {
        validate_meta(&self.metadata)?;
        for acquisition in &self.acquisitions {
            acquisition.validate()?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Software {
    pub name: String,
    pub version: String,
    pub cv_terms: CVTermList,
}
impl Software {
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
        let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
        let days = match month {
            1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
            4 | 6 | 9 | 11 => 30,
            2 => {
                if leap {
                    29
                } else {
                    28
                }
            }
            _ => 0,
        };
        if year == 0 || day == 0 || day > days || hour > 23 || minute > 59 || second > 59 {
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
    pub ion_mobility_format: IonMobilityFormat,
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

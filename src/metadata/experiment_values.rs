// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Owned sample and instrument values. These store the full source scalar
//! domain; validation of physical values and serialization belongs to consumers.

use super::{MetaInfo, Polarity, Software, value::hash_float};
use crate::{Error, Result};
use std::{
    fmt,
    hash::{Hash, Hasher},
    str::FromStr,
};

// Keep the existing metadata enum API: source order, exact names, no SIZE sentinel.
macro_rules! named_enum {
    ($(#[$doc:meta])* $name:ident { $($variant:ident => $label:literal),+ $(,)? }) => {
        $(#[$doc])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $name { #[default] Unknown, $($variant),+ }
        impl $name {
            /// Every valid value in source enum order, without allocating.
            pub const ALL: &'static [Self] = &[Self::Unknown, $(Self::$variant),+];
            /// Exact source spelling, including historical capitalization/typos.
            pub const fn name(self) -> &'static str {
                match self { Self::Unknown => "Unknown", $(Self::$variant => $label),+ }
            }
        }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str(self.name()) }
        }
        impl FromStr for $name {
            type Err = Error;
            fn from_str(value: &str) -> Result<Self> {
                Self::ALL.iter().copied().find(|item| item.name() == value)
                    .ok_or_else(|| Error::InvalidValue(concat!("unknown ", stringify!($name)).into()))
            }
        }
    };
}

named_enum!(/// Aggregation state of the sample.
    SampleState {
    Solid => "solid", Liquid => "liquid", Gas => "gas", Solution => "solution",
    Emulsion => "emulsion", Suspension => "suspension"
});
named_enum!(/// Ion-source inlet types in source order.
    InletType {
    Direct => "Direct", Batch => "Batch", Chromatography => "Chromatography",
    ParticleBeam => "Particle beam", MembraneSeparator => "Membrane sparator",
    OpenSplit => "Open split", JetSeparator => "Jet separator", Septum => "Septum",
    Reservoir => "Reservoir", MovingBelt => "Moving belt", MovingWire => "Moving wire",
    FlowInjectionAnalysis => "Flow injection analysis", ElectrosprayInlet => "Electro spray",
    ThermosprayInlet => "Thermo spray", Infusion => "Infusion",
    ContinuousFlowFastAtomBombardment => "Continuous flow fast atom bombardment",
    InductivelyCoupledPlasma => "Inductively coupled plasma", Membrane => "Membrane inlet",
    Nanospray => "Nanospray inlet"
});
named_enum!(/// Source ionization methods; abbreviations retain the source distinctions.
    IonizationMethod {
    Esi => "Electrospray ionisation", Ei => "Electron ionization", Ci => "Chemical ionisation",
    Fab => "Fast atom bombardment", Tsp => "Thermospray", Ld => "Laser desorption",
    Fd => "Field desorption", Fi => "Flame ionization", Pd => "Plasma desorption",
    Si => "Secondary ion MS", Ti => "Thermal ionization", Api => "Atmospheric pressure ionisation",
    Isi => "ISI", Cid => "Collsion induced decomposition", Cad => "Collsiona activated decomposition",
    Hn => "HN", Apci => "Atmospheric pressure chemical ionization", Appi => "Atmospheric pressure photo ionization",
    Icp => "Inductively coupled plasma", Nesi => "Nano electrospray ionization", Mesi => "Micro electrospray ionization",
    Seldi => "Surface enhanced laser desorption ionization", Send => "Surface enhanced neat desorption",
    Fib => "Fast ion bombardment", Maldi => "Matrix-assisted laser desorption ionization",
    Mpi => "Multiphoton ionization", Di => "Desorption ionization", Fa => "Flowing afterglow",
    Fii => "Field ionization", GdMs => "Glow discharge ionization", Nici => "Negative ion chemical ionization",
    Nrms => "Neutralization reionization mass spectrometry", Pi => "Photoionization", Pyms => "Pyrolysis mass spectrometry",
    Rempi => "Resonance enhanced multiphoton ionization", Ai => "Adiabatic ionization", Asi => "Associative ionization",
    Ad => "Autodetachment", Aui => "Autoionization", Cei => "Charge exchange ionization", Chemi => "Chemi-ionization",
    Dissi => "Dissociative ionization", Lsi => "Liquid secondary ionization", Pei => "Penning ionization",
    Soi => "Soft ionization", Spi => "Spark ionization", Sui => "Surface ionization", Vi => "Vertical ionization",
    ApMaldi => "Atmospheric pressure matrix-assisted laser desorption ionization", Sili => "Desorption/ionization on silicon",
    Saldi => "Surface-assisted laser desorption ionization"
});
named_enum!(/// Mass-analyzer types in source order.
    AnalyzerType {
    Quadrupole => "Quadrupole", PaulIonTrap => "Quadrupole ion trap",
    RadialEjectionLinearIonTrap => "Radial ejection linear ion trap", AxialEjectionLinearIonTrap => "Axial ejection linear ion trap",
    Tof => "Time-of-flight", Sector => "Magnetic sector",
    FourierTransform => "Fourier transform ion cyclotron resonance mass spectrometer",
    IonStorage => "Ion storage", Esa => "Electrostatic energy analyzer", It => "Ion trap",
    Swift => "Stored waveform inverse fourier transform", Cyclotron => "Cyclotron", Orbitrap => "Orbitrap",
    Lit => "Linear ion trap"
});
named_enum!(/// Standard measure used to determine mass resolution.
    ResolutionMethod { Fwhm => "Full width at half max", TenPercentValley => "Ten percent valley", Baseline => "Baseline" });
named_enum!(/// Dependence of the resolution on mass.
    ResolutionType { Constant => "Constant", Proportional => "Proportional" });
named_enum!(/// Direction of the mass scan.
    ScanDirection { Up => "Up", Down => "Down" });
named_enum!(/// Scan law; the source display spelling for Linear is `Linar`.
    ScanLaw { Exponential => "Exponential", Linear => "Linar", Quadratic => "Quadratic" });
named_enum!(/// Reflectron state, including the explicit source None value.
    ReflectronState { On => "On", Off => "Off", None => "None" });
named_enum!(/// Ion detector type, distinct from its acquisition mode.
    DetectorType {
    ElectronMultiplier => "Electron multiplier", Photomultiplier => "Photo multiplier",
    FocalPlaneArray => "Focal plane array", FaradayCup => "Faraday cup",
    ConversionDynodeElectronMultiplier => "Conversion dynode electron multiplier",
    ConversionDynodePhotomultiplier => "Conversion dynode photo multiplier",
    MultiCollector => "Multi-collector", ChannelElectronMultiplier => "Channel electron multiplier",
    Channeltron => "channeltron", DalyDetector => "daly detector", MicrochannelPlateDetector => "microchannel plate detector",
    ArrayDetector => "array detector", ConversionDynode => "conversion dynode", Dynode => "dynode",
    FocalPlaneCollector => "focal plane collector", IonToPhotonDetector => "ion-to-photon detector",
    PointCollector => "point collector", PostaccelerationDetector => "postacceleration detector",
    PhotodiodeArrayDetector => "photodiode array detector", InductiveDetector => "inductive detector",
    ElectronMultiplierTube => "electron multiplier tube"
});
named_enum!(/// Detector electronics acquisition mode, not mzML scan materialization.
    DetectorAcquisitionMode {
    PulseCounting => "Pulse counting", Adc => "Analog-digital converter", Tdc => "Time-digital converter",
    TransientRecorder => "Transient recorder"
});
named_enum!(/// Instrument ion-optics types in source order.
    IonOpticsType {
    MagneticDeflection => "magnetic deflection", DelayedExtraction => "delayed extraction",
    CollisionQuadrupole => "collision quadrupole", SelectedIonFlowTube => "selected ion flow tube",
    TimeLagFocusing => "time lag focusing", Reflectron => "reflectron", EinzelLens => "einzel lens",
    FirstStabilityRegion => "first stability region", FringingField => "fringing field",
    KineticEnergyAnalyzer => "kinetic energy analyzer", StaticField => "static field"
});

/// Sample description with independently owned, ordered subsamples.
/// All floating values are stored unchanged, including signed zero, infinity and
/// NaN. IEEE equality therefore makes a sample containing NaN unequal to itself.
/// Field assignment, recursive Clone/equality/Drop and Vec operations have normal
/// Rust costs; this value type does not impose physical or tree-size restrictions.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Sample {
    pub name: String,
    pub organism: String,
    pub number: String,
    pub comment: String,
    pub state: SampleState,
    /// Grams.
    pub mass: f64,
    /// Milliliters.
    pub volume: f64,
    /// Grams per liter.
    pub concentration: f64,
    pub subsamples: Vec<Sample>,
    pub metadata: MetaInfo,
}

/// One source component. Vector position and the signed instrument order are
/// separate; setting order does not reorder any containing instrument.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IonSource {
    pub inlet_type: InletType,
    pub ionization_method: IonizationMethod,
    pub polarity: Polarity,
    pub order: i32,
    pub metadata: MetaInfo,
}
impl Hash for IonSource {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.order.hash(state);
        self.inlet_type.hash(state);
        self.ionization_method.hash(state);
        (self.polarity as u8).hash(state);
        self.metadata.hash(state);
    }
}

/// Complete analyzer value. Scalars retain their source f64 domain and equality;
/// a consuming format/algorithm is responsible for its own validity constraints.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassAnalyzer {
    pub analyzer_type: AnalyzerType,
    pub resolution_method: ResolutionMethod,
    pub resolution_type: ResolutionType,
    pub scan_direction: ScanDirection,
    pub scan_law: ScanLaw,
    pub reflectron_state: ReflectronState,
    /// Maximum m/z at which peaks can be resolved under the selected method.
    pub resolution: f64,
    /// Parts per million.
    pub accuracy: f64,
    /// Source-documented seconds.
    pub scan_rate: f64,
    /// Seconds per scan.
    pub scan_time: f64,
    /// Meters.
    pub tof_total_path_length: f64,
    /// Precursor isolation width in m/z.
    pub isolation_width: f64,
    pub final_ms_exponent: i32,
    /// Tesla.
    pub magnetic_field_strength: f64,
    pub order: i32,
    pub metadata: MetaInfo,
}
// Source specialization deliberately omits metadata, even though equality uses it.
impl Hash for MassAnalyzer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.analyzer_type.hash(state);
        self.resolution_method.hash(state);
        self.resolution_type.hash(state);
        self.scan_direction.hash(state);
        self.scan_law.hash(state);
        self.reflectron_state.hash(state);
        for value in [
            self.resolution,
            self.accuracy,
            self.scan_rate,
            self.scan_time,
            self.tof_total_path_length,
            self.isolation_width,
            self.magnetic_field_strength,
        ] {
            hash_float(value, state);
        }
        self.final_ms_exponent.hash(state);
        self.order.hash(state);
    }
}

/// Complete detector value. Floating fields retain their original bits; IEEE
/// comparison treats both zeros as equal and any stored NaN as unequal.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IonDetector {
    pub detector_type: DetectorType,
    pub acquisition_mode: DetectorAcquisitionMode,
    pub resolution: f64,
    /// Analog-to-digital converter sampling frequency.
    pub adc_sampling_frequency: f64,
    pub order: i32,
    pub metadata: MetaInfo,
}
// As in the source specialization, metadata is excluded from this hash.
impl Hash for IonDetector {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.detector_type.hash(state);
        self.acquisition_mode.hash(state);
        hash_float(self.resolution, state);
        hash_float(self.adc_sampling_frequency, state);
        self.order.hash(state);
    }
}

/// Complete instrument description with ordered component families and owned
/// software/metadata. Native ownership replaces source getters/setters and copy
/// operations. No component count/order or model-name validation is imposed here.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Instrument {
    pub name: String,
    pub vendor: String,
    pub model: String,
    pub customizations: String,
    pub ion_sources: Vec<IonSource>,
    pub mass_analyzers: Vec<MassAnalyzer>,
    pub ion_detectors: Vec<IonDetector>,
    pub software: Software,
    pub ion_optics: IonOpticsType,
    pub metadata: MetaInfo,
}

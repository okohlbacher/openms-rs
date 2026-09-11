// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Literal source enum/accession routing, independent of dictionary names.
//! Source reader/writer pairs are audited in instrument_enum_pairs.json.

use crate::{Error, Result, metadata::*};

macro_rules! terms {
    ($read:ident, $write:ident, $ty:ident, $unknown:expr, {$($variant:ident => $id:literal),* $(,)?}) => {
        pub(super) fn $read(id: &str) -> Option<$ty> {
            if Some(id) == $unknown { return Some($ty::Unknown); }
            match id { $($id => Some($ty::$variant),)* _ => None }
        }
        pub(super) fn $write(value: $ty) -> Result<Option<&'static str>> {
            if value == $ty::Unknown { return Ok($unknown); }
            match value { $($ty::$variant => Ok(Some($id)),)*
                _ => Err(Error::Unsupported(concat!("mzML cannot preserve ", stringify!($ty), " value").into())) }
        }
    };
}

terms!(ion_optics, write_ion_optics, IonOpticsType, None, {
    DelayedExtraction => "MS:1000246",
    MagneticDeflection => "MS:1000221",
    CollisionQuadrupole => "MS:1000275",
    SelectedIonFlowTube => "MS:1000281",
    TimeLagFocusing => "MS:1000286",
    Reflectron => "MS:1000300",
    EinzelLens => "MS:1000307",
    FirstStabilityRegion => "MS:1000309",
    FringingField => "MS:1000310",
    KineticEnergyAnalyzer => "MS:1000311",
    StaticField => "MS:1000320",
});

terms!(inlet, write_inlet, InletType, None, {
    ContinuousFlowFastAtomBombardment => "MS:1000055",
    Direct => "MS:1000056",
    ElectrosprayInlet => "MS:1000057",
    FlowInjectionAnalysis => "MS:1000058",
    InductivelyCoupledPlasma => "MS:1000059",
    Infusion => "MS:1000060",
    JetSeparator => "MS:1000061",
    MembraneSeparator => "MS:1000062",
    MovingBelt => "MS:1000063",
    MovingWire => "MS:1000064",
    OpenSplit => "MS:1000065",
    ParticleBeam => "MS:1000066",
    Reservoir => "MS:1000067",
    Septum => "MS:1000068",
    ThermosprayInlet => "MS:1000069",
    Batch => "MS:1000248",
    Chromatography => "MS:1000249",
    Membrane => "MS:1000396",
    Nanospray => "MS:1000485",
});

terms!(ionization, write_ionization, IonizationMethod, Some("MS:1000008"), {
    Ci => "MS:1000071",
    Esi => "MS:1000073",
    Fab => "MS:1000074",
    Mpi => "MS:1000227",
    Api => "MS:1000240",
    Di => "MS:1000247",
    Fa => "MS:1000255",
    Fii => "MS:1000258",
    GdMs => "MS:1000259",
    Nici => "MS:1000271",
    Nrms => "MS:1000272",
    Pi => "MS:1000273",
    Pyms => "MS:1000274",
    Rempi => "MS:1000276",
    Ai => "MS:1000380",
    Asi => "MS:1000381",
    Ad => "MS:1000383",
    Aui => "MS:1000384",
    Cei => "MS:1000385",
    Chemi => "MS:1000386",
    Dissi => "MS:1000388",
    Ei => "MS:1000389",
    Lsi => "MS:1000395",
    Pei => "MS:1000399",
    Pd => "MS:1000400",
    Si => "MS:1000402",
    Soi => "MS:1000403",
    Spi => "MS:1000404",
    Sui => "MS:1000406",
    Ti => "MS:1000407",
    Vi => "MS:1000408",
    Fib => "MS:1000446",
    Apci => "MS:1000070",
    ApMaldi => "MS:1000239",
    Appi => "MS:1000382",
    Maldi => "MS:1000075",
    Fd => "MS:1000257",
    Sili => "MS:1000387",
    Ld => "MS:1000393",
    Saldi => "MS:1000405",
    Mesi => "MS:1000397",
    Nesi => "MS:1000398",
    Seldi => "MS:1000278",
    Send => "MS:1000279",
});

terms!(analyzer, write_analyzer, AnalyzerType, Some("MS:1000443"), {
    FourierTransform => "MS:1000079",
    Sector => "MS:1000080",
    Quadrupole => "MS:1000081",
    Tof => "MS:1000084",
    Esa => "MS:1000254",
    It => "MS:1000264",
    Swift => "MS:1000284",
    Cyclotron => "MS:1000288",
    Orbitrap => "MS:1000484",
    AxialEjectionLinearIonTrap => "MS:1000078",
    PaulIonTrap => "MS:1000082",
    RadialEjectionLinearIonTrap => "MS:1000083",
    Lit => "MS:1000291",
});

terms!(reflectron, write_reflectron, ReflectronState, None, {
    Off => "MS:1000105",
    On => "MS:1000106",
});

terms!(detector, write_detector, DetectorType, Some("MS:1000026"), {
    Channeltron => "MS:1000107",
    DalyDetector => "MS:1000110",
    FaradayCup => "MS:1000112",
    MicrochannelPlateDetector => "MS:1000114",
    MultiCollector => "MS:1000115",
    Photomultiplier => "MS:1000116",
    ElectronMultiplier => "MS:1000253",
    ArrayDetector => "MS:1000345",
    ConversionDynode => "MS:1000346",
    Dynode => "MS:1000347",
    FocalPlaneCollector => "MS:1000348",
    IonToPhotonDetector => "MS:1000349",
    PointCollector => "MS:1000350",
    PostaccelerationDetector => "MS:1000351",
    PhotodiodeArrayDetector => "MS:1000621",
    InductiveDetector => "MS:1000624",
    ConversionDynodeElectronMultiplier => "MS:1000108",
    ConversionDynodePhotomultiplier => "MS:1000109",
    ElectronMultiplierTube => "MS:1000111",
    FocalPlaneArray => "MS:1000113",
});

terms!(acquisition, write_acquisition, DetectorAcquisitionMode, None, {
    Adc => "MS:1000117",
    PulseCounting => "MS:1000118",
    Tdc => "MS:1000119",
    TransientRecorder => "MS:1000120",
});

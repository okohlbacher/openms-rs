// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Every `START_SECTION` of `IMTypes_test.cpp` at Core SDK `bc9cc12`, the three
//! ion-mobility sections of `MSSpectrum_test.cpp` that earlier waves could only
//! map, the one `SpectrumSettings_test.cpp` section whose members this work
//! package added, and the native boundaries of `src/metadata/im_types.rs`.
//!
//! Every expectation carrying a section reference is a transcribed C++ literal
//! (tier 3, source review). The finiteness and buffer-gas checks, the ceiling
//! and the name-table agreement are native invariants (tier 4): the source has
//! no analogue for them.

use openms::kernel::DataArray;
use openms::kernel::SpectrumType;
use openms::metadata::{
    DriftTimeUnit, ImTypes, IonMobilityFormat, IonMobilityPeakType, NAMES_OF_DRIFT_TIME_UNIT,
    NAMES_OF_IM_FORMAT, NAMES_OF_IM_PEAK_TYPE, SpectrumSettings, drift_time_unit_to_string,
    im_format_to_string, im_peak_type_to_string, to_drift_time_unit, to_im_format, to_im_peak_type,
};
use openms::{MSExperiment, MSSpectrum};

/// The name `IMDataArrayUtils::setIMUnit` writes for `DriftTimeUnit::VSSC`: the
/// PSI-MS name of `MS:1003008 ! raw inverse reduced ion mobility array`
/// (`IMDataArrayUtils.cpp:28`). `reshapeIMFrameToSingle` goes through that
/// function, so this is the array name `IMTypes_test.cpp`'s `IMwithFDA` fixture
/// carries.
const IM_VSSC_ARRAY: &str = "raw inverse reduced ion mobility array";

/// The generic `UserParam` array name of `IMTypes_test.cpp:193`.
const IM_USER_PARAM_ARRAY: &str = "Ion Mobility";

/// `IMwithDrift`, `IMTypes_test.cpp:93-98`: one drift time of 123.4 for the
/// whole spectrum, in volt-seconds per square centimetre.
fn im_with_drift() -> MSSpectrum {
    let mut spectrum = MSSpectrum::new();
    spectrum.set_drift_time(Some(123.4)).unwrap();
    spectrum.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    spectrum
}

/// `IMwithFDA`, `IMTypes_test.cpp:100-106`: the same spectrum after
/// `IMDataConverter::reshapeIMFrameToSingle`, i.e. the per-spectrum drift time
/// replaced by a per-peak ion-mobility float data array.
///
/// `IMDataConverter` is not ported, so the reshaped spectrum is built directly:
/// one float data array under the name `setIMUnit` writes for the fixture's
/// `VSSC` unit, and no per-spectrum drift time. That is exactly the shape the
/// two `determineIMFormat` overloads examine — `containsIMData()` and
/// `getDriftTime()` — and the substitution is recorded in
/// `docs/IM_TYPES_SUPPORT.md`.
fn im_with_float_array() -> MSSpectrum {
    MSSpectrum {
        float_data_arrays: vec![DataArray::new(IM_VSSC_ARRAY, vec![123.4_f32])],
        ..MSSpectrum::new()
    }
}

/// Construct a `T` through its `Default` implementation, written generically so
/// that clippy's `default_constructed_unit_structs` does not fire on the one
/// type whose `Default` this file exists to exercise.
fn default_of<T: Default>() -> T {
    T::default()
}

// Section 1: IMTypes(). The section constructs one with `new IMTypes` and
// checks the pointer is not null. The Rust counterpart is a zero-sized unit
// struct, so construction cannot fail or allocate; both are asserted, together
// with the `Default` implementation standing in for the default constructor.
#[test]
fn im_types_default_constructor() {
    let constructed: ImTypes = default_of();
    assert_eq!(constructed, ImTypes);
    assert_eq!(std::mem::size_of::<ImTypes>(), 0);
}

// Section 2: ~IMTypes(). `delete e_ptr` in C++. The Rust type has no drop glue
// at all, which is the closest observable property and stronger than "the
// destructor ran without crashing".
#[test]
fn im_types_destructor() {
    assert!(!std::mem::needs_drop::<ImTypes>());
}

// Section 3: DriftTimeUnit toDriftTimeUnit(const std::string&).
#[test]
fn to_drift_time_unit_covers_every_name_and_rejects_others() {
    assert_eq!(to_drift_time_unit("<NONE>").unwrap(), DriftTimeUnit::None);
    for (index, name) in NAMES_OF_DRIFT_TIME_UNIT.iter().enumerate() {
        assert_eq!(to_drift_time_unit(name).unwrap(), DriftTimeUnit::ALL[index]);
    }
    // TEST_EXCEPTION(Exception::InvalidValue, toDriftTimeUnit("haha"))
    assert!(to_drift_time_unit("haha").is_err());
}

// Section 4: const std::string& driftTimeUnitToString(const DriftTimeUnit).
// The section's third assertion passes SIZE_OF_DRIFTTIMEUNIT and expects
// InvalidValue; DriftTimeUnit has no such variant, so there is no value to pass.
#[test]
fn drift_time_unit_to_string_covers_every_variant() {
    assert_eq!(drift_time_unit_to_string(DriftTimeUnit::None), "<NONE>");
    for (index, name) in NAMES_OF_DRIFT_TIME_UNIT.iter().enumerate() {
        assert_eq!(drift_time_unit_to_string(DriftTimeUnit::ALL[index]), *name);
    }
    assert_eq!(DriftTimeUnit::ALL.len(), NAMES_OF_DRIFT_TIME_UNIT.len());
}

// Section 5: IMFormat toIMFormat(const std::string&).
#[test]
fn to_im_format_covers_every_name_and_rejects_others() {
    assert_eq!(to_im_format("none").unwrap(), IonMobilityFormat::None);
    for (index, name) in NAMES_OF_IM_FORMAT.iter().enumerate() {
        assert_eq!(to_im_format(name).unwrap(), IonMobilityFormat::ALL[index]);
    }
    // TEST_EXCEPTION(Exception::InvalidValue, toIMFormat("haha"))
    assert!(to_im_format("haha").is_err());
}

// Section 6: const std::string& imFormatToString(const IMFormat).
#[test]
fn im_format_to_string_covers_every_variant() {
    assert_eq!(im_format_to_string(IonMobilityFormat::None), "none");
    for (index, name) in NAMES_OF_IM_FORMAT.iter().enumerate() {
        assert_eq!(im_format_to_string(IonMobilityFormat::ALL[index]), *name);
    }
    assert_eq!(IonMobilityFormat::ALL.len(), NAMES_OF_IM_FORMAT.len());
}

// Section 7: IMPeakType string conversions, IMTypes_test.cpp:79-89.
#[test]
fn im_peak_type_string_conversions() {
    assert_eq!(
        to_im_peak_type("im_profile").unwrap(),
        IonMobilityPeakType::Profile
    );
    assert_eq!(
        to_im_peak_type("im_centroided").unwrap(),
        IonMobilityPeakType::Centroid
    );
    assert_eq!(
        to_im_peak_type("unknown").unwrap(),
        IonMobilityPeakType::Unknown
    );
    assert_eq!(
        im_peak_type_to_string(IonMobilityPeakType::Profile),
        "im_profile"
    );
    assert_eq!(
        im_peak_type_to_string(IonMobilityPeakType::Centroid),
        "im_centroided"
    );
    assert_eq!(
        im_peak_type_to_string(IonMobilityPeakType::Unknown),
        "unknown"
    );
    // TEST_EXCEPTION(Exception::InvalidValue, toIMPeakType("garbage"))
    assert!(to_im_peak_type("garbage").is_err());
}

// Section 8: static IMFormat determineIMFormat(const MSExperiment&, int).
#[test]
fn determine_im_format_per_ms_level() {
    // empty experiment
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&MSExperiment::new(), 1).unwrap(),
        IonMobilityFormat::None
    );

    let mut plain = MSExperiment::new();
    plain.spectra.push(MSSpectrum::new()); // default MS level = 1
    plain.spectra.push(MSSpectrum::new());
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&plain, 1).unwrap(),
        IonMobilityFormat::None
    );

    let mut with_drift = MSExperiment::new();
    with_drift.spectra.push(MSSpectrum::new());
    with_drift.spectra.push(im_with_drift());
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&with_drift, 1).unwrap(),
        IonMobilityFormat::PerSpectrum
    );

    let mut with_array = MSExperiment::new();
    with_array.spectra.push(MSSpectrum::new());
    with_array.spectra.push(im_with_float_array());
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&with_array, 1).unwrap(),
        IonMobilityFormat::PerPeak
    );

    // MS1 = IM_PEAK, MS2 = IM_SPECTRUM: per-level queries work independently.
    let mut ms1_peak = im_with_float_array();
    ms1_peak.ms_level = 1;
    let mut ms2_drift = im_with_drift();
    ms2_drift.ms_level = 2;
    let mut per_level = MSExperiment::new();
    per_level.spectra.push(ms1_peak);
    per_level.spectra.push(ms2_drift);
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&per_level, 1).unwrap(),
        IonMobilityFormat::PerPeak
    );
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&per_level, 2).unwrap(),
        IonMobilityFormat::PerSpectrum
    );

    // no spectra of the requested level
    let mut only_ms1 = MSExperiment::new();
    only_ms1.spectra.push(im_with_drift());
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&only_ms1, 2).unwrap(),
        IonMobilityFormat::None
    );

    // mixed formats within the same MS level throw
    let mut mixed = MSExperiment::new();
    mixed.spectra.push(im_with_drift());
    mixed.spectra.push(im_with_float_array());
    assert!(ImTypes::determine_im_format_for_ms_level(&mixed, 1).is_err());
}

// Section 9: static IMFormat determineIMFormat(const MSSpectrum&).
#[test]
fn determine_im_format_per_spectrum() {
    assert_eq!(
        ImTypes::determine_im_format(&MSSpectrum::new()),
        IonMobilityFormat::None
    );
    // single IM value for the whole spectrum
    assert_eq!(
        ImTypes::determine_im_format(&im_with_drift()),
        IonMobilityFormat::PerSpectrum
    );
    // IM frame with a float meta-data array
    assert_eq!(
        ImTypes::determine_im_format(&im_with_float_array()),
        IonMobilityFormat::PerPeak
    );
    // setting both is valid (typically concatenated plus some average value)
    let mut both = im_with_float_array();
    both.set_drift_time(Some(123.4)).unwrap();
    assert_eq!(
        ImTypes::determine_im_format(&both),
        IonMobilityFormat::PerPeak
    );
}

// Section 10: determineIMFormat returns IM_PEAK for centroided IM data,
// IMTypes_test.cpp:189-197. The section also calls setIMPeakType, which the
// determination never reads; the array name alone decides.
#[test]
fn determine_im_format_for_centroided_im_data() {
    let spectrum = MSSpectrum {
        float_data_arrays: vec![DataArray::new(IM_USER_PARAM_ARRAY, Vec::<f32>::new())],
        ..MSSpectrum::new()
    };
    let mut settings = SpectrumSettings::default();
    settings.set_im_peak_type(IonMobilityPeakType::Centroid);
    assert_eq!(
        ImTypes::determine_im_format_with_stored(settings.im_format(), &spectrum),
        IonMobilityFormat::PerPeak
    );
}

// Section 11: static double oneOverK0ToCCS(double, double, int, double).
#[test]
fn one_over_k0_to_ccs() {
    // Reserpine [M+H]+ (m/z 609.28, z=1): 1/K0 ~1.196 gives a CCS near the
    // published N2 value of ~245 Angstrom^2. TOLERANCE_ABSOLUTE(0.01).
    let ccs = ImTypes::one_over_k0_to_ccs(1.196, 609.28, 1).unwrap();
    assert!(
        (ccs - 244.9402).abs() < 0.01,
        "CCS {ccs} is not within 0.01 of 244.9402"
    );

    // the charge sign must not matter (|z| is used)
    assert_eq!(
        ImTypes::one_over_k0_to_ccs(0.9, 300.0, -1).unwrap(),
        ImTypes::one_over_k0_to_ccs(0.9, 300.0, 1).unwrap()
    );

    // larger 1/K0 gives a larger CCS (monotonic)
    assert!(
        ImTypes::one_over_k0_to_ccs(1.2, 300.0, 1).unwrap()
            > ImTypes::one_over_k0_to_ccs(0.9, 300.0, 1).unwrap()
    );

    // invalid inputs throw
    assert!(ImTypes::one_over_k0_to_ccs(0.0, 300.0, 1).is_err());
    assert!(ImTypes::one_over_k0_to_ccs(-1.0, 300.0, 1).is_err());
    assert!(ImTypes::one_over_k0_to_ccs(0.9, 0.0, 1).is_err());
    assert!(ImTypes::one_over_k0_to_ccs(0.9, 300.0, 0).is_err());
}

// Section 12: static double ccsToOneOverK0(double, double, int, double).
#[test]
fn ccs_to_one_over_k0_round_trips() {
    // 1/K0 -> CCS -> 1/K0 must recover the original value
    let one_over_k0 = 0.95;
    let ccs = ImTypes::one_over_k0_to_ccs(one_over_k0, 412.5, 1).unwrap();
    let back = ImTypes::ccs_to_one_over_k0(ccs, 412.5, 1).unwrap();
    assert!(
        (back - one_over_k0).abs() < 1e-9,
        "round trip gave {back}, not {one_over_k0}"
    );

    // round trip for a multiply charged ion
    let ok0_2 = 0.62;
    let ccs2 = ImTypes::one_over_k0_to_ccs(ok0_2, 524.3, 2).unwrap();
    let back2 = ImTypes::ccs_to_one_over_k0(ccs2, 524.3, 2).unwrap();
    assert!(
        (back2 - ok0_2).abs() < 1e-9,
        "round trip gave {back2}, not {ok0_2}"
    );

    // invalid inputs throw
    assert!(ImTypes::ccs_to_one_over_k0(0.0, 300.0, 1).is_err());
    assert!(ImTypes::ccs_to_one_over_k0(200.0, 300.0, 0).is_err());
}

// MSSpectrum_test.cpp section 65: void setIMFormat(IMFormat imf),
// MSSpectrum_test.cpp:1453-1461. The source calls it on an MSSpectrum, which
// inherits SpectrumSettings; the Rust MSSpectrum flattens SpectrumSettings and
// carries no copy of the two fields, so the assertions are made on the
// declaring type. The gap is recorded in docs/IM_TYPES_SUPPORT.md.
#[test]
fn set_im_format_round_trips() {
    let mut settings = SpectrumSettings::default();
    settings.set_im_format(IonMobilityFormat::PerPeak);
    assert_eq!(settings.im_format(), IonMobilityFormat::PerPeak);
    settings.set_im_format(IonMobilityFormat::None);
    assert_eq!(settings.im_format(), IonMobilityFormat::None);
}

// MSSpectrum_test.cpp section 66: IMPeakType getIMPeakType() const,
// MSSpectrum_test.cpp:1463-1468.
#[test]
fn get_im_peak_type_defaults_to_unknown() {
    let settings = SpectrumSettings::default();
    assert_eq!(settings.im_peak_type(), IonMobilityPeakType::Unknown);
}

// MSSpectrum_test.cpp section 67: void setIMPeakType(IMPeakType),
// MSSpectrum_test.cpp:1470-1478.
#[test]
fn set_im_peak_type_round_trips() {
    let mut settings = SpectrumSettings::default();
    settings.set_im_peak_type(IonMobilityPeakType::Centroid);
    assert_eq!(settings.im_peak_type(), IonMobilityPeakType::Centroid);
    settings.set_im_peak_type(IonMobilityPeakType::Profile);
    assert_eq!(settings.im_peak_type(), IonMobilityPeakType::Profile);
}

// Native: getIMFormat has no section of its own; the header documents the
// default and the note that an UNKNOWN value should be resolved from the data.
#[test]
fn get_im_format_defaults_to_unknown_and_directs_to_the_determination() {
    let settings = SpectrumSettings::default();
    assert_eq!(settings.im_format(), IonMobilityFormat::Unknown);
    // The note: an UNKNOWN stored format falls through to the data.
    assert_eq!(
        ImTypes::determine_im_format_with_stored(settings.im_format(), &im_with_drift()),
        IonMobilityFormat::PerSpectrum
    );
}

// Native: the stored-format short circuit of determineIMFormat(spec). Any
// stored value other than UNKNOWN wins outright, including NONE on a spectrum
// that does carry ion mobility (IMTypesExperiment.cpp:56-60).
#[test]
fn stored_format_short_circuits_the_determination() {
    let spectrum = im_with_float_array();
    for stored in [
        IonMobilityFormat::None,
        IonMobilityFormat::PerPeak,
        IonMobilityFormat::PerSpectrum,
    ] {
        assert_eq!(
            ImTypes::determine_im_format_with_stored(stored, &spectrum),
            stored
        );
    }
    assert_eq!(
        ImTypes::determine_im_format_with_stored(IonMobilityFormat::Unknown, &spectrum),
        IonMobilityFormat::PerPeak
    );
}

// Native: the three name tables are the source arrays and must stay in step
// with the enums' own labels, which are the port's single definition.
#[test]
fn name_tables_match_the_enum_labels() {
    let units: Vec<&str> = DriftTimeUnit::ALL.iter().map(|u| u.name()).collect();
    assert_eq!(units, NAMES_OF_DRIFT_TIME_UNIT.to_vec());
    let formats: Vec<&str> = IonMobilityFormat::ALL.iter().map(|f| f.name()).collect();
    assert_eq!(formats, NAMES_OF_IM_FORMAT.to_vec());
    let peak_types: Vec<&str> = IonMobilityPeakType::ALL.iter().map(|p| p.name()).collect();
    assert_eq!(peak_types, NAMES_OF_IM_PEAK_TYPE.to_vec());
}

// Native: the unset drift time sentinel and the default buffer gas mass are the
// source constants, and the spectrum sentinel agrees with the kernel accessor.
#[test]
fn constants_match_the_source() {
    assert_eq!(ImTypes::DRIFTTIME_NOT_SET, -1.0);
    assert_eq!(ImTypes::N2_BUFFER_GAS_MASS, 28.0);
    assert_eq!(ImTypes::MAX_SPECTRA, 100_000_000);
    let spectrum = MSSpectrum::new();
    assert_eq!(spectrum.drift_time, ImTypes::DRIFTTIME_NOT_SET);
    assert!(!spectrum.has_drift_time());
}

// Native: a negative MS level matches no spectrum. The source compares with
// std::cmp_not_equal against an unsigned level, so the sign is not wrapped.
#[test]
fn negative_ms_level_matches_nothing() {
    let mut experiment = MSExperiment::new();
    experiment.spectra.push(im_with_drift());
    assert_eq!(
        ImTypes::determine_im_format_for_ms_level(&experiment, -1).unwrap(),
        IonMobilityFormat::None
    );
}

// Native: the conversions reject what the source accepts. A NaN passes every
// `<= 0.0` comparison in C++ and yields NaN; a non-positive buffer gas mass is
// not examined at all there and makes the reduced mass negative.
#[test]
fn conversions_reject_nonfinite_and_bad_buffer_gas() {
    assert!(ImTypes::one_over_k0_to_ccs(f64::NAN, 300.0, 1).is_err());
    assert!(ImTypes::one_over_k0_to_ccs(0.9, f64::NAN, 1).is_err());
    assert!(ImTypes::one_over_k0_to_ccs(f64::INFINITY, 300.0, 1).is_err());
    assert!(ImTypes::ccs_to_one_over_k0(f64::NAN, 300.0, 1).is_err());
    assert!(ImTypes::ccs_to_one_over_k0(f64::INFINITY, 300.0, 1).is_err());
    assert!(
        ImTypes::one_over_k0_to_ccs_with_buffer_gas(0.9, 300.0, 1, 0.0).is_err(),
        "a zero buffer gas mass must be rejected"
    );
    assert!(ImTypes::one_over_k0_to_ccs_with_buffer_gas(0.9, 300.0, 1, -28.0).is_err());
    assert!(ImTypes::ccs_to_one_over_k0_with_buffer_gas(200.0, 300.0, 1, f64::NAN).is_err());
    // An explicit N2 mass is the same call as the default.
    assert_eq!(
        ImTypes::one_over_k0_to_ccs_with_buffer_gas(0.9, 300.0, 1, ImTypes::N2_BUFFER_GAS_MASS)
            .unwrap(),
        ImTypes::one_over_k0_to_ccs(0.9, 300.0, 1).unwrap()
    );
    // A heavier drift gas raises the reduced mass, and the cross section
    // divides by its square root, so the same 1/K0 maps to a smaller CCS.
    assert!(
        ImTypes::one_over_k0_to_ccs_with_buffer_gas(0.9, 300.0, 1, 40.0).unwrap()
            < ImTypes::one_over_k0_to_ccs(0.9, 300.0, 1).unwrap()
    );
}

// Native: i32::MIN is a usable charge. `std::abs(INT_MIN)` is undefined in C++;
// `i32::unsigned_abs` is not, so the port has no undefined corner here.
#[test]
fn extreme_charge_is_defined() {
    let most_negative = ImTypes::one_over_k0_to_ccs(0.9, 300.0, i32::MIN).unwrap();
    let most_positive = ImTypes::one_over_k0_to_ccs(0.9, 300.0, i32::MAX).unwrap();
    assert!(
        most_negative.is_finite() && most_negative > 0.0,
        "CCS for i32::MIN was {most_negative}"
    );
    // |i32::MIN| is one larger than i32::MAX, and the cross section grows with
    // the absolute charge once the reduced mass has saturated at the gas mass.
    assert!(most_negative > most_positive);
}

// SpectrumSettings_test.cpp section 32: static StringList
// getAllNamesOfSpectrumType(), SpectrumSettings_test.cpp:428-433. Ported here
// because this work package added the three spectrum-type conversions the
// METADATA package had left out.
#[test]
fn all_names_of_spectrum_type() {
    let names = SpectrumSettings::all_names_of_spectrum_type();
    assert_eq!(names.len(), SpectrumSettings::NAMES_OF_SPECTRUM_TYPE.len());
    assert_eq!(names[SpectrumType::Centroid as usize], "Centroid");
    assert_eq!(names[SpectrumType::Profile as usize], "Profile");
}

// Native: spectrumTypeToString and toSpectrumType have no class-test section.
// Every name round-trips, and an unknown name is rejected as the source's
// std::find over NamesOfSpectrumType is.
#[test]
fn spectrum_type_string_conversions_round_trip() {
    for (index, name) in SpectrumSettings::NAMES_OF_SPECTRUM_TYPE.iter().enumerate() {
        let spectrum_type = SpectrumSettings::to_spectrum_type(name).unwrap();
        assert_eq!(spectrum_type as usize, index);
        assert_eq!(
            SpectrumSettings::spectrum_type_to_string(spectrum_type),
            *name
        );
    }
    assert_eq!(
        SpectrumSettings::spectrum_type_to_string(SpectrumType::Unknown),
        "Unknown"
    );
    assert!(SpectrumSettings::to_spectrum_type("centroid").is_err());
    assert!(SpectrumSettings::to_spectrum_type("").is_err());
}

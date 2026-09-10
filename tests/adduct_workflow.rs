// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AdductInfo, C13C12_MASSDIFF_U, CoarseIsotopePatternGenerator, CoarseMassMode, ELECTRON_MASS_U,
    EmpiricalFormula,
};
use openms::kernel::DataArray;
use openms::{MSSpectrum, Peak1D, Precursor};

fn sodium_spectrum() -> MSSpectrum {
    let molecule = EmpiricalFormula::parse("C6H12O6").unwrap();
    let adduct = AdductInfo::parse("M+2Na;2+").unwrap();
    assert!(adduct.is_compatible(&molecule));
    let generator =
        CoarseIsotopePatternGenerator::new(Some(4), CoarseMassMode::Approximate).unwrap();
    let neutral = generator.run(&molecule).unwrap();
    let ion_atoms = molecule.checked_add(adduct.empirical_formula()).unwrap();
    let ion = generator.run(&ion_atoms).unwrap();
    assert_eq!(neutral.len(), 4);
    let peaks = neutral
        .peaks()
        .iter()
        .zip(ion.peaks())
        .map(|(neutral, ion)| {
            // Natural sodium has one isotope: adding it changes the mass origin,
            // but preserves this coarse probability vector. Charge removes electrons.
            let mz = adduct.mz(neutral.mass).unwrap();
            assert!((mz - (ion.mass - 2.0 * ELECTRON_MASS_U) / 2.0).abs() < 1e-12);
            assert!((adduct.neutral_mass(mz).unwrap() - neutral.mass).abs() < 1e-12);
            assert!((neutral.probability - ion.probability).abs() < 1e-14);
            Peak1D::new(mz, neutral.probability as f32)
        })
        .collect();
    let mut spectrum = MSSpectrum::from_peaks(peaks);
    spectrum.ms_level = 2;
    spectrum.rt = 42.0;
    spectrum.native_id = "scan=1".into();
    spectrum.precursors.push(Precursor::new(
        adduct.mz(molecule.mono_mass()).unwrap(),
        adduct.charge(),
    ));
    spectrum
        .metadata
        .insert("neutral_formula".into(), "C6H12O6".into());
    spectrum
        .metadata
        .insert("adduct".into(), adduct.name().into());
    spectrum.integer_data_arrays.push(DataArray::new(
        "adduct charge",
        vec![adduct.charge(); spectrum.len()],
    ));
    spectrum.validate().unwrap();
    spectrum
}

#[test]
fn sodium_adduct_masses_match_full_ion_composition_and_half_spacing() {
    let spectrum = sodium_spectrum();
    assert!(spectrum.is_sorted());
    for peaks in spectrum.peaks.windows(2) {
        assert!((peaks[1].mz - peaks[0].mz - C13C12_MASSDIFF_U / 2.0).abs() < 1e-12);
    }
    assert_eq!(spectrum.integer_data_arrays[0].data, [2; 4]);
}

#[cfg(feature = "mzml")]
#[test]
fn adduct_metadata_charge_and_isotope_masses_roundtrip_through_mzml() {
    use openms::MSExperiment;
    use openms::format::mzml::{self, WriteOptions};
    let spectrum = sodium_spectrum();
    let experiment = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    for zlib_compression in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(&mut bytes, &experiment, &WriteOptions { zlib_compression })
            .unwrap();
        let restored = mzml::read(bytes.as_slice()).unwrap();
        assert_eq!(restored.spectra[0], experiment.spectra[0]);
        let adduct = AdductInfo::parse(&restored.spectra[0].metadata["adduct"]).unwrap();
        assert_eq!(adduct.charge(), restored.spectra[0].precursors[0].charge);
    }
}

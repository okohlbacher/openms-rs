// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    CoarseIsotopePatternGenerator, CoarseMassMode, ELECTRON_MASS_U, EmpiricalFormula,
    FineIsotopePatternGenerator, NAFragmentType, NASequence, Ribonucleotide, RibonucleotideDB,
    RibonucleotideRecord,
};
use openms::kernel::DataArray;
use openms::{MSSpectrum, Peak1D, Precursor};

fn sulfur_rna_spectrum() -> MSSpectrum {
    let rna = NASequence::parse("pA[C*]Gp").unwrap();
    // A+C+G = C29H39N13O14. Two linkages add H-2P2O3S; terminal
    // phosphates add H2P2O6. Negative charge removes two natural H atoms.
    let ion_atoms = EmpiricalFormula::parse("C29H37N13O23P4S").unwrap();
    assert_eq!(rna.formula(NAFragmentType::Full, -2).unwrap(), ion_atoms);
    assert_eq!(ion_atoms.charge(), 0);
    let ion_mass = rna.mono_mass(NAFragmentType::Full, -2).unwrap();
    assert_eq!(ion_mass, ion_atoms.mono_mass() + 2.0 * ELECTRON_MASS_U);
    let generator =
        CoarseIsotopePatternGenerator::new(Some(5), CoarseMassMode::Approximate).unwrap();
    let distribution = generator.run(&ion_atoms).unwrap();
    let peaks = distribution
        .peaks()
        .iter()
        .map(|peak| {
            // Formula isotope peaks have atomic masses; add the two retained
            // electrons before converting to the magnitude of negative m/z.
            Peak1D::new(
                (peak.mass + 2.0 * ELECTRON_MASS_U) / 2.0,
                peak.probability as f32,
            )
        })
        .collect();
    let mut spectrum = MSSpectrum::from_peaks(peaks);
    spectrum.native_id = "scan=12".into();
    spectrum.rt = 15.0;
    spectrum.ms_level = 1;
    spectrum.precursors.push(Precursor::new(ion_mass / 2.0, -2));
    spectrum
        .metadata
        .insert("RNA sequence".into(), rna.to_string().into());
    spectrum
        .metadata
        .insert("RNA ion formula".into(), ion_atoms.to_string().into());
    spectrum
        .integer_data_arrays
        .push(DataArray::new("RNA ion charge", vec![-2; spectrum.len()]));
    spectrum.validate().unwrap();
    spectrum
}

#[test]
fn sulfur_rna_composition_produces_negative_ion_isotope_coordinates() {
    let spectrum = sulfur_rna_spectrum();
    assert_eq!(spectrum.len(), 5);
    assert!(spectrum.is_sorted());
    assert_eq!(spectrum.peaks[0].mz, spectrum.precursors[0].mz);
    assert_eq!(spectrum.integer_data_arrays[0].data, [-2; 5]);
    let sequence = NASequence::parse(spectrum.metadata["RNA sequence"].as_str().unwrap()).unwrap();
    let suffix = sequence.suffix(1).unwrap();
    assert_eq!(suffix.to_string(), "*Gp");
    // Cutting after C* transfers the sulfur linkage to the retained 5' end.
    assert_eq!(
        suffix.formula(NAFragmentType::Full, 0).unwrap(),
        EmpiricalFormula::parse("C10H15N5O10P2S").unwrap()
    );
}

#[test]
fn owned_isotope_labeled_rna_keeps_formula_mass_and_registry_identity() {
    let heavy = Ribonucleotide::from_record(RibonucleotideRecord {
        name: "Carbon-13 adenosine".into(),
        code: "A".into(),
        origin: 'A',
        formula: EmpiricalFormula::parse("(13)C10H13N5O4").unwrap(),
        mono_mass: 42.0, // Independent declared field; sequence uses composition.
        average_mass: 43.0,
        ..Default::default()
    })
    .unwrap();
    let guanosine = RibonucleotideDB::global().get("G").unwrap();
    let registry = RibonucleotideDB::from_records(vec![heavy, (*guanosine).clone()]).unwrap();
    let labeled = NASequence::parse_with_registry("AG", &registry).unwrap();
    let ordinary = NASequence::parse("AG").unwrap();
    assert_eq!(labeled.to_string(), ordinary.to_string());
    assert_ne!(labeled, ordinary);
    assert_eq!(
        labeled.checked_string_with_registry(&registry).unwrap(),
        "AG"
    );
    assert!(
        labeled
            .checked_string_with_registry(RibonucleotideDB::global())
            .is_err()
    );
    drop(registry);
    let formula = EmpiricalFormula::parse("(13)C10C10H25N10O11P").unwrap();
    assert_eq!(labeled.formula(NAFragmentType::Full, 0).unwrap(), formula);
    assert_eq!(
        labeled.mono_mass(NAFragmentType::Full, 0).unwrap(),
        formula.mono_mass()
    );
    assert_eq!(labeled.residues()[0].mono_mass(), 42.0);
    // Preserve the source carbon-13 mass; shortening this literal changes its bits.
    #[allow(clippy::excessive_precision)]
    let expected_shift = 10.0 * (13.003_355_000_000_001 - 12.0);
    assert!(
        (labeled.mono_mass(NAFragmentType::Full, 0).unwrap()
            - ordinary.mono_mass(NAFragmentType::Full, 0).unwrap()
            - expected_shift)
            .abs()
            < 1e-10
    );

    let fine = FineIsotopePatternGenerator::default();
    assert_eq!(
        fine.run(&labeled.formula(NAFragmentType::Full, 0).unwrap())
            .unwrap(),
        fine.run(&formula).unwrap()
    );
}

#[cfg(feature = "mzml")]
#[test]
fn rna_isotope_coordinates_and_negative_charge_roundtrip_through_mzml() {
    use openms::MSExperiment;
    use openms::format::mzml::{self, WriteOptions};
    let experiment = MSExperiment {
        spectra: vec![sulfur_rna_spectrum()],
        ..Default::default()
    };
    for zlib_compression in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(&mut bytes, &experiment, &WriteOptions { zlib_compression })
            .unwrap();
        let restored = mzml::read(bytes.as_slice()).unwrap();
        assert_eq!(restored.spectra[0], experiment.spectra[0]);
        let rna = NASequence::parse(
            restored.spectra[0].metadata["RNA sequence"]
                .as_str()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(
            rna.mono_mass(NAFragmentType::Full, -2).unwrap() / 2.0,
            restored.spectra[0].precursors[0].mz
        );
    }
}

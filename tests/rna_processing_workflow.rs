// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::MSSpectrum;
use openms::chemistry::{
    ModifiedNASequenceGenerator, NASequence, NucleicAcidSpectrumGenerator, RNaseDigestion,
    RibonucleotideDB,
};

fn generator() -> NucleicAcidSpectrumGenerator {
    NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        ..Default::default()
    }
}
fn ion(spectrum: &MSSpectrum, name: &str) -> f64 {
    let index = spectrum.string_data_arrays[0]
        .data
        .iter()
        .position(|n| n == name)
        .unwrap();
    spectrum.peaks[index].mz
}

#[test]
fn installed_ribose_and_base_modifications_have_distinct_digestion_effects() {
    let database = RibonucleotideDB::global();
    let input = NASequence::parse("AGGCA").unwrap();
    let digestion = RNaseDigestion::default();
    let modifications = ModifiedNASequenceGenerator::default();
    let mut ribose = input.clone();
    modifications
        .apply_fixed_modifications(&[database.get("Gm").unwrap()], &mut ribose)
        .unwrap();
    assert_eq!(ribose.to_string(), "A[Gm][Gm]CA");
    assert_eq!(digestion.digest(&ribose).unwrap(), [ribose.clone()]);
    let mut base = input;
    modifications
        .apply_fixed_modifications(&[database.get("m1G").unwrap()], &mut base)
        .unwrap();
    let products = digestion.digest_with_positions(&base).unwrap();
    assert_eq!(
        products
            .iter()
            .map(|p| p.sequence.to_string())
            .collect::<Vec<_>>(),
        ["A[m1G]p", "[m1G]p", "CA"]
    );
    assert_eq!(
        products
            .iter()
            .map(|p| (p.start, p.end))
            .collect::<Vec<_>>(),
        [(0, 2), (2, 3), (3, 5)]
    );
    for product in products {
        let spectrum = generator().generate(&product.sequence, -1, -1).unwrap();
        spectrum.validate().unwrap();
        assert_eq!(
            spectrum.integer_data_arrays[0].data,
            vec![-1; spectrum.len()]
        );
    }
}

fn processed_spectra() -> Vec<MSSpectrum> {
    let database = RibonucleotideDB::global();
    let source = NASequence::parse("AGUC").unwrap();
    let digest = RNaseDigestion::default()
        .digest_with_positions(&source)
        .unwrap();
    assert_eq!(digest[1].sequence.to_string(), "UC");
    let alternatives = [database.get("s4U").unwrap(), database.get("m3U").unwrap()];
    let variants = ModifiedNASequenceGenerator::default()
        .variable_modifications(&alternatives, &digest[1].sequence, 1, true)
        .unwrap();
    assert_eq!(
        variants.iter().map(ToString::to_string).collect::<Vec<_>>(),
        ["UC", "[s4U]C", "[m3U]C"]
    );
    variants
        .iter()
        .enumerate()
        .map(|(i, variant)| {
            let mut spectrum = generator().generate(variant, -1, -1).unwrap();
            spectrum.ms_level = 2;
            spectrum.native_id = format!("scan={}", i + 1);
            spectrum
                .metadata
                .insert("RNA sequence".into(), variant.to_string());
            spectrum
                .metadata
                .insert("RNA parent sequence".into(), source.to_string());
            spectrum
                .metadata
                .insert("RNA parent start".into(), digest[1].start.to_string());
            spectrum
                .metadata
                .insert("RNA parent end exclusive".into(), digest[1].end.to_string());
            spectrum
        })
        .collect()
}

#[test]
fn digest_variants_shift_declared_fragment_masses_and_keep_unaffected_ions() {
    let spectra = processed_spectra();
    let database = RibonucleotideDB::global();
    let unmodified = database.get("U").unwrap().mono_mass();
    for (i, code) in [(1, "s4U"), (2, "m3U")] {
        let delta = database.get(code).unwrap().mono_mass() - unmodified;
        assert!((ion(&spectra[i], "b1") - ion(&spectra[0], "b1") - delta).abs() < 1e-10);
        assert_eq!(ion(&spectra[i], "y1"), ion(&spectra[0], "y1"));
        assert_eq!(spectra[i].integer_data_arrays[0].data, [-1, -1]);
    }
}

#[cfg(feature = "mzml")]
#[test]
fn annotated_rna_products_roundtrip_with_coordinates_and_both_compression_modes() {
    use openms::{
        MSExperiment,
        format::mzml::{self, WriteOptions},
    };
    let experiment = MSExperiment {
        spectra: processed_spectra(),
        ..Default::default()
    };
    for zlib_compression in [false, true] {
        let mut bytes = Vec::new();
        mzml::write_with_options(&mut bytes, &experiment, &WriteOptions { zlib_compression })
            .unwrap();
        let restored = mzml::read(bytes.as_slice()).unwrap();
        assert_eq!(restored.spectra, experiment.spectra);
        for spectrum in restored.spectra {
            let sequence = NASequence::parse(&spectrum.metadata["RNA sequence"]).unwrap();
            let regenerated = generator().generate(&sequence, -1, -1).unwrap();
            assert_eq!(regenerated.peaks, spectrum.peaks);
            assert_eq!(regenerated.string_data_arrays, spectrum.string_data_arrays);
            assert_eq!(
                regenerated.integer_data_arrays,
                spectrum.integer_data_arrays
            );
        }
    }
}

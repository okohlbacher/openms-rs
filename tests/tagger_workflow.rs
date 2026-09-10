// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::analysis::peptide_indexing::PeptideIndexing;
use openms::chemistry::{
    AASequence, DecoyGenerator, ProteaseDigestion, Tagger, TaggerOptions, TheoreticalIonSeries,
    TheoreticalSpectrumGenerator,
};
use openms::comparison::Tolerance;
use openms::{MSSpectrum, Peak1D};

fn observed() -> MSSpectrum {
    // Rounded source residue increments define measured inputs independently of
    // Tagger's mass table and the native theoretical-spectrum generator.
    let mut position = 150.0;
    let mut peaks = vec![Peak1D::new(position, 10.0)];
    for mass in [97.0527, 129.0426, 97.0527, 101.0477] {
        position += mass;
        peaks.push(Peak1D::new(position, 20.0));
    }
    let mut spectrum = MSSpectrum::from_peaks(peaks);
    spectrum.ms_level = 2;
    spectrum.rt = 12.0;
    spectrum.native_id = "scan=8".into();
    spectrum
}

fn tagger() -> Tagger {
    let mut options = TaggerOptions::new(2, Tolerance::Absolute(0.02));
    options.max_tag_length = 4;
    Tagger::new(options).unwrap()
}

#[test]
fn measured_tags_distinguish_target_from_reversed_decoy_by_exact_substring_matching() {
    let spectrum = observed();
    let before = spectrum.clone();
    let tags = tagger().get_spectrum_tags(&spectrum).unwrap();
    assert_eq!(tags, ["EP", "EPT", "PE", "PEP", "PEPT", "PT"]);
    assert_eq!(spectrum, before);
    let longest = tags.iter().max_by_key(|tag| tag.len()).unwrap();
    let target = AASequence::parse("MPEPTIDER").unwrap();
    let decoy = DecoyGenerator::with_seed(4711)
        .reverse_protein(&target)
        .unwrap();
    let matcher = PeptideIndexing::default();
    assert_eq!(matcher.find_matches(longest, target.as_str()).unwrap(), [1]);
    assert!(
        matcher
            .find_matches(longest, decoy.as_str())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn digested_modified_peptide_gaps_recover_parent_letters_with_fixed_and_variable_tables() {
    let protein = AASequence::parse("(Acetyl)AC(Carbamidomethyl)M(Oxidation)KPEP").unwrap();
    let digestion = ProteaseDigestion {
        enzyme: openms::chemistry::Protease::TrypsinP,
        ..Default::default()
    };
    let products = digestion.digest(&protein).unwrap();
    let peptide = &products[0].sequence;
    assert_eq!(peptide.as_str(), "ACMK");
    assert!(peptide.n_terminal_modification().is_some());
    let generator = TheoreticalSpectrumGenerator {
        ion_series: vec![TheoreticalIonSeries::B],
        add_first_prefix_ion: true,
        add_metainfo: true,
        ..Default::default()
    };
    let spectrum = generator.generate(peptide, 1, 1, None).unwrap();
    assert_eq!(spectrum.len(), 3);
    let mut options = TaggerOptions::new(2, Tolerance::Ppm(10.0));
    options.max_tag_length = 2;
    assert!(
        Tagger::new(options.clone())
            .unwrap()
            .get_spectrum_tags(&spectrum)
            .unwrap()
            .is_empty()
    );
    options.fixed_mods.push("Carbamidomethyl (C)".into());
    options.variable_mods.push("Oxidation (M)".into());
    let tags = Tagger::new(options)
        .unwrap()
        .get_spectrum_tags(&spectrum)
        .unwrap();
    // The N-terminal acetyl mass cancels between prefix peaks. Output retains
    // parent letters rather than serializing their modification names.
    assert_eq!(tags, ["CM"]);
    assert_eq!(
        PeptideIndexing::default()
            .find_matches(&tags[0], protein.as_str())
            .unwrap(),
        [1]
    );
}

#[cfg(feature = "mzml")]
#[test]
fn centroided_mzml_roundtrip_keeps_the_same_tag_set_in_both_compression_modes() {
    use openms::MSExperiment;
    use openms::format::mzml::{self, WriteOptions};
    let spectrum = observed();
    let expected = tagger().get_spectrum_tags(&spectrum).unwrap();
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
        assert_eq!(
            tagger().get_spectrum_tags(&restored.spectra[0]).unwrap(),
            expected
        );
    }
}

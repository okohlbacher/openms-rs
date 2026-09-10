// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{
    AASequence, EmpiricalFormula, ModificationRecord, ModificationsDB, PROTON_MASS_U,
    ProteaseDigestion, ResidueModification, SpectrumAnnotator, TheoreticalSpectrumGenerator,
    ion_naming,
};
use openms::comparison::{SpectrumAlignment, Tolerance};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::kernel::DataArray;
use openms::metadata::MetaValue;
use openms::{MSSpectrum, Peak1D, Precursor};

fn generator() -> TheoreticalSpectrumGenerator {
    TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    }
}
fn alignment() -> SpectrumAlignment {
    SpectrumAlignment {
        tolerance: Tolerance::Absolute(1e-5),
        ..Default::default()
    }
}
fn modified_input() -> (PeptideIdentification, MSSpectrum) {
    let protein = AASequence::parse("(Acetyl)AC(Carbamidomethyl)M(Oxidation)KAGHIK").unwrap();
    let products = ProteaseDigestion::default().digest(&protein).unwrap();
    assert_eq!((products[0].start, products[0].end), (0, 4));
    let peptide = products[0].sequence.clone();
    let hit = PeptideHit::new(1.0, 1, 2, peptide.clone()).unwrap();
    // Independent chemical formulas for acetylated/carbamidomethylated AC b2,
    // and intact terminal lysine y1. No theoretical spectrum generates these inputs.
    let b2 = EmpiricalFormula::parse("C10H15N3O4S").unwrap().mono_mass() + PROTON_MASS_U;
    let y1 = EmpiricalFormula::parse("C6H14N2O2").unwrap().mono_mass() + PROTON_MASS_U;
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(b2, 50.0),
        Peak1D::new(999.0, 0.5),
        Peak1D::new(y1, 20.0),
    ]);
    spectrum.ms_level = 2;
    spectrum.rt = 23.0;
    spectrum.native_id = "scan=7".into();
    spectrum
        .precursors
        .push(Precursor::new(peptide.mz(2).unwrap(), 2));
    spectrum
        .metadata
        .insert("sample".into(), "digestion workflow".into());
    spectrum
        .float_data_arrays
        .push(DataArray::new("old float", vec![3.0; 3]));
    spectrum
        .integer_data_arrays
        .push(DataArray::new("old integer", vec![9; 3]));
    spectrum
        .string_data_arrays
        .push(DataArray::new("old string", vec!["previous".into(); 3]));
    let identification = PeptideIdentification {
        identifier: "annotation".into(),
        hits: vec![hit],
        rt: Some(23.0),
        mz: Some(peptide.mz(2).unwrap()),
        ..Default::default()
    };
    (identification, spectrum)
}

fn annotate() -> (PeptideIdentification, MSSpectrum) {
    let (mut id, mut spectrum) = modified_input();
    let annotator = SpectrumAnnotator::default();
    annotator
        .add_ion_match_statistics(&mut id, &mut spectrum, &generator(), &alignment())
        .unwrap();
    annotator
        .add_peak_annotations(&mut id.hits[0], &spectrum, &generator(), &alignment(), true)
        .unwrap();
    (id, spectrum)
}

#[test]
fn digested_modified_peptide_annotates_independent_fragment_formulas_and_retains_acquisition() {
    let (original_id, original) = modified_input();
    let (id, spectrum) = annotate();
    assert_eq!(id.hits[0].sequence, original_id.hits[0].sequence);
    assert_eq!(spectrum.precursors, original.precursors);
    assert_eq!(spectrum.rt, original.rt);
    assert_eq!(spectrum.native_id, original.native_id);
    assert_eq!(spectrum.metadata["sample"], original.metadata["sample"]);
    assert_eq!(spectrum.string_data_arrays.len(), 1);
    assert_eq!(spectrum.float_data_arrays[0].name, "IonMatchError");
    assert_eq!(spectrum.integer_data_arrays[0].name, "Charges");
    assert_eq!(spectrum.string_data_arrays[0].data, ["y1+", "b2+", ""]);
    assert_eq!(spectrum.integer_data_arrays[0].data, [1, 1, 0]);
    assert_eq!(
        id.hits[0].metadata["matched_ion_number"],
        MetaValue::from(2_i64)
    );
    assert_eq!(
        id.hits[0].metadata["matched_intensity"],
        MetaValue::try_from(70.0).unwrap()
    );
    assert_eq!(id.hits[0].metadata["peak_number"], MetaValue::from(3_i64));
    assert_eq!(
        id.hits[0].metadata["precursor_in_ms2"],
        MetaValue::from(0_i64)
    );
    assert_eq!(spectrum.metadata["fragment_mass_tolerance_ppm"], "0");
    for (peak, annotation) in spectrum.peaks.iter().zip(&id.hits[0].peak_annotations) {
        assert_eq!(annotation.mz, peak.mz);
        assert_eq!(annotation.intensity, f64::from(peak.intensity));
        assert_eq!(
            ion_naming::with_charge(&annotation.annotation, annotation.charge).unwrap(),
            annotation.annotation
        );
    }
    assert_eq!(
        ion_naming::ordinal_from_name(&id.hits[0].peak_annotations[1].annotation),
        2
    );
}

#[test]
fn equal_sequence_text_with_distinct_owned_chemistry_produces_distinct_matches() {
    let peptide = |formula: &str| {
        let formula = EmpiricalFormula::parse(formula).unwrap();
        let record = ResidueModification::from_record(ModificationRecord {
            name: "Lab".into(),
            full_name: "Laboratory annotation".into(),
            origin: Some('M'),
            diff_mono_mass: formula.mono_mass(),
            diff_average_mass: formula.average_mass(),
            diff_formula: formula,
            ..Default::default()
        })
        .unwrap();
        let database = ModificationsDB::from_records(vec![record]).unwrap();
        AASequence::parse_with_registry("AM(Lab)A", &database).unwrap()
    };
    let first = peptide("O");
    let second = peptide("O2");
    assert_eq!(first.to_string(), second.to_string());
    assert_ne!(first, second);
    let mz = EmpiricalFormula::parse("C8H14N2O3S").unwrap().mono_mass() + PROTON_MASS_U;
    let observed = MSSpectrum::from_peaks(vec![Peak1D::new(mz, 10.0)]);
    let mut first = PeptideHit::new(1.0, 1, 2, first).unwrap();
    let mut second = PeptideHit::new(1.0, 1, 2, second).unwrap();
    for hit in [&mut first, &mut second] {
        SpectrumAnnotator::default()
            .add_peak_annotations(hit, &observed, &generator(), &alignment(), false)
            .unwrap();
    }
    assert_eq!(first.peak_annotations.len(), 1);
    assert_eq!(first.peak_annotations[0].annotation, "b2+");
    assert!(second.peak_annotations.is_empty());
}

#[cfg(feature = "idxml")]
#[test]
fn modified_chemistry_peak_annotations_and_statistics_roundtrip_through_idxml() {
    use openms::format::idxml::{self, IdXmlDocument};
    use openms::identification::ProteinIdentification;
    let (mut id, _) = annotate();
    let pi = openms::chemistry::IsoelectricPoint::default()
        .compute_pi(&id.hits[0].sequence)
        .unwrap();
    id.hits[0]
        .metadata
        .insert("Lehninger:pI".into(), MetaValue::try_from(pi).unwrap());
    let document = IdXmlDocument {
        protein_identifications: vec![ProteinIdentification {
            identifier: "annotation".into(),
            date_time: Some("2026-09-10T12:00:00".into()),
            search_engine: "annotation demonstration".into(),
            ..Default::default()
        }],
        peptide_identifications: vec![id],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    idxml::write(&mut bytes, &document).unwrap();
    assert_eq!(idxml::read(bytes.as_slice()).unwrap(), document);
}

#[cfg(feature = "mzml")]
#[test]
fn measured_annotation_arrays_and_precursor_survive_both_mzml_compression_modes() {
    use openms::MSExperiment;
    use openms::format::mzml::{self, WriteOptions};
    let (_, spectrum) = annotate();
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
    }
}

use openms::metadata::{MetaInfo, MetaValue, Unit};
use openms::processing::{Normalizer, SpectrumFilter};
use openms::{MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
fn metadata() -> MetaInfo {
    [
        ("empty".into(), MetaValue::default()),
        ("string".into(), "42".into()),
        ("integer".into(), 42_i64.into()),
        (
            "float".into(),
            MetaValue::try_from(2.5)
                .unwrap()
                .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
                .unwrap(),
        ),
        (
            "strings".into(),
            vec![String::from("a"), String::from("b")].into(),
        ),
        ("integers".into(), vec![1_i64, 2].into()),
        ("floats".into(), MetaValue::try_from(vec![1., 2.]).unwrap()),
    ]
    .into()
}
#[test]
fn all_seven_types_and_units_survive_record_clone_and_batched_processing() {
    let s = MSSpectrum {
        metadata: metadata(),
        peaks: vec![Peak1D::new(2., 4.), Peak1D::new(1., 2.)],
        ..Default::default()
    };
    let c = MSChromatogram {
        metadata: metadata(),
        ..Default::default()
    };
    let mut e = MSExperiment {
        spectra: vec![s.clone()],
        chromatograms: vec![c.clone()],
        ..Default::default()
    };
    Normalizer::default().filter_experiment(&mut e).unwrap();
    assert_eq!(e.spectra[0].metadata, s.metadata);
    assert_eq!(e.chromatograms[0], c);
    assert_eq!(e.spectra[0].peaks[0].intensity, 1.);
    e.spectra[0]
        .metadata
        .insert("string".into(), "changed".into());
    assert_eq!(s.metadata["string"].as_str().unwrap(), "42");
}
#[test]
fn flat_formats_preserve_old_strings_and_reject_new_unrepresentable_types_atomically() {
    use openms::format::{dta, mgf};
    use std::io::Cursor;
    let e = mgf::read(Cursor::new(
        "BEGIN IONS\nPEPMASS=500\nCHARGE=2+\nCUSTOM=42\n100 2\nEND IONS\n",
    ))
    .unwrap();
    assert_eq!(e.spectra[0].metadata["CUSTOM"].as_str().unwrap(), "42");
    let mut out = Vec::new();
    mgf::write(&mut out, &e).unwrap();
    assert_eq!(
        mgf::read(Cursor::new(out)).unwrap().spectra[0].metadata,
        e.spectra[0].metadata
    );
    for value in [
        42_i64.into(),
        MetaValue::try_from(42.).unwrap(),
        MetaValue::default(),
        vec![1_i64].into(),
        MetaValue::from("text")
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    ] {
        let mut e = e.clone();
        e.spectra[0].metadata.insert("CUSTOM".into(), value);
        let mut out = b"unchanged".to_vec();
        assert!(mgf::write(&mut out, &e).is_err());
        assert_eq!(out, b"unchanged");
    }
    let mut out = b"unchanged".to_vec();
    assert!(dta::write(&mut out, &e.spectra[0]).is_err());
    assert_eq!(out, b"unchanged");
}

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::PROTON_MASS_U;
use openms::format::{dta, fasta, mgf};
use openms::{Error, MSExperiment, MSSpectrum, Peak1D, Precursor};

#[test]
fn dta_upstream_fixture_and_precursor_conversion() {
    let spec = dta::read(&include_bytes!("data/Transformers_tests.dta")[..]).unwrap();
    assert_eq!(spec.len(), 121);
    assert_eq!(spec.ms_level, 2);
    assert!(
        (spec.precursors[0].mz - ((739.771 - PROTON_MASS_U) / 2.0 + PROTON_MASS_U)).abs() < 1e-12
    );
    assert_eq!(spec.peaks[0], Peak1D::new(104.0541, 3.5));
    let mut text = Vec::new();
    dta::write(&mut text, &spec).unwrap();
    let copy = dta::read(text.as_slice()).unwrap();
    assert_eq!(spec, copy);
}

#[test]
fn dta_exact_and_legacy_mass_conventions() {
    for charge in [0, 1, 2, 3, -2] {
        let spec = MSSpectrum {
            ms_level: 2,
            precursors: vec![Precursor::new(501.25, charge)],
            ..Default::default()
        };
        let mut bytes = Vec::new();
        dta::write(&mut bytes, &spec).unwrap();
        let copy = dta::read(bytes.as_slice()).unwrap();
        assert!((copy.precursors[0].mz - 501.25).abs() < 1e-12);
    }
    let spec = MSSpectrum {
        precursors: vec![Precursor::new(500.0, 2)],
        ..Default::default()
    };
    let mut bytes = Vec::new();
    dta::write_with_convention(&mut bytes, &spec, dta::MassConvention::LegacyOpenMS).unwrap();
    assert_eq!(String::from_utf8(bytes).unwrap(), "999 2\n");
}

#[test]
fn dta_rejects_corrupt_data() {
    for data in [
        "",
        "bad 2",
        "500 2 extra",
        "500 x",
        "500 2\n100 NaN",
        "500 2\n100 1e50",
        "500 2\n100",
        "inf 2",
    ] {
        assert!(dta::read(data.as_bytes()).is_err(), "accepted {data}");
    }
    let error = dta::read(b"500 2\n100 3\ninvalid\n".as_slice()).unwrap_err();
    assert!(matches!(error, Error::Parse { line: 3, .. }));
}

#[test]
fn fasta_streaming_whitespace_and_roundtrip() {
    let input = "\u{feff}>sp|P1| protein one\r\nPEP TIDE\r\nK\r\n\r\n>P2\nGAS\n";
    let entries = fasta::read(input.as_bytes()).unwrap();
    assert_eq!(entries.len(), 2);
    assert_eq!(entries[0].identifier, "sp|P1|");
    assert_eq!(entries[0].description, "protein one");
    assert_eq!(entries[0].sequence, "PEPTIDEK");
    let mut output = Vec::new();
    fasta::write(&mut output, &entries).unwrap();
    assert_eq!(fasta::read(output.as_slice()).unwrap(), entries);
    let mut iter = fasta::FastaReader::new(input.as_bytes());
    assert_eq!(iter.next().unwrap().unwrap(), entries[0]);
    assert_eq!(iter.next().unwrap().unwrap(), entries[1]);
    assert!(iter.next().is_none());
}

#[test]
fn fasta_rejects_malformed_entries_and_stops_after_error() {
    for input in ["PEPTIDE", ">\nABC", ">empty\n", ">p", ">p description"] {
        assert!(fasta::read(input.as_bytes()).is_err(), "accepted {input}");
    }
    let mut reader = fasta::FastaReader::new(b">p\n\xff\n>q\nABC".as_slice());
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}

#[test]
fn fasta_writer_wraps_and_rejects_header_injection_before_writing() {
    let entry = fasta::FASTAEntry {
        identifier: "p".into(),
        description: "".into(),
        sequence: "A".repeat(161),
    };
    let mut output = Vec::new();
    fasta::write(&mut output, std::slice::from_ref(&entry)).unwrap();
    assert_eq!(
        String::from_utf8(output)
            .unwrap()
            .lines()
            .map(str::len)
            .collect::<Vec<_>>(),
        [3, 80, 80, 1]
    );
    let bad = fasta::FASTAEntry {
        description: "injected\n>evil".into(),
        ..entry.clone()
    };
    let mut output = Vec::new();
    assert!(fasta::write(&mut output, &[entry, bad]).is_err());
    assert!(output.is_empty());
}

const MGF: &str = "# comment\nCOM=global\nCHARGE=2+\nBEGIN IONS\nTITLE=peptide spectrum\nPEPMASS=500.25 1000\nRTINSECONDS=30.5\nSCANS=42\n100.1 20\n200.2 40\nEND IONS\n\nBEGIN IONS\nCHARGE=1-\nPEPMASS=600\nMSLEVEL=3\n50 2\nEND IONS\n";

#[test]
fn mgf_streaming_global_parameters_and_charge() {
    let exp = mgf::read(MGF.as_bytes()).unwrap();
    assert_eq!(exp.spectra.len(), 2);
    let spec = &exp.spectra[0];
    assert_eq!(
        spec.precursors[0],
        Precursor {
            mz: 500.25,
            intensity: 1000.0,
            charge: 2,
            ..Precursor::default()
        }
    );
    assert_eq!(spec.rt, 30.5);
    assert_eq!(spec.name, "peptide spectrum");
    assert_eq!(spec.native_id, "index=0");
    assert_eq!(spec.metadata["SCANS"], "42");
    assert_eq!(spec.metadata["COM"], "global");
    assert_eq!(exp.spectra[1].precursors[0].charge, -1);
    assert_eq!(exp.spectra[1].ms_level, 3);
    let mut output = Vec::new();
    mgf::write(&mut output, &exp).unwrap();
    assert_eq!(mgf::read(output.as_slice()).unwrap(), exp);
}

#[test]
fn mgf_rejects_truncation_ambiguous_charge_and_nonfinite_peaks() {
    for text in [
        "BEGIN IONS\n100 20",
        "END IONS",
        "100 20",
        "BEGIN IONS\nBEGIN IONS",
        "BEGIN IONS\nCHARGE=2+ and 3+\nEND IONS",
        "BEGIN IONS\nCHARGE=-2-\nEND IONS",
        "BEGIN IONS\n100 NaN\nEND IONS",
        "BEGIN IONS\n100 1 annotation\nEND IONS",
        "BEGIN IONS\nTITLE=a\nTITLE=b\nEND IONS",
        "BEGIN IONS\nMSLEVEL=0\nEND IONS",
    ] {
        assert!(mgf::read(text.as_bytes()).is_err(), "accepted {text}");
    }
}

#[test]
fn mgf_writer_semantic_error_precedes_any_output() {
    let good = mgf::read(MGF.as_bytes()).unwrap().spectra.remove(0);
    let mut bad = good.clone();
    bad.name = "injected\nEND IONS".into();
    let exp = MSExperiment {
        spectra: vec![good, bad],
        ..Default::default()
    };
    let mut output = Vec::new();
    assert!(mgf::write(&mut output, &exp).is_err());
    assert!(output.is_empty());
}

#[test]
fn io_failures_propagate() {
    struct Broken;
    impl std::io::Write for Broken {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("test failure"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    assert!(matches!(
        dta::write(Broken, &MSSpectrum::default()),
        Err(Error::Io(_))
    ));
}

#[test]
fn text_writers_propagate_final_flush_errors() {
    struct FlushFailure;
    impl std::io::Write for FlushFailure {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("flush failed"))
        }
    }
    assert!(matches!(
        dta::write(FlushFailure, &MSSpectrum::default()),
        Err(Error::Io(_))
    ));
    assert!(matches!(fasta::write(FlushFailure, &[]), Err(Error::Io(_))));
    assert!(matches!(
        mgf::write(FlushFailure, &MSExperiment::default()),
        Err(Error::Io(_))
    ));
}

#[test]
fn short_malformed_ascii_inputs_never_panic() {
    // Bounded deterministic corpus exercising parser state transitions.
    let alphabet = *b">\n 0-=A\t";
    for a in alphabet {
        for b in alphabet {
            for c in alphabet {
                let bytes = [a, b, c];
                let _ = dta::read(bytes.as_slice());
                let _ = fasta::read(bytes.as_slice());
                let _ = mgf::read(bytes.as_slice());
            }
        }
    }
}

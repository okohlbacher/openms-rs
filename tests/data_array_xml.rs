// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(any(
    feature = "mzml",
    feature = "consensusxml",
    feature = "idxml",
    feature = "featurexml"
))]
//! Transport regressions for newly represented DataArray descriptions.
use openms::kernel::DataArray;
use openms::metadata::DataProcessing;
use std::io::{self, Write};
use std::sync::Arc;
#[derive(Default)]
struct NoOutput {
    writes: usize,
    flushes: usize,
}
impl Write for NoOutput {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        self.writes += 1;
        panic!("unsupported array metadata reached output")
    }
    fn flush(&mut self) -> io::Result<()> {
        self.flushes += 1;
        panic!("unsupported array metadata reached flush")
    }
}
fn annotate<T>(array: &mut DataArray<T>, processing: bool) {
    if processing {
        array
            .data_processing
            .push(Arc::new(DataProcessing::default()));
    } else {
        array
            .metadata
            .insert("description".into(), "scientific ownership".into());
    }
}
fn unsupported<T>(result: openms::Result<T>, output: &NoOutput) {
    assert!(
        matches!(result,Err(openms::Error::Unsupported(s)) if s.contains("array description metadata"))
    );
    assert_eq!((output.writes, output.flushes), (0, 0));
}

#[cfg(feature = "mzml")]
#[test]
fn mzml_preserves_all_description_kinds_on_spectra_and_chromatograms() {
    use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
    for chromatogram in [false, true] {
        for processing in [false, true] {
            for empty in [false, true] {
                for kind in 0..3 {
                    let mut exp = MSExperiment::new();
                    exp.spectra
                        .push(MSSpectrum::from_peaks(vec![Peak1D::new(1., 2.)]));
                    exp.chromatograms.push(MSChromatogram::from_peaks(vec![
                        ChromatogramPeak::new(1., 2.),
                    ]));
                    let (floats, integers, strings) = if chromatogram {
                        let c = &mut exp.chromatograms[0];
                        (
                            &mut c.float_data_arrays,
                            &mut c.integer_data_arrays,
                            &mut c.string_data_arrays,
                        )
                    } else {
                        let s = &mut exp.spectra[0];
                        (
                            &mut s.float_data_arrays,
                            &mut s.integer_data_arrays,
                            &mut s.string_data_arrays,
                        )
                    };
                    match kind {
                        0 => {
                            let mut a =
                                DataArray::new("values", if empty { vec![] } else { vec![1_f32] });
                            annotate(&mut a, processing);
                            floats.push(a);
                        }
                        1 => {
                            let mut a =
                                DataArray::new("values", if empty { vec![] } else { vec![1_i32] });
                            annotate(&mut a, processing);
                            integers.push(a);
                        }
                        _ => {
                            let mut a = DataArray::new(
                                "values",
                                if empty { vec![] } else { vec!["a".to_string()] },
                            );
                            annotate(&mut a, processing);
                            strings.push(a);
                        }
                    }
                    let before = exp.clone();
                    let mut output = Vec::new();
                    openms::format::mzml::write(&mut output, &exp).unwrap();
                    let read = openms::format::mzml::read(std::io::Cursor::new(output)).unwrap();
                    let (actual, expected) = if chromatogram {
                        let a = &read.chromatograms[0];
                        let b = &exp.chromatograms[0];
                        (
                            (
                                &a.float_data_arrays,
                                &a.integer_data_arrays,
                                &a.string_data_arrays,
                            ),
                            (
                                &b.float_data_arrays,
                                &b.integer_data_arrays,
                                &b.string_data_arrays,
                            ),
                        )
                    } else {
                        let a = &read.spectra[0];
                        let b = &exp.spectra[0];
                        (
                            (
                                &a.float_data_arrays,
                                &a.integer_data_arrays,
                                &a.string_data_arrays,
                            ),
                            (
                                &b.float_data_arrays,
                                &b.integer_data_arrays,
                                &b.string_data_arrays,
                            ),
                        )
                    };
                    assert_eq!(actual, expected);
                    assert_eq!(exp, before);
                }
            }
        }
    }
}

#[cfg(any(feature = "consensusxml", feature = "idxml", feature = "featurexml"))]
fn run(
    kind: usize,
    processing: bool,
    empty: bool,
) -> openms::identification::ProteinIdentification {
    use openms::identification::{ProteinGroup, ProteinIdentification};
    let mut run = ProteinIdentification {
        identifier: "run".into(),
        ..Default::default()
    };
    let mut group = ProteinGroup {
        probability: 1.,
        accessions: vec!["protein".into()],
        ..Default::default()
    };
    match kind {
        0 => {
            let mut a = DataArray::new("psm_count", if empty { vec![] } else { vec![0_f32] });
            annotate(&mut a, processing);
            group.float_data_arrays.push(a);
        }
        1 => {
            let mut a = DataArray::new("psm_count", if empty { vec![] } else { vec![0_i32] });
            annotate(&mut a, processing);
            group.integer_data_arrays.push(a);
        }
        _ => {
            let mut a = DataArray::new("label", if empty { vec![] } else { vec!["a".to_string()] });
            annotate(&mut a, processing);
            group.string_data_arrays.push(a);
        }
    }
    if kind == 2 {
        run.protein_groups.push(group);
    } else {
        run.indistinguishable_groups.push(group);
    }
    run
}

#[cfg(feature = "consensusxml")]
#[test]
fn consensus_quantities_reject_description_before_empty_or_zero_array_omission() {
    for processing in [false, true] {
        for empty in [false, true] {
            for kind in 0..3 {
                let mut map = openms::kernel::ConsensusMap::default();
                map.protein_identifications
                    .push(run(kind, processing, empty));
                let before = map.clone();
                let mut output = NoOutput::default();
                unsupported(
                    openms::format::consensusxml::write(&mut output, &map),
                    &output,
                );
                assert_eq!(map, before);
            }
        }
    }
}

#[cfg(feature = "idxml")]
#[test]
fn idxml_group_preflight_rejects_description_before_rendering_or_metadata_clones() {
    for kind in 0..3 {
        let document = openms::format::idxml::IdXmlDocument {
            protein_identifications: vec![run(kind, true, true)],
            ..Default::default()
        };
        let mut output = NoOutput::default();
        unsupported(
            openms::format::idxml::write(&mut output, &document),
            &output,
        );
    }
}

#[cfg(feature = "featurexml")]
#[test]
fn featurexml_group_preflight_rejects_description_before_rendering_or_metadata_clones() {
    for kind in 0..3 {
        let mut map = openms::kernel::FeatureMap::default();
        map.protein_identifications.push(run(kind, false, true));
        let mut output = NoOutput::default();
        assert!(
            matches!(openms::format::featurexml::write(&mut output, &map),
            Err(openms::Error::Unsupported(message)) if message.contains("structured protein groups"))
        );
        assert_eq!((output.writes, output.flushes), (0, 0));
    }
}

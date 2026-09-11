#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::format::{
    mzml::{self, MzMLCounts, ReadOptions},
    peak_options::PeakFileOptions,
};
use openms::kernel::NumericRange;
use std::io::{BufRead, BufReader, Cursor, Read, Write};

fn cv(id: &str, value: &str) -> String {
    format!("<cvParam accession=\"MS:{id}\" name=\"count fixture\" value=\"{value}\"/>")
}
fn scan(value: &str) -> String {
    format!("<scanList><scan>{}</scan></scanList>", cv("1000016", value))
}
fn spectrum(body: &str) -> String {
    format!("<spectrum id=\"s\" defaultArrayLength=\"0\">{body}</spectrum>")
}
fn lists(body: &str) -> String {
    format!(
        "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\">{body}</spectrumList><chromatogramList count=\"0\" defaultDataProcessingRef=\"dp\"/>"
    )
}
fn document(header: &str, body: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">{header}<run defaultInstrumentConfigurationRef=\"ic\">{body}</run></mzML>"
    )
}
fn count(xml: &str, opt: &PeakFileOptions) -> openms::Result<MzMLCounts> {
    mzml::read_size_with_options(Cursor::new(xml), opt, &ReadOptions::default())
}
fn filtered() -> PeakFileOptions {
    let mut opt = PeakFileOptions::default();
    opt.add_ms_level(1).unwrap();
    opt
}
fn one(body: &str, opt: &PeakFileOptions) -> usize {
    count(&document("", &lists(&spectrum(body))), opt)
        .unwrap()
        .spectra
}

#[test]
fn all_five_literal_source_count_pairs_on_unchanged_source_fixture() {
    let xml = include_str!("data/mzml_load_source_original.mzML");
    let mut opt = PeakFileOptions::default();
    let mut observed = vec![count(xml, &opt).unwrap()];
    opt.add_ms_level(2).unwrap();
    observed.push(count(xml, &opt).unwrap());
    opt.add_ms_level(1).unwrap();
    observed.push(count(xml, &opt).unwrap());
    opt.clear_ms_levels();
    opt.set_rt_range(NumericRange {
        min: 5.0,
        max: 5.25,
    });
    observed.push(count(xml, &opt).unwrap());
    opt.add_ms_level(1).unwrap();
    observed.push(count(xml, &opt).unwrap());
    let expected: Vec<_> = include_str!("data/mzml_counts_source_goldens.tsv")
        .lines()
        .skip(1)
        .map(|row| {
            let columns: Vec<_> = row.split('\t').collect();
            MzMLCounts {
                spectra: columns[1].parse().unwrap(),
                chromatograms: columns[2].parse().unwrap(),
            }
        })
        .collect();
    assert_eq!(observed, expected);
}

#[test]
fn raw_counts_trust_headers_and_stop_before_unconsumed_malformed_tail() {
    let prefix = "<spectrumList count=\"12\" defaultDataProcessingRef=\"dp\"><spectrum nonsense=\"yes\"><binary>not base64!</binary></spectrum></spectrumList><chromatogramList count=\"34\" defaultDataProcessingRef=\"dp\">";
    let xml = document("", &format!("{prefix}<broken"));
    assert_eq!(
        mzml::read_size(Cursor::new(&xml)).unwrap(),
        MzMLCounts {
            spectra: 12,
            chromatograms: 34
        }
    );
    let only = document(
        "",
        "<spectrumList count=\"12\" defaultDataProcessingRef=\"dp\"><binary>junk!</binary></spectrumList>",
    );
    assert_eq!(
        mzml::read_size(Cursor::new(&only)).unwrap(),
        MzMLCounts {
            spectra: 12,
            chromatograms: 0
        }
    );
    assert!(mzml::read_size(Cursor::new(only.replace("</run>", "<broken</run>"))).is_err());
    let reverse = document(
        "",
        "<chromatogramList count=\"3\" defaultDataProcessingRef=\"dp\"/><spectrumList count=\"2\" defaultDataProcessingRef=\"dp\"><broken",
    );
    assert_eq!(
        mzml::read_size(Cursor::new(reverse)).unwrap(),
        MzMLCounts {
            spectra: 2,
            chromatograms: 3
        }
    );
    assert_eq!(
        mzml::read_size(Cursor::new(document("", ""))).unwrap(),
        MzMLCounts::default()
    );
}

#[test]
fn source_signed_count_sentinel_and_checked_numeric_boundaries() {
    let negative = document(
        "",
        "<spectrumList count=\"-2\" defaultDataProcessingRef=\"dp\"/><chromatogramList count=\"3\" defaultDataProcessingRef=\"dp\"><broken",
    );
    assert_eq!(
        mzml::read_size(Cursor::new(negative)).unwrap(),
        MzMLCounts {
            spectra: 0,
            chromatograms: 3
        }
    );
    let sentinel = document(
        "",
        "<spectrumList count=\"-1\" defaultDataProcessingRef=\"dp\"/><chromatogramList count=\"3\" defaultDataProcessingRef=\"dp\"><broken",
    );
    assert!(mzml::read_size(Cursor::new(sentinel)).is_err());
    for n in ["2147483648", "-2147483649", "++1", "1.2"] {
        assert!(
            mzml::read_size(Cursor::new(document(
                "",
                &format!("<spectrumList count=\"{n}\" defaultDataProcessingRef=\"dp\"/>")
            )))
            .is_err()
        );
    }
}

#[test]
fn metadata_only_stops_before_count_but_after_required_reference() {
    let mut opt = PeakFileOptions::default();
    opt.metadata_only = true;
    let xml = document(
        "",
        "<spectrumList count=\"not a number\" defaultDataProcessingRef=\"dp\"><broken",
    );
    assert_eq!(count(&xml, &opt).unwrap(), MzMLCounts::default());
    assert!(count(&xml.replace(" defaultDataProcessingRef=\"dp\"", ""), &opt).is_err());
}

#[test]
fn count_dispatch_ignores_peak_and_write_options_but_preserves_rt_order() {
    let mut opt = PeakFileOptions::default();
    opt.fill_data = false;
    opt.skip_xml_checks = true;
    opt.always_append_data = true;
    opt.max_data_pool_size = 0;
    opt.set_mz_range(NumericRange {
        min: 100.0,
        max: 101.0,
    });
    opt.set_intensity_range(NumericRange {
        min: 10.0,
        max: 11.0,
    });
    let xml = document(
        "",
        "<spectrumList count=\"9\" defaultDataProcessingRef=\"dp\"><broken",
    );
    // Only one list requires a well-formed remainder, but not record fields.
    assert!(count(&xml, &opt).is_err());
    assert_eq!(count(&document("", &lists("")), &opt).unwrap().spectra, 1);
    let opt = filtered();
    assert_eq!(
        one(&format!("{}{}", cv("1000511", "2"), scan("5")), &opt),
        0
    );
    assert_eq!(
        one(&format!("{}{}", scan("5"), cv("1000511", "2")), &opt),
        1
    );
    assert_eq!(one(&scan("5"), &opt), 1); // Missing MS-level is not rejection.
    assert_eq!(one(&cv("1000511", "1"), &opt), 0); // Missing RT is not counted.
    assert_eq!(
        one(&format!("<scan>{}</scan>", cv("1000826", "5")), &opt),
        0
    );
}

#[test]
fn rt_half_open_minutes_and_checked_nonfinite_values() {
    let mut opt = PeakFileOptions::default();
    opt.set_rt_range(NumericRange {
        min: 60.0,
        max: 120.0,
    });
    assert_eq!(one(&scan("60"), &opt), 1);
    assert_eq!(one(&scan("120"), &opt), 0);
    assert_eq!(
        one(
            "<scan><cvParam accession=\"MS:1000016\" name=\"rt\" value=\"1\" unitAccession=\"UO:0000031\"/></scan>",
            &opt
        ),
        1
    );
    assert_eq!(
        one(
            "<scan><cvParam accession=\"MS:1000016\" name=\"rt\" value=\"60\" unitAccession=\"unknown\"/></scan>",
            &opt
        ),
        1
    );
    for rt in ["NaN", "inf", "1e999", "1e-999"] {
        assert!(count(&document("", &lists(&spectrum(&scan(rt)))), &opt).is_err());
    }
    let overflow = "<scan><cvParam accession=\"MS:1000016\" name=\"rt\" value=\"1e308\" unitAccession=\"UO:0000031\"/></scan>";
    assert!(count(&document("", &lists(&spectrum(overflow))), &opt).is_err());
}

fn group(terms: &str) -> String {
    format!(
        "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"g\">{terms}</referenceableParamGroup></referenceableParamGroupList>"
    )
}
#[test]
fn source_reference_expansion_can_count_more_than_one_rt_per_record() {
    let opt = filtered();
    let terms = format!("{}{}", cv("1000016", "5"), cv("1000016", "6"));
    let header = group(&terms);
    let scans = "<scan><referenceableParamGroupRef ref=\"g\"/><referenceableParamGroupRef ref=\"missing_after_skip\"/></scan>";
    assert_eq!(
        count(&document(&header, &lists(&spectrum(scans))), &opt)
            .unwrap()
            .spectra,
        2
    );
    assert_eq!(one(&format!("<scan>{terms}</scan>"), &opt), 1);
    let mut opt = PeakFileOptions::default();
    opt.set_rt_range(NumericRange { min: 5.5, max: 7.0 });
    assert_eq!(
        count(&document(&header, &lists(&spectrum(scans))), &opt)
            .unwrap()
            .spectra,
        1
    );
    assert_eq!(one(&format!("<scan>{terms}</scan>"), &opt), 0);
}

#[test]
fn parameter_references_validate_consumed_uses_and_forward_header_references() {
    let opt = filtered();
    let header = format!(
        "<fileDescription><fileContent><referenceableParamGroupRef ref=\"g\"/></fileContent></fileDescription>{}",
        group(&cv("1000016", "5"))
    );
    let body = lists(&spectrum(
        "<scan><referenceableParamGroupRef ref=\"g\"/></scan>",
    ));
    assert_eq!(count(&document(&header, &body), &opt).unwrap().spectra, 1);
    assert!(count(&document("", &body), &opt).is_err());
    assert!(
        count(
            &document(&header.replace("ref=\"g\"", "ref=\"missing\""), &body),
            &opt
        )
        .is_err()
    );
    assert!(
        count(
            &document(&group("").replace("count=\"1\"", "count=\"2\""), ""),
            &opt
        )
        .is_err()
    );
    let duplicate = "<referenceableParamGroupList count=\"2\"><referenceableParamGroup id=\"g\"/><referenceableParamGroup id=\"g\"/></referenceableParamGroupList>";
    assert!(count(&document(duplicate, ""), &opt).is_err());
    let recursive = group("<referenceableParamGroupRef ref=\"g\"/>");
    assert!(count(&document(&recursive, ""), &opt).is_err());
}

fn precursor(target: &str, selected: &[&str]) -> String {
    format!(
        "<precursorList><precursor><isolationWindow>{}</isolationWindow><selectedIonList>{}</selectedIonList></precursor></precursorList>",
        cv("1000827", target),
        selected.iter().fold(String::new(), |mut text, n| {
            text.push_str("<selectedIon>");
            text.push_str(&cv("1000744", n));
            text.push_str("</selectedIon>");
            text
        })
    )
}
#[test]
fn precursor_mode_equality_first_ion_and_event_order_are_source_derived() {
    let mut opt = PeakFileOptions::default();
    opt.set_precursor_mz_range(NumericRange {
        min: 200.0,
        max: 300.0,
    });
    assert_eq!(
        one(
            &format!("{}{}", precursor("100", &["100"]), scan("5")),
            &opt
        ),
        1
    );
    assert_eq!(
        one(
            &format!("{}{}", precursor("100", &["300", "250"]), scan("5")),
            &opt
        ),
        0
    );
    assert_eq!(
        one(
            &format!("{}{}", precursor("100", &["250", "300"]), scan("5")),
            &opt
        ),
        1
    );
    assert_eq!(
        one(
            &format!("{}{}", scan("5"), precursor("100", &["300"])),
            &opt
        ),
        1
    );
    opt.precursor_mz_selected_ion = false;
    assert_eq!(
        one(
            &format!("{}{}", precursor("100", &["250"]), scan("5")),
            &opt
        ),
        0
    );
    assert_eq!(
        one(
            &format!("{}{}", precursor("250", &["300"]), scan("5")),
            &opt
        ),
        1
    );
}

#[test]
fn approved_skipping_and_missing_rt_corrections_do_not_decode_binary() {
    let body = format!(
        "<spectrumList count=\"2\" defaultDataProcessingRef=\"dp\">{}{}</spectrumList><chromatogramList count=\"1\" defaultDataProcessingRef=\"dp\"><chromatogram id=\"c\" defaultArrayLength=\"0\"><binary>invalid!</binary></chromatogram></chromatogramList>",
        spectrum(&scan("5")),
        spectrum(
            "<binaryDataArrayList count=\"nonsense\"><binary>invalid numpress!</binary></binaryDataArrayList>"
        )
    );
    let xml = document("", &body);
    let mut opt = filtered();
    assert_eq!(
        count(&xml, &opt).unwrap(),
        MzMLCounts {
            spectra: 1,
            chromatograms: 1
        }
    );
    opt.skip_chromatograms = true;
    assert_eq!(
        count(&xml, &opt).unwrap(),
        MzMLCounts {
            spectra: 1,
            chromatograms: 0
        }
    );
    opt.clear_ms_levels();
    assert_eq!(
        count(&xml, &opt).unwrap(),
        MzMLCounts {
            spectra: 2,
            chromatograms: 0
        }
    );
    assert!(count(&xml.replace("id=\"c\" ", ""), &filtered()).is_err());
}

struct Tiny<R> {
    inner: R,
    width: usize,
}
impl<R: BufRead> Read for Tiny<R> {
    fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
        let n = {
            let a = self.fill_buf()?;
            let n = a.len().min(b.len());
            b[..n].copy_from_slice(&a[..n]);
            n
        };
        self.consume(n);
        Ok(n)
    }
}
impl<R: BufRead> BufRead for Tiny<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let b = self.inner.fill_buf()?;
        Ok(&b[..b.len().min(self.width)])
    }
    fn consume(&mut self, n: usize) {
        self.inner.consume(n)
    }
}
#[test]
fn bounded_text_discard_handles_large_arrays_and_arbitrary_utf8_chunks() {
    let prefix = "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run><spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><binary>";
    let suffix = "</binary></spectrumList></run></mzML>";
    let stream = Cursor::new(prefix)
        .chain(std::io::repeat(b'A').take(3 * 1024 * 1024))
        .chain(Cursor::new(suffix));
    let limits = ReadOptions {
        max_param_bytes: 128 * 1024,
        max_array_bytes: 0,
        max_total_arrays: 0,
        max_total_peaks: 0,
        ..Default::default()
    };
    assert_eq!(
        mzml::read_size_with_options(
            BufReader::with_capacity(64, stream),
            &PeakFileOptions::default(),
            &limits
        )
        .unwrap()
        .spectra,
        1
    );
    let text = document(
        "",
        "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><binary>α😀]]&gt; &amp; &#945; <![CDATA[α & >]]></binary></spectrumList>",
    );
    for width in 1..=7 {
        assert_eq!(
            mzml::read_size(Tiny {
                inner: Cursor::new(&text),
                width
            })
            .unwrap()
            .spectra,
            1
        );
    }
    let mut input = Cursor::new(document("", &lists("")));
    assert_eq!(mzml::read_size(&mut input).unwrap().spectra, 1);
    assert!((input.position() as usize) < input.get_ref().len());
}

#[test]
fn text_legality_namespaces_and_encoding_are_checked_before_the_stop() {
    let raw = document(
        "",
        "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><binary>BODY</binary></spectrumList>",
    );
    for bad in ["]]>", "\u{0}", "&unknown;", "&#0;", "&amp"] {
        for width in 1..=4 {
            assert!(
                mzml::read_size(Tiny {
                    inner: Cursor::new(raw.replace("BODY", bad)),
                    width
                })
                .is_err(),
                "{bad:?}"
            );
        }
    }
    let mut invalid_utf8 = raw.as_bytes().to_vec();
    let at = raw.find("BODY").unwrap();
    invalid_utf8[at] = 0xff;
    assert!(mzml::read_size(Cursor::new(invalid_utf8)).is_err());
    for xml in [
        format!("x{raw}"),
        format!("{raw}x"),
        raw.replace("http://psi.hupo.org/ms/mzml", "wrong"),
        raw.replace("</binary>", "</wrong>"),
    ] {
        assert!(mzml::read_size(Cursor::new(xml)).is_err());
    }
    let ascii = format!(
        "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>{}",
        raw.replace("BODY", "α")
    );
    assert!(mzml::read_size(Cursor::new(ascii)).is_err());
    assert!(
        mzml::read_size(Cursor::new(format!(
            "<?xml version=\"1.0\" encoding=\"UTF-16\"?>{raw}"
        )))
        .is_err()
    );
    assert!(mzml::read_size(Cursor::new(format!("<!DOCTYPE mzML>{raw}"))).is_err());
}

#[test]
fn xml_record_parameter_and_event_limits_are_enforced() {
    let xml = document("", &lists(&spectrum(&scan("5"))));
    for limits in [
        ReadOptions {
            max_xml_bytes: 10,
            ..Default::default()
        },
        ReadOptions {
            max_records: 0,
            ..Default::default()
        },
        ReadOptions {
            max_total_params: 0,
            ..Default::default()
        },
        ReadOptions {
            max_param_bytes: 0,
            ..Default::default()
        },
    ] {
        assert!(mzml::read_size_with_options(Cursor::new(&xml), &filtered(), &limits).is_err());
    }
    let huge = "x".repeat(mzml::MAX_COUNT_EVENT_BYTES + 1);
    for body in [
        format!("<!--{huge}-->"),
        format!("<x value=\"{huge}\"/>"),
        format!("<![CDATA[{huge}]]>"),
        format!("&{huge};"),
    ] {
        assert!(mzml::read_size(Cursor::new(document("", &body))).is_err());
    }
    let nested = format!(
        "{}{}",
        "<x>".repeat(mzml::MAX_COUNT_XML_DEPTH),
        "</x>".repeat(mzml::MAX_COUNT_XML_DEPTH)
    );
    assert!(mzml::read_size(Cursor::new(document("", &nested))).is_err());
    let header = group(&format!("{}{}", cv("1000016", "5"), cv("1000016", "6")));
    let xml = document(
        &header,
        &lists(&spectrum(
            "<scan><referenceableParamGroupRef ref=\"g\"/></scan>",
        )),
    );
    let limits = ReadOptions {
        max_records: 1,
        ..Default::default()
    };
    assert!(mzml::read_size_with_options(Cursor::new(xml), &filtered(), &limits).is_err());
}

struct Temp(std::path::PathBuf);
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn all_path_entrypoints_use_shared_magic_and_compression_policy() {
    let dir = Temp(std::env::temp_dir().join(format!("openms-counts-{}", std::process::id())));
    std::fs::create_dir_all(&dir.0).unwrap();
    let xml = document("", &lists(&spectrum(&scan("5"))));
    let plain = dir.0.join("plain.gz");
    std::fs::write(&plain, &xml).unwrap();
    assert_eq!(mzml::load_size(&plain).unwrap().spectra, 1);
    assert_eq!(
        mzml::load_size_with_options(&plain, &filtered(), &ReadOptions::default())
            .unwrap()
            .spectra,
        1
    );
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(xml.as_bytes()).unwrap();
    let gzip = dir.0.join("gzip.unrelated");
    std::fs::write(&gzip, gz.finish().unwrap()).unwrap();
    assert_eq!(mzml::load_size(&gzip).unwrap().spectra, 1);
    let mut bz = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bz.write_all(xml.as_bytes()).unwrap();
    let bzip = dir.0.join("bzip.xml");
    std::fs::write(&bzip, bz.finish().unwrap()).unwrap();
    assert_eq!(mzml::load_size(&bzip).unwrap().spectra, 1);
    let limits = ReadOptions {
        max_xml_bytes: 20,
        ..Default::default()
    };
    assert!(mzml::load_size_with_options(&gzip, &filtered(), &limits).is_err());
    let zip = dir.0.join("zip.mzML");
    std::fs::write(&zip, b"PK\x03\x04junk").unwrap();
    assert!(matches!(
        mzml::load_size(zip),
        Err(openms::Error::Unsupported(_))
    ));
    assert!(mzml::load_size(dir.0.join("missing")).is_err());
}

#[test]
fn bom_declaration_entity_and_markup_state_survive_one_byte_chunks() {
    let xml = format!(
        "\u{feff}<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>{}",
        document(
            "",
            "<spectrumList count=\"2\" defaultDataProcessingRef=\"dp\"><binary>start &amp; α after &lt; end<child/>after child</binary></spectrumList>"
        )
    );
    for width in 1..=9 {
        assert_eq!(
            mzml::read_size(Tiny {
                inner: Cursor::new(&xml),
                width
            })
            .unwrap()
            .spectra,
            2
        );
    }
    for declaration in [
        "<?xml version=\"1.0\" version=\"1.0\"?>",
        "<?xml version=\"1.0\" standalone=\"invalid\"?>",
        "<?xml version=\"1.0\" standalone=\"yes\" encoding=\"UTF-8\"?>",
        "<?XML version=\"1.0\"?>",
        "<?1invalid?>",
    ] {
        assert!(
            mzml::read_size(Cursor::new(format!("{declaration}{}", document("", "")))).is_err()
        );
    }
    for prefix in [
        "&#32;",
        "<![CDATA[ ]]> ",
        "<?xml version=\"1.0\"?> <?xml version=\"1.0\"?>",
    ] {
        assert!(mzml::read_size(Cursor::new(format!("{prefix}{}", document("", "")))).is_err());
    }
}

#[test]
fn ignored_payload_xml_still_checks_structure_but_not_parameter_content() {
    let xml = document(
        "",
        "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><spectrum><cvParam>ignored semantic text</cvParam><referenceableParamGroupRef ref=\"absent\"/></spectrum></spectrumList>",
    );
    assert_eq!(mzml::read_size(Cursor::new(&xml)).unwrap().spectra, 1);
    assert!(mzml::read_size(Cursor::new(xml.replace("</cvParam>", "</wrong>"))).is_err());
    let valid = format!(
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>{}",
        document("", &lists(""))
    );
    assert_eq!(mzml::read_size(Cursor::new(&valid)).unwrap().spectra, 1);
    let invalid = valid.replace("<run", "<!--é--><run");
    assert!(matches!(
        mzml::read_size(Cursor::new(invalid)),
        Err(openms::Error::Unsupported(_))
    ));
    let sign = document(
        "",
        "<spectrumList count=\"+-2\" defaultDataProcessingRef=\"dp\"/>",
    );
    assert_eq!(mzml::read_size(Cursor::new(sign)).unwrap().spectra, 0);
}

#[test]
fn repeated_parameter_use_and_ms_membership_share_operation_budgets() {
    let header = group(&format!("{}{}", cv("1000016", "5"), cv("1000016", "6")));
    let xml = document(
        &header,
        &lists(&spectrum(
            "<scan><referenceableParamGroupRef ref=\"g\"/></scan>",
        )),
    );
    let limits = ReadOptions {
        max_total_params: 13,
        ..Default::default()
    };
    assert_eq!(
        mzml::read_size_with_options(Cursor::new(&xml), &filtered(), &limits)
            .unwrap()
            .spectra,
        2
    );
    let limits = ReadOptions {
        max_total_params: 12,
        ..limits
    };
    assert!(mzml::read_size_with_options(Cursor::new(&xml), &filtered(), &limits).is_err());
    let mut opt = PeakFileOptions::default();
    opt.set_ms_levels(&vec![1; 1_000_000]).unwrap();
    let body = spectrum(&format!("{}{}", cv("1000511", "1"), scan("5"))).repeat(51);
    let xml = document(
        "",
        &format!(
            "<spectrumList count=\"51\" defaultDataProcessingRef=\"dp\">{body}</spectrumList>"
        ),
    );
    assert!(
        count(&xml, &opt)
            .unwrap_err()
            .to_string()
            .contains("work limit")
    );
}

struct ErrorTail;
impl Read for ErrorTail {
    fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("unconsumed tail was read"))
    }
}
#[test]
fn successful_early_stop_does_not_drain_io_or_compressed_checksums() {
    let prefix = "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run><spectrumList count=\"7\" defaultDataProcessingRef=\"dp\"/><chromatogramList count=\"9\" defaultDataProcessingRef=\"dp\">";
    let stream = BufReader::with_capacity(1, Cursor::new(prefix).chain(ErrorTail));
    assert_eq!(
        mzml::read_size(stream).unwrap(),
        MzMLCounts {
            spectra: 7,
            chromatograms: 9
        }
    );
    let dir = Temp(std::env::temp_dir().join(format!("openms-counts-tail-{}", std::process::id())));
    std::fs::create_dir_all(&dir.0).unwrap();
    let tail = "A".repeat(100_000);
    let two = format!("{prefix}<binary>{tail}</binary></chromatogramList></run></mzML>");
    let one = document(
        "",
        &format!(
            "<spectrumList count=\"7\" defaultDataProcessingRef=\"dp\"><binary>{tail}</binary></spectrumList>"
        ),
    );
    for (name, xml, early) in [("two", two, true), ("one", one, false)] {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(xml.as_bytes()).unwrap();
        let mut encoded = encoder.finish().unwrap();
        let crc = encoded.len() - 8;
        encoded[crc] ^= 1;
        // Independently prove the test member has a failing trailer.
        assert!(
            flate2::read::MultiGzDecoder::new(encoded.as_slice())
                .read_to_end(&mut Vec::new())
                .is_err()
        );
        let path = dir.0.join(name);
        std::fs::write(&path, encoded).unwrap();
        if early {
            assert_eq!(mzml::load_size(&path).unwrap().spectra, 7);
        } else {
            assert!(mzml::load_size(&path).is_err());
        }
    }
}

#[test]
fn attribute_lexical_validation_also_applies_in_skipped_bodies() {
    for body in [
        "<x a=\"v\"b=\"w\"/>",
        "<x invalid&amp;key=\"v\"/>",
        "<x 1bad=\"v\"/>",
        "<x a=\"<\"/>",
        "<x a:b:c=\"v\"/>",
    ] {
        let xml = document(
            "",
            &format!(
                "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\">{body}</spectrumList>"
            ),
        );
        assert!(mzml::read_size(Cursor::new(xml)).is_err(), "{body}");
    }
    let xml = document(
        "",
        "<spectrumList count=\"1\" defaultDataProcessingRef=\"dp\"><x a=\"&lt;\" b='&quot;' xmlns:n=\"urn:attribute\" n:name=\"v\"/></spectrumList>",
    );
    assert_eq!(mzml::read_size(Cursor::new(xml)).unwrap().spectra, 1);
    let bad = format!(
        "<?xml version=\"1.0\"encoding=\"UTF-8\"?>{}",
        document("", "")
    );
    assert!(mzml::read_size(Cursor::new(bad)).is_err());
}

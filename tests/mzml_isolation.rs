// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use openms::format::mzml::{self, LoadOptions, ReadOptions, TransformOptions};
use openms::interfaces::MSDataConsumer;
use openms::kernel::NumericRange;
use openms::metadata::{ExperimentalSettings, MetaValueData};
use openms::{MSChromatogram, MSExperiment, MSSpectrum, Result};
use std::{
    io::{Cursor, Write},
    ops::ControlFlow,
};
fn cv(id: &str, value: &str) -> String {
    format!("<cvParam accession=\"MS:{id}\" name=\"fixture\" value=\"{value}\"/>")
}
fn target(value: &str) -> String {
    format!(
        "<isolationWindow>{}</isolationWindow>",
        cv("1000827", value)
    )
}
fn selected(ions: &[String]) -> String {
    format!(
        "<selectedIonList count=\"{}\">{}</selectedIonList>",
        ions.len(),
        ions.iter().fold(String::new(), |mut text, s| {
            text.push_str("<selectedIon>");
            text.push_str(s);
            text.push_str("</selectedIon>");
            text
        })
    )
}
fn precursor(contents: &str) -> String {
    format!("<precursor>{contents}<activation/></precursor>")
}
fn parts(t: Option<&str>, s: Option<&str>) -> String {
    precursor(&format!(
        "{}{}",
        t.map(target).unwrap_or_default(),
        selected(&s.map(|s| vec![cv("1000744", s)]).unwrap_or_default())
    ))
}
fn scan() -> String {
    "<scanList count=\"1\"><scan><cvParam accession=\"MS:1000016\" name=\"scan start time\" value=\"1\" unitAccession=\"UO:0000010\"/></scan></scanList>".into()
}
fn spec(id: usize, precursors: &[String], rt_first: bool) -> String {
    let p = format!(
        "<precursorList count=\"{}\">{}</precursorList>",
        precursors.len(),
        precursors.concat()
    );
    let fields = if rt_first {
        format!("{}{p}", scan())
    } else {
        format!("{p}{}", scan())
    };
    format!(
        "<spectrum id=\"s{id}\" index=\"{id}\" defaultArrayLength=\"0\">{}{fields}</spectrum>",
        cv("1000511", "2")
    )
}
fn chrom(p: &str) -> String {
    format!("<chromatogram id=\"c\" index=\"0\" defaultArrayLength=\"0\">{p}</chromatogram>")
}
fn doc(spectra: &[String], chromatograms: &[String], groups: &str) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">{groups}<fileDescription><fileContent/></fileDescription><softwareList count=\"1\"><software id=\"sw\" version=\"1\">{}</software></softwareList><instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\"/></instrumentConfigurationList><dataProcessingList count=\"1\"><dataProcessing id=\"dp\"><processingMethod order=\"0\" softwareRef=\"sw\">{}</processingMethod></dataProcessing></dataProcessingList><run id=\"r\" defaultInstrumentConfigurationRef=\"ic\"><spectrumList count=\"{}\" defaultDataProcessingRef=\"dp\">{}</spectrumList><chromatogramList count=\"{}\" defaultDataProcessingRef=\"dp\">{}</chromatogramList></run></mzML>",
        cv("1000799", "test"),
        cv("1000544", ""),
        spectra.len(),
        spectra.concat(),
        chromatograms.len(),
        chromatograms.concat()
    )
}
fn options() -> LoadOptions {
    let mut o = LoadOptions::default();
    o.scientific.precursor_mz_selected_ion = false;
    o
}
fn filtered() -> LoadOptions {
    let mut o = options();
    o.scientific.set_precursor_mz_range(NumericRange {
        min: 200.,
        max: 300.,
    });
    o
}
fn read(xml: &str, o: &LoadOptions) -> Result<MSExperiment> {
    mzml::read_with_load_options(Cursor::new(xml), o, &ReadOptions::default())
}
fn extra(p: &openms::Precursor) -> Option<f64> {
    p.cv_terms.metadata.get("selected ion m/z").map(|v| {
        assert!(matches!(v.data(), MetaValueData::Float(_)));
        assert!(v.unit().is_none());
        v.as_f64().unwrap()
    })
}
#[test]
fn literal_source_precursor_values_and_absent_equal_cases() {
    let fixture = include_str!("data/mzml_isolation_target.mzML");
    let e = read(fixture, &options()).unwrap();
    assert_eq!(e.spectra[0].precursors[0].mz, 500.);
    assert_eq!(extra(&e.spectra[0].precursors[0]), Some(499.5));
    let alias = fixture.replace("accession=\"MS:1000744\"", "accession=\"MS:1000040\"");
    assert_eq!(read(&alias, &options()).unwrap(), e);
    // Source test constructs target500/selected499.5; false-mode expectation is
    // independently derived from MzMLHandler's target/selected event branches.
    for (t, s, mz, other) in [
        (Some("500"), Some("499.5"), 500., Some(499.5)),
        (Some("500"), Some("500"), 500., None),
        (None, Some("499.5"), 0., Some(499.5)),
        (Some("500"), None, 500., None),
        (None, None, 0., None),
        (None, Some("-0"), 0., None),
    ] {
        let xml = doc(&[spec(0, &[parts(t, s)], false)], &[], "");
        let result = read(&xml, &options()).unwrap();
        let p = &result.spectra[0].precursors[0];
        assert_eq!(p.mz, mz);
        assert_eq!(extra(p), other);
        assert_eq!(p.isolation_target_mz, None);
    }
    let xml = doc(
        &[spec(0, &[parts(Some("500"), Some("499.5"))], false)],
        &[],
        "",
    );
    let legacy = mzml::read(Cursor::new(&xml)).unwrap();
    assert_eq!(legacy.spectra[0].precursors[0].mz, 499.5);
    assert_eq!(
        legacy.spectra[0].precursors[0].isolation_target_mz,
        Some(500.)
    );
    assert_eq!(extra(&legacy.spectra[0].precursors[0]), None);
    assert_eq!(read(&xml, &LoadOptions::default()).unwrap(), legacy);
}
#[test]
fn target_events_alone_filter_with_half_open_endpoints_and_missing_cv_gating() {
    for (t, s, keep) in [
        (Some("199.9"), Some("250"), false),
        (Some("200"), Some("999"), true),
        (Some("299.9"), Some("1"), true),
        (Some("300"), Some("250"), false),
        (None, Some("999"), true),
        (Some("300"), Some("300"), false),
        (Some("300"), None, false),
    ] {
        let xml = doc(&[spec(0, &[parts(t, s)], false)], &[], "");
        assert_eq!(
            read(&xml, &filtered()).unwrap().spectra.len(),
            usize::from(keep)
        );
        assert_eq!(
            mzml::read_size_with_options(
                Cursor::new(&xml),
                &filtered().scientific,
                &ReadOptions::default()
            )
            .unwrap()
            .spectra,
            usize::from(keep)
        );
    }
}
#[test]
fn repeated_target_events_are_sticky_even_when_the_final_value_passes() {
    for values in [["100", "250"], ["250", "100"], ["300", "200"]] {
        let p = precursor(&format!(
            "<isolationWindow>{}{}</isolationWindow>{}",
            cv("1000827", values[0]),
            cv("1000827", values[1]),
            selected(&[cv("1000744", "250")])
        ));
        let xml = doc(&[spec(0, &[p], false)], &[], "");
        assert!(read(&xml, &filtered()).unwrap().spectra.is_empty());
        assert_eq!(
            mzml::read_size_with_options(
                Cursor::new(&xml),
                &filtered().scientific,
                &ReadOptions::default()
            )
            .unwrap()
            .spectra,
            0
        );
        assert_eq!(
            read(&xml, &options()).unwrap().spectra[0].precursors[0].mz,
            values[1].parse::<f64>().unwrap()
        );
        assert!(read(&xml, &LoadOptions::default()).is_err());
    }
}
#[test]
fn repeated_selected_events_update_only_on_difference_and_target_order_is_literal() {
    for (events, expected) in [
        (
            format!("{}{}", cv("1000744", "260"), cv("1000744", "250")),
            Some(260.),
        ),
        (
            format!("{}{}", cv("1000744", "260"), cv("1000744", "270")),
            Some(270.),
        ),
    ] {
        let p = precursor(&format!("{}{}", target("250"), selected(&[events])));
        let e = read(&doc(&[spec(0, &[p], false)], &[], ""), &options()).unwrap();
        assert_eq!(e.spectra[0].precursors[0].mz, 250.);
        assert_eq!(extra(&e.spectra[0].precursors[0]), expected);
    }
    let p = precursor(&format!(
        "{}{}",
        selected(&[cv("1000744", "250")]),
        target("250")
    ));
    let e = read(&doc(&[spec(0, &[p], false)], &[], ""), &filtered()).unwrap();
    assert_eq!(extra(&e.spectra[0].precursors[0]), Some(250.));
}
#[test]
fn later_selected_ions_are_ignored_but_markup_groups_and_counts_are_checked() {
    let ignored = format!(
        "{}{}<userParam name=\"ignored\" type=\"xsd:unknown\" value=\"abc\"/>",
        cv("1000744", "NaN"),
        cv("1000041", "not-an-int")
    );
    let p = precursor(&format!(
        "{}{}",
        target("250"),
        selected(&[cv("1000744", "260"), ignored])
    ));
    let xml = doc(&[spec(0, &[p], false)], &[], "");
    let e = read(&xml, &filtered()).unwrap();
    let p = &e.spectra[0].precursors[0];
    assert_eq!(p.mz, 250.);
    assert_eq!(extra(p), Some(260.));
    assert_eq!(p.charge, 0);
    assert!(!p.cv_terms.metadata.contains_key("ignored"));
    assert!(read(&xml, &LoadOptions::default()).is_err());
    assert!(
        read(
            &xml.replace("type=\"xsd:unknown\"", "type=\"bad<xml\""),
            &filtered()
        )
        .is_err()
    );
    assert!(
        read(
            &xml.replace(
                "<userParam name=\"ignored\" type=\"xsd:unknown\" value=\"abc\"/>",
                "<referenceableParamGroupRef ref=\"missing\"/>"
            ),
            &filtered()
        )
        .is_err()
    );
    assert!(
        read(
            &xml.replace("selectedIonList count=\"2\"", "selectedIonList count=\"3\""),
            &filtered()
        )
        .is_err()
    );
}
#[test]
fn reference_groups_share_event_order_and_every_precursor_can_exclude_the_spectrum() {
    let groups = format!(
        "<referenceableParamGroupList count=\"3\"><referenceableParamGroup id=\"bad\">{}</referenceableParamGroup><referenceableParamGroup id=\"good\">{}</referenceableParamGroup><referenceableParamGroup id=\"sel\">{}</referenceableParamGroup></referenceableParamGroupList>",
        cv("1000827", "100"),
        cv("1000827", "250"),
        cv("1000744", "999")
    );
    let good = precursor(
        "<isolationWindow><referenceableParamGroupRef ref=\"good\"/></isolationWindow><selectedIonList count=\"1\"><selectedIon><referenceableParamGroupRef ref=\"sel\"/></selectedIon></selectedIonList>",
    );
    let bad = good.replace("ref=\"good\"", "ref=\"bad\"");
    for ps in [vec![good.clone(), bad.clone()], vec![bad, good.clone()]] {
        let xml = doc(&[spec(0, &ps, false)], &[], &groups);
        assert!(read(&xml, &filtered()).unwrap().spectra.is_empty());
        assert_eq!(
            mzml::read_size_with_options(
                Cursor::new(&xml),
                &filtered().scientific,
                &ReadOptions::default()
            )
            .unwrap()
            .spectra,
            0
        );
    }
    let e = read(
        &doc(&[spec(0, &[good.clone(), good], false)], &[], &groups),
        &filtered(),
    )
    .unwrap();
    assert_eq!(e.spectra[0].precursors.len(), 2);
    assert_eq!(extra(&e.spectra[0].precursors[1]), Some(999.));
}
#[test]
fn chromatogram_precursors_choose_target_without_whole_record_filtering() {
    let p = parts(Some("100"), Some("999"));
    let e = read(&doc(&[], &[chrom(&p)], ""), &filtered()).unwrap();
    assert_eq!(e.chromatograms.len(), 1);
    assert_eq!(e.chromatograms[0].precursor.mz, 100.);
    assert_eq!(extra(&e.chromatograms[0].precursor), Some(999.));
}
#[test]
fn count_reader_preserves_source_early_rt_count_while_full_reader_filters_later_target() {
    let xml = doc(
        &[spec(0, &[parts(Some("100"), Some("250"))], true)],
        &[],
        "",
    );
    assert!(read(&xml, &filtered()).unwrap().spectra.is_empty());
    assert_eq!(
        mzml::read_size_with_options(
            Cursor::new(&xml),
            &filtered().scientific,
            &ReadOptions::default()
        )
        .unwrap()
        .spectra,
        1
    );
}
#[derive(Default)]
struct Sink {
    spectra: Vec<MSSpectrum>,
    chroms: Vec<MSChromatogram>,
    size: Option<(usize, usize)>,
}
impl MSDataConsumer for Sink {
    fn set_expected_size(&mut self, s: usize, c: usize) -> Result<()> {
        self.size = Some((s, c));
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }
    fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.spectra.push(s.clone());
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, c: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        self.chroms.push(c.clone());
        Ok(ControlFlow::Continue(()))
    }
}
#[test]
fn consumer_pools_and_population_modes_reuse_the_same_target_selection() {
    let xml = doc(
        &[
            spec(0, &[parts(Some("100"), Some("250"))], false),
            spec(1, &[parts(Some("250"), Some("999"))], false),
        ],
        &[chrom(&parts(Some("100"), Some("999")))],
        "",
    );
    for fill in [false, true] {
        for pool in [1, 2] {
            let mut o = TransformOptions {
                load: filtered(),
                ..Default::default()
            };
            o.load.scientific.fill_data = fill;
            o.load.scientific.max_data_pool_size = pool;
            let mut s = Sink::default();
            let mut target = MSExperiment::default();
            let report =
                mzml::transform_from_into(|| Ok(Cursor::new(&xml)), &mut s, &mut target, &o)
                    .unwrap();
            assert_eq!(s.size, Some((2, 1)));
            assert_eq!(
                (report.delivered.spectra, report.delivered.chromatograms),
                (1, 1)
            );
            assert_eq!(s.spectra, target.spectra);
            assert_eq!(s.chroms, target.chromatograms);
            assert_eq!(s.spectra[0].precursors[0].mz, 250.);
            assert_eq!(extra(&s.spectra[0].precursors[0]), Some(999.));
        }
    }
}
#[test]
fn selection_and_parameter_limits_remain_cumulative_and_retained_failure_is_atomic() {
    let p = parts(Some("250"), Some("999"));
    let xml = doc(
        &[
            spec(0, std::slice::from_ref(&p), false),
            spec(1, std::slice::from_ref(&p), false),
        ],
        &[],
        "",
    );
    let mut o = filtered();
    o.max_selection_work = 1;
    let error = read(&xml, &o).unwrap_err();
    assert!(error.to_string().contains("selection resource"));
    let mut o = TransformOptions {
        load: filtered(),
        skip_first_pass: true,
        ..Default::default()
    };
    o.load.scientific.max_data_pool_size = 1;
    let mut s = Sink::default();
    let mut dst = MSExperiment {
        spectra: vec![MSSpectrum {
            name: "old".into(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let before = dst.clone();
    let invalid = xml.replace("id=\"s1\"", "id=\"s0\"");
    assert!(mzml::transform_from_into(|| Ok(Cursor::new(&invalid)), &mut s, &mut dst, &o).is_err());
    assert_eq!(s.spectra.len(), 1);
    assert_eq!(dst, before);
}
#[test]
fn plain_and_compressed_path_load_count_and_transform_use_target_mode() {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let dir = std::env::temp_dir().join(format!(
        "openms-isolation-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    struct Guard(std::path::PathBuf);
    impl Drop for Guard {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let _guard = Guard(dir.clone());
    let xml = doc(
        &[spec(0, &[parts(Some("250"), Some("999"))], false)],
        &[],
        "",
    );
    let mut gz = flate2::write::GzEncoder::new(vec![], flate2::Compression::default());
    gz.write_all(xml.as_bytes()).unwrap();
    let mut bz = bzip2::write::BzEncoder::new(vec![], bzip2::Compression::default());
    bz.write_all(xml.as_bytes()).unwrap();
    for (i, bytes) in [
        xml.as_bytes().to_vec(),
        gz.finish().unwrap(),
        bz.finish().unwrap(),
    ]
    .into_iter()
    .enumerate()
    {
        let path = dir.join(format!("case{i}.mzML"));
        std::fs::write(&path, bytes).unwrap();
        let e = mzml::load_with_options(&path, &filtered(), &ReadOptions::default()).unwrap();
        assert_eq!(e.spectra[0].precursors[0].mz, 250.);
        assert_eq!(
            mzml::load_size_with_options(&path, &filtered().scientific, &ReadOptions::default())
                .unwrap()
                .spectra,
            1
        );
        let mut sink = Sink::default();
        mzml::transform_with_options(
            &path,
            &mut sink,
            &TransformOptions {
                load: filtered(),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(extra(&sink.spectra[0].precursors[0]), Some(999.));
    }
}

#[test]
fn every_writer_preserves_isolation_and_selected_values_without_activation_shadow() {
    let synthetic = doc(
        &[spec(0, &[parts(Some("250"), Some("999"))], false)],
        &[chrom(&parts(Some("250"), Some("999")))],
        "",
    );
    for xml in [
        synthetic.as_str(),
        include_str!("data/mzml_isolation_target.mzML"),
    ] {
        let mut e = read(xml, &options()).unwrap();
        // Synthetic read fixtures use arbitrary IDs; writer's established native
        // ID contract requires the key=value form. Only the fixture ID changes.
        for (index, spectrum) in e.spectra.iter_mut().enumerate() {
            spectrum.native_id = format!("scan={index}");
        }
        let target = e.spectra[0].precursors[0].mz;
        let ion = extra(&e.spectra[0].precursors[0]).unwrap();
        for mode in 0..4 {
            let mut bytes = Vec::new();
            match mode {
                0 => mzml::write(&mut bytes, &e).unwrap(),
                1 => {
                    mzml::write_with_numpress(&mut bytes, &e, &Default::default()).unwrap();
                }
                _ => {
                    let mut options = openms::format::peak_options::PeakFileOptions::default();
                    options.force_tpp_compatibility = mode == 3;
                    mzml::write_with_peak_options(&mut bytes, &e, &options).unwrap();
                }
            }
            let encoded = String::from_utf8(bytes).unwrap();
            assert!(!encoded.contains("<userParam name=\"selected ion m/z\""));
            let default = read(&encoded, &LoadOptions::default()).unwrap();
            let isolation = read(&encoded, &options()).unwrap();
            for p in std::iter::once(&default.spectra[0].precursors[0])
                .chain(default.chromatograms.iter().map(|c| &c.precursor))
            {
                assert_eq!(p.mz, ion);
                assert_eq!(
                    p.isolation_target_mz,
                    if mode == 3 { None } else { Some(target) }
                );
                assert_eq!(extra(p), None);
            }
            for p in std::iter::once(&isolation.spectra[0].precursors[0])
                .chain(isolation.chromatograms.iter().map(|c| &c.precursor))
            {
                assert_eq!(p.mz, if mode == 3 { 0. } else { target });
                assert_eq!(extra(p), Some(ion));
            }
        }
    }
}

#[test]
fn repeated_group_expansion_respects_cumulative_parameter_limits() {
    let groups = format!(
        "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"s\">{}</referenceableParamGroup></referenceableParamGroupList>",
        cv("1000744", "999")
    );
    let one = "<referenceableParamGroupRef ref=\"s\"/>";
    let build = |n: usize| {
        doc(
            &[spec(
                0,
                &[precursor(&format!(
                    "{}{}",
                    target("250"),
                    selected(&[one.repeat(n)])
                ))],
                false,
            )],
            &[],
            &groups,
        )
    };
    let small = build(1);
    let limits = ReadOptions {
        max_total_params: 40,
        ..Default::default()
    };
    assert!(mzml::read_with_load_options(Cursor::new(small), &options(), &limits).is_ok());
    let error =
        mzml::read_with_load_options(Cursor::new(build(40)), &options(), &limits).unwrap_err();
    assert!(error.to_string().contains("parameter count"));
}

#[test]
fn raw_attribute_guard_applies_to_group_definitions_but_escaped_less_than_is_valid() {
    let make = |value: &str| {
        let groups = format!(
            "<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"s\"><userParam name=\"ignored\" value=\"{value}\"/></referenceableParamGroup></referenceableParamGroupList>"
        );
        doc(
            &[spec(
                0,
                &[precursor(&format!(
                    "{}{}",
                    target("250"),
                    selected(&[
                        cv("1000744", "260"),
                        "<referenceableParamGroupRef ref=\"s\"/>".into()
                    ])
                ))],
                false,
            )],
            &[],
            &groups,
        )
    };
    assert!(read(&make("a<b"), &options()).is_err());
    assert!(read(&make("a&lt;b"), &options()).is_ok());
    let inline = doc(
        &[spec(
            0,
            &[precursor(&format!(
                "{}{}",
                target("250"),
                selected(&[
                    cv("1000744", "260"),
                    "<userParam name=\"ignored\" value=\"a&lt;b\"/>".into()
                ])
            ))],
            false,
        )],
        &[],
        "",
    );
    assert!(read(&inline, &options()).is_ok());
}
#[test]
fn selected_metadata_storage_is_precharged_and_filtered_records_remain_validated() {
    let make = |selected: &str| {
        doc(
            &[spec(0, &[parts(Some("250"), Some(selected))], false)],
            &[],
            "",
        )
    };
    let minimum = |xml: &str| {
        let (mut low, mut high) = (0, 1_000_000);
        while low < high {
            let mid = low + (high - low) / 2;
            let limits = ReadOptions {
                max_param_bytes: mid,
                ..Default::default()
            };
            if mzml::read_with_load_options(Cursor::new(xml), &options(), &limits).is_ok() {
                high = mid
            } else {
                low = mid + 1
            }
        }
        low
    };
    let equal = minimum(&make("250"));
    let different = minimum(&make("999"));
    assert!(different > equal + 1000);
    assert!(
        mzml::read_with_load_options(
            Cursor::new(make("999")),
            &options(),
            &ReadOptions {
                max_param_bytes: different - 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    let bad = doc(
        &[spec(0, &[parts(Some("100"), Some("NaN"))], false)],
        &[],
        "",
    );
    assert!(read(&bad, &filtered()).is_err());
}

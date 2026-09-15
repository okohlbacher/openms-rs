// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use base64::{Engine, engine::general_purpose::STANDARD};
use openms::format::mzml::{self, LoadOptions, ReadOptions};
use openms::kernel::NumericRange;
use openms::{Error, MSExperiment};
use std::io::Cursor;

const SOURCE: &str = include_str!("data/mzml_load_source_projection.mzML");
fn read(xml: &str, options: &LoadOptions) -> openms::Result<MSExperiment> {
    mzml::read_with_load_options(Cursor::new(xml), options, &ReadOptions::default())
}
fn range(min: f64, max: f64) -> NumericRange {
    NumericRange { min, max }
}
fn cv(id: &str, value: &str) -> String {
    let unit = if matches!(id, "MS:1000016" | "MS:1000595") {
        " unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\""
    } else {
        ""
    };
    format!("<cvParam accession=\"{id}\" name=\"fixture\" value=\"{value}\"{unit}/>")
}

fn array(kind: &str, encoding: &str, data: &[u8]) -> String {
    let text = STANDARD.encode(data);
    format!(
        "<binaryDataArray encodedLength=\"{}\">{}{}{}<binary>{text}</binary></binaryDataArray>",
        text.len(),
        cv(kind, ""),
        cv(encoding, ""),
        cv("MS:1000576", "")
    )
}
fn float_array(kind: &str, values: &[f64]) -> String {
    array(
        kind,
        "MS:1000523",
        &values
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<_>>(),
    )
}
fn record(
    kind: &str,
    id: usize,
    position: &[f64],
    intensity: &[f64],
    fields: &str,
    extras: &[String],
) -> String {
    let coordinate = if kind == "spectrum" {
        "MS:1000514"
    } else {
        "MS:1000595"
    };
    format!(
        "<{kind} id=\"scan={id}\" index=\"{id}\" defaultArrayLength=\"{}\">{fields}<binaryDataArrayList count=\"{}\">{}{}{}</binaryDataArrayList></{kind}>",
        position.len(),
        2 + extras.len(),
        float_array(coordinate, position),
        float_array("MS:1000515", intensity),
        extras.concat()
    )
}
fn document(spectra: &[String], chromatograms: &[String]) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\"><run id=\"r\"><userParam name=\"project\" value=\"preserved\"/><spectrumList count=\"{}\">{}</spectrumList><chromatogramList count=\"{}\">{}</chromatogramList></run></mzML>",
        spectra.len(),
        spectra.concat(),
        chromatograms.len(),
        chromatograms.concat()
    )
}
fn rt(value: f64) -> String {
    format!(
        "<scanList count=\"1\"><scan>{}</scan></scanList>",
        cv("MS:1000016", &value.to_string())
    )
}
fn precursor(isolation: Option<f64>, selected: f64) -> String {
    format!(
        "<precursor>{}<selectedIonList count=\"1\"><selectedIon>{}</selectedIon></selectedIonList><activation/></precursor>",
        isolation
            .map(|v| format!(
                "<isolationWindow>{}</isolationWindow>",
                cv("MS:1000827", &v.to_string())
            ))
            .unwrap_or_default(),
        cv("MS:1000744", &selected.to_string())
    )
}

#[test]
fn source_ms_level_rt_mz_and_intensity_literals() {
    let mut o = LoadOptions::default();
    o.scientific.add_ms_level(1).unwrap();
    let exp = read(SOURCE, &o).unwrap();
    assert_eq!(
        exp.spectra.iter().map(|s| s.rt).collect::<Vec<_>>(),
        [5.1, 5.3, 5.4]
    );
    o.scientific.clear_ms_levels();
    o.scientific.set_rt_range(range(5.15, 5.35));
    assert_eq!(
        read(SOURCE, &o)
            .unwrap()
            .spectra
            .iter()
            .map(|s| s.rt)
            .collect::<Vec<_>>(),
        [5.2, 5.3]
    );
    o = LoadOptions::default();
    o.scientific.set_mz_range(range(6.5, 9.5));
    assert_eq!(
        read(SOURCE, &o)
            .unwrap()
            .spectra
            .iter()
            .map(|s| s.peaks.iter().map(|p| p.mz).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        [vec![7., 8., 9.], vec![8.], vec![7., 8., 9.], vec![]]
    );
    o = LoadOptions::default();
    o.scientific.set_intensity_range(range(6.5, 9.5));
    assert_eq!(
        read(SOURCE, &o)
            .unwrap()
            .spectra
            .iter()
            .map(|s| s.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        [vec![9., 8., 7.], vec![8.], vec![9., 8., 7.], vec![]]
    );
}

#[test]
fn all_source_primary_values_match_independent_binary64_projection() {
    let exp = read(SOURCE, &LoadOptions::default()).unwrap();
    assert_eq!((exp.spectra.len(), exp.chromatograms.len()), (4, 2));
    let mut count = 0;
    for row in include_str!("data/mzml_load_source_values.tsv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split('\t').collect();
        let r: usize = f[1].parse().unwrap();
        let i: usize = f[2].parse().unwrap();
        let p = f64::from_bits(u64::from_str_radix(f[3], 16).unwrap());
        let intensity = f64::from_bits(u64::from_str_radix(f[4], 16).unwrap()) as f32;
        if f[0] == "spectrum" {
            assert_eq!(
                (
                    exp.spectra[r].peaks[i].mz,
                    exp.spectra[r].peaks[i].intensity
                ),
                (p, intensity)
            );
        } else {
            assert_eq!(
                (
                    exp.chromatograms[r].peaks[i].rt,
                    exp.chromatograms[r].peaks[i].intensity
                ),
                (p, intensity)
            );
        }
        count += 1;
    }
    assert_eq!(count, 65);
    assert_eq!(
        exp.spectra[1].float_data_arrays[0].name,
        "signal to noise array"
    );
    assert_eq!(
        exp.spectra[1].float_data_arrays[1].name,
        "user-defined name"
    );
}

#[test]
fn raw_double_filters_precede_narrowing_and_use_half_open_bounds() {
    let xml = document(
        &[record(
            "spectrum",
            0,
            &[1., 2., 3., 4.],
            &[1., 1. + 1e-9, 2., f64::MAX],
            "",
            &[],
        )],
        &[],
    );
    let mut o = LoadOptions::default();
    o.scientific.set_intensity_range(range(1. + 5e-10, 2.));
    let exp = read(&xml, &o).unwrap();
    assert_eq!(exp.spectra[0].peaks.len(), 1);
    assert_eq!(
        (
            exp.spectra[0].peaks[0].mz,
            exp.spectra[0].peaks[0].intensity
        ),
        (2., 1.)
    );
    // f64::MAX is legal decoded input; filtering it out avoids f32 overflow.
    assert!(mzml::read(Cursor::new(&xml)).is_err());
    o.scientific.set_mz_range(range(2., 2.));
    assert!(read(&xml, &o).unwrap().spectra[0].is_empty());
}

#[test]
fn sorts_and_filters_all_array_kinds_without_losing_metadata_or_ties() {
    let float = float_array("MS:1000517", &[30., 10., 20., 11.]);
    let integer = array(
        "MS:1000516",
        "MS:1000519",
        &[3_i32, 1, 2, 4]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect::<Vec<_>>(),
    );
    let string = array("MS:1000786", "MS:1001479", b"c\0a\0b\0d\0").replacen(
        "value=\"\"",
        "value=\"labels\"",
        1,
    );
    let fields = format!(
        "{}{}<userParam name=\"note\" value=\"unchanged\"/>",
        cv("MS:1000511", "2"),
        rt(8.)
    );
    let xml = document(
        &[record(
            "spectrum",
            0,
            &[3., 1., 2., 1.],
            &[30., 10., 20., 11.],
            &fields,
            &[float, integer, string],
        )],
        &[],
    );
    let mut o = LoadOptions::default();
    o.scientific.set_mz_range(range(1., 3.));
    let exp = read(&xml, &o).unwrap();
    let s = &exp.spectra[0];
    assert_eq!(
        s.peaks.iter().map(|p| p.intensity).collect::<Vec<_>>(),
        [10., 11., 20.]
    );
    assert_eq!(s.float_data_arrays[0].data, [10., 11., 20.]);
    assert_eq!(s.integer_data_arrays[0].data, [1, 4, 2]);
    assert_eq!(s.string_data_arrays[0].data, ["a", "d", "b"]);
    assert_eq!((s.rt, s.ms_level), (8., 2));
    assert_eq!(s.metadata["note"].as_str().unwrap(), "unchanged");
    assert_eq!(
        exp.settings.metadata["project"].as_str().unwrap(),
        "preserved"
    );
    o.scientific.sort_spectra_by_mz = false;
    assert_eq!(
        read(&xml, &o).unwrap().spectra[0]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [10., 20., 11.]
    );
    assert_eq!(
        mzml::read(Cursor::new(xml)).unwrap().spectra[0].peaks[0].mz,
        3.
    );
}

#[test]
fn chromatogram_rt_and_intensity_filters_do_not_apply_spectrum_mz_ranges() {
    let xml = document(
        &[],
        &[record(
            "chromatogram",
            0,
            &[4., 1., 3., 2.],
            &[40., 10., 30., 20.],
            "",
            &[float_array("MS:1000517", &[4., 1., 3., 2.])],
        )],
    );
    let mut o = LoadOptions::default();
    o.scientific.set_rt_range(range(2., 4.));
    o.scientific.set_mz_range(range(100., 200.));
    o.scientific.set_intensity_range(range(0., 30.));
    let exp = read(&xml, &o).unwrap();
    assert_eq!(exp.chromatograms[0].peaks.len(), 1);
    assert_eq!(exp.chromatograms[0].peaks[0].rt, 2.);
    assert_eq!(exp.chromatograms[0].float_data_arrays[0].data, [2.]);
    o = LoadOptions::default();
    o.scientific.sort_chromatograms_by_rt = false;
    assert_eq!(read(&xml, &o).unwrap().chromatograms[0].peaks[0].rt, 4.);
}

#[test]
fn missing_ms_level_and_rt_do_not_activate_source_cv_filters() {
    let xml = document(
        &[
            record("spectrum", 0, &[], &[], "", &[]),
            record(
                "spectrum",
                1,
                &[],
                &[],
                &format!("{}{}", cv("MS:1000511", "2"), rt(4.)),
                &[],
            ),
        ],
        &[],
    );
    let mut o = LoadOptions::default();
    o.scientific.add_ms_level(9).unwrap();
    o.scientific.set_rt_range(range(100., 200.));
    assert_eq!(read(&xml, &o).unwrap().spectra.len(), 1);
    let duplicate = xml.replace(&cv("MS:1000511", "2"), &cv("MS:1000511", "2").repeat(2));
    assert!(read(&duplicate, &o).is_err());
}

#[test]
fn precursor_filter_retains_actual_source_equal_target_and_zero_bypass() {
    let cases = [
        (Some(100.), 100.),
        (Some(99.), 100.),
        (None, 0.),
        (None, 150.),
    ];
    let specs: Vec<_> = cases
        .into_iter()
        .enumerate()
        .map(|(i, (iso, sel))| {
            record(
                "spectrum",
                i,
                &[],
                &[],
                &format!(
                    "{}<precursorList count=\"1\">{}</precursorList>",
                    cv("MS:1000511", "2"),
                    precursor(iso, sel)
                ),
                &[],
            )
        })
        .collect();
    let mut o = LoadOptions::default();
    o.scientific.set_precursor_mz_range(range(140., 160.));
    let exp = read(&document(&specs, &[]), &o).unwrap();
    assert_eq!(
        exp.spectra
            .iter()
            .map(|s| s.native_id.as_str())
            .collect::<Vec<_>>(),
        ["scan=0", "scan=2", "scan=3"]
    );
    let xml = document(
        &[record(
            "spectrum",
            0,
            &[],
            &[],
            &format!(
                "<precursorList count=\"2\">{}{}</precursorList>",
                precursor(None, 150.),
                precursor(None, 100.)
            ),
            &[],
        )],
        &[],
    );
    assert!(read(&xml, &o).unwrap().spectra.is_empty());
}

#[test]
fn whole_record_skips_keep_raw_count_checks_and_validate_excluded_payload() {
    let mut o = LoadOptions {
        skip_spectra: true,
        ..Default::default()
    };
    o.scientific.skip_chromatograms = true;
    let exp = read(SOURCE, &o).unwrap();
    assert!(exp.spectra.is_empty() && exp.chromatograms.is_empty());
    assert!(read(&SOURCE.replacen("<binary>", "<binary>!", 1), &o).is_err());
    let limits = ReadOptions {
        max_total_peaks: 64,
        ..Default::default()
    };
    assert!(mzml::read_with_load_options(Cursor::new(SOURCE), &o, &limits).is_err());
    // The declared record count is advisory even when every record is skipped:
    // MzMLHandler.cpp:965-979 spends it on a progress range and
    // `reserveSpaceSpectra` and never compares it. The raw peak, binary and
    // parameter ceilings above are what still bound a skipped record.
    assert_eq!(
        read(
            &SOURCE.replacen("spectrumList count=\"4\"", "spectrumList count=\"3\"", 1),
            &o
        )
        .unwrap(),
        exp
    );
    assert!(
        read(
            &SOURCE.replacen("spectrumList count=\"4\"", "spectrumList count=\"\"", 1),
            &o
        )
        .is_err()
    );
}

#[test]
fn selection_budget_is_cumulative_across_records_and_independent_of_binary_caps() {
    let xml = document(
        &[
            record("spectrum", 0, &[2., 1.], &[2., 1.], "", &[]),
            record("spectrum", 1, &[2., 1.], &[2., 1.], "", &[]),
        ],
        &[],
    );
    let per = 2 * 3 * std::mem::size_of::<usize>();
    let mut o = LoadOptions {
        max_selection_bytes: per,
        ..Default::default()
    };
    assert!(
        matches!(read(&xml,&o),Err(Error::InvalidValue(s)) if s.contains("selection resource"))
    );
    o.max_selection_bytes = 2 * per;
    assert_eq!(read(&xml, &o).unwrap().spectra.len(), 2);
    o.max_selection_work = 0;
    assert!(read(&xml, &o).is_err());
}

#[test]
fn write_flags_are_ignored_on_load() {
    let mut o = LoadOptions::default();
    o.scientific.zlib_compression = true;
    o.scientific.mz_32_bit = true;
    o.scientific.intensity_32_bit = false;
    o.scientific.force_tpp_compatibility = true;
    o.scientific.force_mq_compatibility = true;
    o.scientific.write_supplemental_data = false;
    o.scientific.max_data_pool_size = 0;
    o.scientific.always_append_data = true;
    assert_eq!(
        read(SOURCE, &o).unwrap(),
        read(SOURCE, &LoadOptions::default()).unwrap()
    );
}

#[test]
fn all_pinned_canonical_names_and_declared_kinds_round_trip() {
    let mut arrays = Vec::new();
    let mut expected = Vec::new();
    for row in include_str!("data/mzml_load_canonical_arrays.tsv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        expected.push((fields[0], fields[1]));
        arrays.push(match fields[2] {
            "1" => array(fields[0], "MS:1000519", &3_i32.to_le_bytes()),
            "2" => array(fields[0], "MS:1000521", &3_f32.to_le_bytes()),
            _ => float_array(fields[0], &[3.]),
        });
    }
    assert_eq!(expected.len(), 26);
    let xml = document(&[record("spectrum", 0, &[3.], &[4.], "", &arrays)], &[]);
    let parsed = read(&xml, &LoadOptions::default()).unwrap();
    let mut output = Vec::new();
    mzml::write(&mut output, &parsed).unwrap();
    let output = String::from_utf8(output).unwrap();
    for (id, name) in expected {
        assert!(output.contains(&format!("accession=\"{id}\" name=\"{name}\"")));
    }
    assert_eq!(mzml::read(Cursor::new(output)).unwrap(), parsed);
}

#[test]
fn canonical_collisions_units_and_wrong_kinds_fail_before_writer_output() {
    let named = float_array("MS:1000786", &[1.]).replacen(
        "value=\"\"",
        "value=\"signal to noise array\"",
        1,
    );
    let wrong = array("MS:1000517", "MS:1001479", b"a\0");
    let units = float_array("MS:1000517", &[1.]).replacen(
        "name=\"fixture\"",
        "name=\"fixture\" unitName=\"counts\"",
        1,
    );
    for extra in [named, wrong, units] {
        assert!(matches!(
            read(
                &document(&[record("spectrum", 0, &[1.], &[1.], "", &[extra])], &[]),
                &LoadOptions::default()
            ),
            Err(Error::Unsupported(_))
        ));
    }
    let mut exp = read(
        &document(&[record("spectrum", 0, &[1.], &[1.], "", &[])], &[]),
        &LoadOptions::default(),
    )
    .unwrap();
    exp.spectra[0]
        .string_data_arrays
        .push(openms::kernel::DataArray::new(
            "signal to noise array",
            vec!["a".into()],
        ));
    let mut output = Vec::new();
    assert!(mzml::write(&mut output, &exp).is_err());
    assert!(output.is_empty());
}

#[test]
fn reference_groups_share_filter_events_and_empty_placeholders_remain_empty() {
    let extra = array("MS:1000786", "MS:1001479", b"")
        .replacen("value=\"\"", "value=\"placeholder\"", 1)
        .replacen(
            "<binaryDataArray ",
            "<binaryDataArray arrayLength=\"0\" ",
            1,
        );
    let xml=document(&[record("spectrum",0,&[3.,1.,2.],&[3.,1.,2.],"<referenceableParamGroupRef ref=\"ms2\"/>",&[extra])],&[]).replace("<run id=\"r\">",&format!("<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"ms2\">{}</referenceableParamGroup></referenceableParamGroupList><run id=\"r\">",cv("MS:1000511","2")));
    let mut o = LoadOptions::default();
    o.scientific.add_ms_level(2).unwrap();
    o.scientific.set_mz_range(range(1., 3.));
    let parsed = read(&xml, &o).unwrap();
    assert_eq!(parsed.spectra[0].peaks.len(), 2);
    assert!(parsed.spectra[0].string_data_arrays[0].data.is_empty());
    o.scientific.set_ms_levels(&[1]).unwrap();
    assert!(read(&xml, &o).unwrap().spectra.is_empty());
}

#[test]
fn finite_signed_and_source_nan_range_endpoints_are_not_normalized() {
    let xml = document(
        &[record(
            "spectrum",
            0,
            &[-2., 0., 2.],
            &[-2., 0., 2.],
            "",
            &[],
        )],
        &[],
    );
    let mut o = LoadOptions::default();
    o.scientific.set_mz_range(range(-2., 0.));
    assert_eq!(read(&xml, &o).unwrap().spectra[0].peaks[0].mz, -2.);
    o.scientific.set_mz_range(range(f64::NAN, 0.));
    assert_eq!(read(&xml, &o).unwrap().spectra[0].peaks.len(), 1);
    o.scientific.set_mz_range(range(2., -2.));
    assert!(read(&xml, &o).unwrap().spectra[0].peaks.is_empty());
}

#[test]
fn minute_units_are_scaled_before_both_spectrum_and_chromatogram_rt_filters() {
    let mut o = LoadOptions::default();
    o.scientific.set_rt_range(range(60., 91.));
    let exp = read(include_str!("data/mzml_independent.mzML"), &o).unwrap();
    assert_eq!(exp.spectra[0].rt, 90.);
    assert_eq!(
        exp.chromatograms[0]
            .peaks
            .iter()
            .map(|p| p.rt)
            .collect::<Vec<_>>(),
        [60., 90.]
    );
}

#[test]
fn precursor_selection_work_is_charged_before_inline_or_referenced_cv_execution() {
    let mut o = LoadOptions {
        max_selection_work: 1,
        ..Default::default()
    };
    o.scientific.set_precursor_mz_range(range(1., 2.));
    let fields = format!(
        "<precursorList count=\"1\">{}</precursorList>",
        precursor(None, 1.5)
    );
    let xml = document(&[record("spectrum", 0, &[], &[], &fields, &[])], &[]);
    // One unit initializes selection; the very first selected-ion filter then
    // exhausts before the deliberately malformed binary payload is reached.
    let xml = xml.replacen("<binary>", "<binary>!", 1);
    assert!(
        matches!(read(&xml,&o),Err(Error::InvalidValue(s)) if s.contains("selection resource"))
    );
    let term = cv("MS:1000744", "1.5");
    let xml=xml.replace(&term,"<referenceableParamGroupRef ref=\"p\"/>").replace("<run id=\"r\">",&format!("<referenceableParamGroupList count=\"1\"><referenceableParamGroup id=\"p\">{term}</referenceableParamGroup></referenceableParamGroupList><run id=\"r\">"));
    assert!(
        matches!(read(&xml,&o),Err(Error::InvalidValue(s)) if s.contains("selection resource"))
    );
}

#[test]
fn auxiliary_double_values_are_narrowed_only_after_aligned_peak_selection() {
    let auxiliary = float_array("MS:1000517", &[f64::MAX, 1. + 1e-9]);
    let xml = document(
        &[record(
            "spectrum",
            0,
            &[1., 2.],
            &[1., 2.],
            "",
            &[auxiliary],
        )],
        &[],
    );
    let mut o = LoadOptions::default();
    o.scientific.set_mz_range(range(2., 3.));
    let exp = read(&xml, &o).unwrap();
    assert_eq!(exp.spectra[0].float_data_arrays[0].data, [1.]);
    assert!(mzml::read(Cursor::new(&xml)).is_err());
    o.skip_spectra = true;
    // Native full-record validation still applies when the entire record is
    // skipped, distinct from source's retained-peak conversion above.
    assert!(read(&xml, &o).is_err());
}

#[test]
fn new_loader_rejects_unresolved_array_processing_references_even_if_excluded() {
    let xml = document(&[record("spectrum", 0, &[], &[], "", &[])], &[]).replacen(
        "<binaryDataArray encodedLength",
        "<binaryDataArray dataProcessingRef=\"processing\" encodedLength",
        1,
    );
    let o = LoadOptions {
        skip_spectra: true,
        ..Default::default()
    };
    assert!(
        matches!(read(&xml,&o),Err(Error::Parse{message:s,..}) if s.contains("unresolved dataProcessingRef"))
    );
}

/// A chromatogram whose time array is 32-bit with an explicit unit.
///
/// The intensity array stays 32-bit, as ProteoWizard writes it.
fn timed_chromatogram(minutes: bool, times: &[f32], intensities: &[f32]) -> String {
    let unit = if minutes {
        " unitAccession=\"UO:0000031\" unitName=\"minute\" unitCvRef=\"UO\""
    } else {
        " unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\""
    };
    let bytes = |v: &[f32]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<_>>();
    let array32 = |accession: &str, unit: &str, data: Vec<u8>| {
        let text = STANDARD.encode(&data);
        format!(
            "<binaryDataArray encodedLength=\"{}\"><cvParam accession=\"{accession}\" name=\"fixture\" value=\"\"{unit}/>{}{}<binary>{text}</binary></binaryDataArray>",
            text.len(),
            cv("MS:1000521", ""),
            cv("MS:1000576", "")
        )
    };
    format!(
        "<chromatogram id=\"TIC\" index=\"0\" defaultArrayLength=\"{}\"><binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList></chromatogram>",
        times.len(),
        array32("MS:1000595", unit, bytes(times)),
        array32("MS:1000515", "", bytes(intensities))
    )
}

/// `source_time_array_precision` reproduces the source's in-place minute
/// conversion of a 32-bit time array, and touches nothing else.
///
/// Source `MzMLHandlerHelper::decodeBase64Arrays` multiplies a 32-bit array in
/// place through `float&` (`MzMLHandlerHelper.cpp:217-222`), so the seconds it
/// stores are the `f32` rounding of the `double` product; the 64-bit branch
/// above it keeps the `double`. The executed differential over the real
/// `UK222.mzML` TIC arrays is `tests/peak_picking_experiment.rs`
/// (`chromatogram_time_*`); this pins the switch itself and the shapes it must
/// leave alone.
#[test]
fn source_time_array_precision_narrows_only_a_converted_32_bit_time_array() {
    let times: Vec<f32> = vec![1.0019801_f32, 1.0071758, 68.336365, 73.2];
    let intensities = vec![1.0_f32, 2.0, 3.0, 4.0];
    let read_times = |xml: &str, source: bool| -> Vec<f64> {
        let options = if source {
            ReadOptions::source()
        } else {
            ReadOptions::default()
        };
        mzml::read_with_options(Cursor::new(xml.to_owned()), &options)
            .unwrap()
            .chromatograms[0]
            .peaks
            .iter()
            .map(|peak| peak.rt)
            .collect()
    };

    let minutes = document(&[], &[timed_chromatogram(true, &times, &intensities)]);
    let native = read_times(&minutes, false);
    let source = read_times(&minutes, true);
    let exact: Vec<f64> = times.iter().map(|&t| f64::from(t) * 60.0).collect();
    assert_eq!(native, exact, "the default keeps the f64 product");
    let narrowed: Vec<f64> = exact.iter().map(|&t| f64::from(t as f32)).collect();
    assert_eq!(source, narrowed, "source mode keeps only f32");
    assert_ne!(native, source, "the fixture must exercise the difference");

    // Seconds: the source sets no multiplier, so neither mode narrows.
    let seconds = document(&[], &[timed_chromatogram(false, &times, &intensities)]);
    let expected: Vec<f64> = times.iter().map(|&t| f64::from(t)).collect();
    assert_eq!(read_times(&seconds, false), expected);
    assert_eq!(read_times(&seconds, true), expected);

    // An empty 32-bit minute array has nothing to narrow and still reads.
    let empty = document(&[], &[timed_chromatogram(true, &[], &[])]);
    assert!(read_times(&empty, true).is_empty());
    assert!(read_times(&empty, false).is_empty());
}

/// A 64-bit time array in minutes keeps its `double` product in both modes.
///
/// The switch is keyed on the source's own condition, which is the 32-bit
/// branch of `decodeBase64Arrays`; a 64-bit array never loses precision there.
#[test]
fn source_time_array_precision_leaves_a_64_bit_time_array_alone() {
    let minutes = [1.0019801_f64, 68.336365, 73.2];
    let chromatogram = record(
        "chromatogram",
        0,
        &minutes,
        &[1.0, 2.0, 3.0],
        "",
        &[],
    )
    .replace(
        "accession=\"MS:1000595\" name=\"fixture\" value=\"\" unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\"",
        "accession=\"MS:1000595\" name=\"fixture\" value=\"\" unitAccession=\"UO:0000031\" unitName=\"minute\" unitCvRef=\"UO\"",
    );
    let xml = document(&[], &[chromatogram]);
    let expected: Vec<f64> = minutes.iter().map(|&m| m * 60.0).collect();
    for source in [false, true] {
        let options = if source {
            ReadOptions::source()
        } else {
            ReadOptions::default()
        };
        let exp = mzml::read_with_options(Cursor::new(xml.clone()), &options).unwrap();
        let times: Vec<f64> = exp.chromatograms[0].peaks.iter().map(|p| p.rt).collect();
        assert_eq!(times, expected, "source={source}");
    }
}

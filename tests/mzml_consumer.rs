#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use base64::{Engine as _, engine::general_purpose::STANDARD};
use openms::{
    Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Result,
    format::mzml::{self, MzMLCounts, TransformOptions},
    interfaces::MSDataConsumer,
    metadata::ExperimentalSettings,
};
use std::{io::Cursor, ops::ControlFlow};
const SOURCE: &str = include_str!("data/mzml_consumer/source.mzML");

#[derive(Default)]
struct Log {
    events: Vec<String>,
    spectra: Vec<MSSpectrum>,
    chromatograms: Vec<MSChromatogram>,
    settings: Option<ExperimentalSettings>,
    stop: Option<usize>,
    fail: Option<usize>,
    mutate: bool,
}
impl Log {
    fn event(&mut self, value: String) -> Result<()> {
        self.events.push(value);
        if self.fail == Some(self.events.len()) {
            return Err(Error::InvalidValue("callback failed".into()));
        }
        Ok(())
    }
    fn flow(&self) -> ControlFlow<()> {
        if self.stop == Some(self.spectra.len() + self.chromatograms.len()) {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    }
}
impl MSDataConsumer for Log {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        self.event(format!("size:{spectra}:{chromatograms}"))
    }
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        self.event("settings".into())?;
        self.settings = Some(settings.clone());
        Ok(())
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.event(spectrum.native_id.clone())?;
        if self.mutate {
            spectrum.name = "changed".into();
            spectrum.peaks.push(Peak1D::new(999.0, 17.0));
        }
        self.spectra.push(spectrum.clone());
        Ok(self.flow())
    }
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.event(chromatogram.native_id.clone())?;
        if self.mutate {
            chromatogram.name = "changed".into();
        }
        self.chromatograms.push(chromatogram.clone());
        Ok(self.flow())
    }
}
fn cv(id: &str, value: &str) -> String {
    let unit = if matches!(id, "1000016" | "1000595") {
        " unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\""
    } else {
        ""
    };
    format!("<cvParam accession=\"MS:{id}\" name=\"fixture\" value=\"{value}\"{unit}/>")
}
fn array(role: &str, values: &[f64]) -> String {
    let bytes: Vec<_> = values.iter().flat_map(|n| n.to_le_bytes()).collect();
    let text = STANDARD.encode(bytes);
    format!(
        "<binaryDataArray encodedLength=\"{}\">{}{}{}<binary>{text}</binary></binaryDataArray>",
        text.len(),
        cv("1000523", ""),
        cv("1000576", ""),
        cv(role, "")
    )
}
fn spec(id: usize, rt: f64, level: u32, mz: &[f64], intensity: &[f64]) -> String {
    format!(
        "<spectrum id=\"s{id}\" index=\"{id}\" defaultArrayLength=\"{}\">{}<scanList count=\"1\"><scan>{}</scan></scanList><binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList></spectrum>",
        mz.len(),
        cv("1000511", &level.to_string()),
        cv("1000016", &rt.to_string()),
        array("1000514", mz),
        array("1000515", intensity)
    )
}
fn chrom(id: usize) -> String {
    format!(
        "<chromatogram id=\"c{id}\" index=\"{id}\" defaultArrayLength=\"1\"><binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList></chromatogram>",
        array("1000595", &[id as f64]),
        array("1000515", &[1.0])
    )
}
fn header() -> String {
    format!(
        "<fileDescription><fileContent/></fileDescription><softwareList count=\"1\"><software id=\"sw\" version=\"1\">{}</software></softwareList><instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"ic\"/></instrumentConfigurationList><dataProcessingList count=\"1\"><dataProcessing id=\"dp\"><processingMethod order=\"0\" softwareRef=\"sw\">{}</processingMethod></dataProcessing></dataProcessingList>",
        cv("1000799", "test"),
        cv("1000544", "")
    )
}
fn doc(spectra: &[String], chromatograms: &[String]) -> String {
    format!(
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">{}<run id=\"run\" defaultInstrumentConfigurationRef=\"ic\"><spectrumList count=\"{}\" defaultDataProcessingRef=\"dp\">{}</spectrumList><chromatogramList count=\"{}\" defaultDataProcessingRef=\"dp\">{}</chromatogramList></run></mzML>",
        header(),
        spectra.len(),
        spectra.concat(),
        chromatograms.len(),
        chromatograms.concat()
    )
}
fn grid(s: usize, c: usize) -> String {
    doc(
        &(0..s)
            .map(|i| spec(i, i as f64, 1, &[2.0, 1.0], &[4.0, 3.0]))
            .collect::<Vec<_>>(),
        &(0..c).map(chrom).collect::<Vec<_>>(),
    )
}
fn run(xml: &str, log: &mut Log, options: &TransformOptions) -> Result<mzml::TransformReport> {
    mzml::transform_from(|| Ok(Cursor::new(xml)), log, options)
}
fn into(
    xml: &str,
    log: &mut Log,
    destination: &mut MSExperiment,
    options: &TransformOptions,
) -> Result<mzml::TransformReport> {
    mzml::transform_from_into(|| Ok(Cursor::new(xml)), log, destination, options)
}

#[test]
fn literal_source_transform_counts_peaks_tic_and_both_retention_modes() {
    for retaining in [false, true] {
        let mut log = Log::default();
        let mut destination = MSExperiment::default();
        let options = TransformOptions {
            skip_first_pass: true,
            skip_full_count: true,
            ..Default::default()
        };
        let result = if retaining {
            into(SOURCE, &mut log, &mut destination, &options)
        } else {
            run(SOURCE, &mut log, &options)
        }
        .unwrap();
        assert_eq!(result.delivered.spectra, 4);
        assert_eq!(log.spectra.iter().map(|s| s.peaks.len()).sum::<usize>(), 40);
        assert_eq!(
            log.spectra
                .iter()
                .flat_map(|s| &s.peaks)
                .map(|p| f64::from(p.intensity))
                .sum::<f64>(),
            350.0
        );
        assert_eq!(
            log.chromatograms
                .iter()
                .map(|c| c.peaks.len())
                .collect::<Vec<_>>(),
            [15, 10]
        );
        assert_eq!(destination.spectra.len(), if retaining { 4 } else { 0 });
        assert!(log.settings.is_none());
    }
}
#[test]
fn all_setup_flags_use_raw_counts_and_exactly_one_or_two_factory_opens() {
    let xml = grid(3, 2);
    for skip_first_pass in [false, true] {
        for skip_full_count in [false, true] {
            let mut options = TransformOptions {
                skip_first_pass,
                skip_full_count,
                ..Default::default()
            };
            options.load.scientific.add_ms_level(2).unwrap();
            let mut opens = 0;
            let mut log = Log::default();
            let report = mzml::transform_from(
                || {
                    opens += 1;
                    Ok(Cursor::new(&xml))
                },
                &mut log,
                &options,
            )
            .unwrap();
            assert_eq!(opens, if skip_first_pass { 1 } else { 2 });
            assert_eq!(
                report.expected,
                if skip_first_pass {
                    None
                } else {
                    Some(MzMLCounts {
                        spectra: if skip_full_count { 0 } else { 3 },
                        chromatograms: if skip_full_count { 0 } else { 2 },
                    })
                }
            );
            assert!(log.spectra.is_empty());
            assert_eq!(log.chromatograms.len(), 2);
            if !skip_first_pass {
                assert!(log.events[0].starts_with("size:"));
                assert_eq!(log.events[1], "settings");
            }
        }
    }
}
#[test]
fn independent_pools_and_final_flush_match_source_order() {
    let xml = grid(3, 4);
    for capacity in [0, 1, 2, 100] {
        let mut options = TransformOptions::default();
        options.load.scientific.max_data_pool_size = capacity;
        let mut log = Log::default();
        run(&xml, &mut log, &options).unwrap();
        let expected = if capacity == 2 {
            vec!["s0", "s1", "c0", "c1", "c2", "c3", "s2"]
        } else {
            vec!["s0", "s1", "s2", "c0", "c1", "c2", "c3"]
        };
        assert_eq!(&log.events[2..], expected);
        assert_eq!(log.spectra[0].peaks[0].mz, 1.0);
    }
}
#[test]
fn entire_pool_is_decoded_before_its_first_callback() {
    let bad = spec(1, 1.0, 1, &[1.0], &[1.0]).replace("<binary>", "<binary>!");
    let xml = doc(&[spec(0, 0.0, 1, &[1.0], &[2.0]), bad], &[]);
    for capacity in [1, 2] {
        let mut options = TransformOptions::default();
        options.load.scientific.max_data_pool_size = capacity;
        let mut log = Log::default();
        let mut old = MSExperiment::default();
        old.spectra.push(MSSpectrum::default());
        let before = old.clone();
        assert!(into(&xml, &mut log, &mut old, &options).is_err());
        assert_eq!(old, before);
        assert_eq!(log.spectra.len(), if capacity == 1 { 1 } else { 0 });
    }
}
#[test]
fn soft_stop_excludes_current_record_and_ignores_pending_pool_and_tail() {
    let xml = grid(3, 2).replace("</chromatogramList>", "</chromatogramList><broken");
    let mut options = TransformOptions::default();
    options.load.scientific.max_data_pool_size = 2;
    let mut log = Log {
        stop: Some(2),
        ..Default::default()
    };
    let mut old = MSExperiment::default();
    old.spectra.push(MSSpectrum::default());
    let report = into(&xml, &mut log, &mut old, &options).unwrap();
    assert!(report.stopped);
    assert_eq!(
        report.delivered,
        MzMLCounts {
            spectra: 2,
            chromatograms: 0
        }
    );
    assert_eq!(old.spectra.len(), 2);
    assert_eq!(old.spectra[1].native_id, "s0");
    assert!(old.chromatograms.is_empty());
}
#[test]
fn setup_callback_errors_and_second_open_failure_stop_at_exact_boundary() {
    let xml = grid(1, 0);
    for fail in [1, 2, 3] {
        let mut opens = 0;
        let mut log = Log {
            fail: Some(fail),
            ..Default::default()
        };
        let mut old = MSExperiment::default();
        let result = mzml::transform_from_into(
            || {
                opens += 1;
                Ok(Cursor::new(&xml))
            },
            &mut log,
            &mut old,
            &TransformOptions::default(),
        );
        assert!(result.is_err());
        assert_eq!(log.events.len(), fail);
        assert_eq!(opens, if fail < 3 { 1 } else { 2 });
        assert!(old.spectra.is_empty());
    }
    let mut opens = 0;
    let mut log = Log::default();
    assert!(
        mzml::transform_from(
            || {
                opens += 1;
                if opens == 2 {
                    Err(Error::InvalidValue("second open".into()))
                } else {
                    Ok(Cursor::new(&xml))
                }
            },
            &mut log,
            &TransformOptions::default()
        )
        .is_err()
    );
    assert_eq!(log.events, ["size:1:0", "settings"]);
}
#[test]
fn retaining_preserves_old_peak_buffer_and_applies_callback_changes() {
    let mut old = MSExperiment::default();
    old.spectra.push(MSSpectrum {
        peaks: vec![Peak1D::new(9.0, 8.0)],
        ..Default::default()
    });
    let pointer = old.spectra[0].peaks.as_ptr();
    old.settings.comment = "prior".into();
    old.settings.document.loaded_file_path = "prior.raw".into();
    old.settings.metadata.insert("prior".into(), "kept".into());
    old.settings
        .metadata
        .insert("overwrite".into(), "old".into());
    let xml = grid(1, 1).replace(
        "<spectrumList",
        "<userParam name=\"overwrite\" type=\"xsd:string\" value=\"new\"/><spectrumList",
    );
    let mut log = Log {
        mutate: true,
        ..Default::default()
    };
    into(&xml, &mut log, &mut old, &TransformOptions::default()).unwrap();
    assert_eq!(old.spectra[0].peaks.as_ptr(), pointer);
    assert_eq!(old.spectra[1].name, "changed");
    assert_eq!(old.spectra[1].peaks.len(), 3);
    assert_eq!(old.chromatograms[0].name, "changed");
    assert_eq!(old.settings.comment, "prior");
    assert_eq!(old.settings.document.loaded_file_path, "prior.raw");
    assert_eq!(old.settings.metadata["prior"].as_str().unwrap(), "kept");
    assert_eq!(old.settings.metadata["overwrite"].as_str().unwrap(), "new");
    assert!(log.settings.unwrap().comment.is_empty());
}
#[test]
fn metadata_only_keeps_setup_raw_counts_and_never_delivers_records() {
    let mut options = TransformOptions::default();
    options.load.scientific.metadata_only = true;
    let mut log = Log::default();
    let report = run(&grid(2, 1), &mut log, &options).unwrap();
    assert_eq!(
        report.expected,
        Some(MzMLCounts {
            spectra: 2,
            chromatograms: 1
        })
    );
    assert_eq!(report.delivered, MzMLCounts::default());
    assert_eq!(log.events.len(), 2);
}
#[test]
fn fill_data_false_skips_decode_and_numeric_limits_but_retains_metadata() {
    let xml = grid(1, 1).replace("<binary>", "<binary>!").replace(
        "<spectrumList",
        "<userParam name=\"run-key\" value=\"v\"/><spectrumList",
    );
    let mut options = TransformOptions::default();
    options.load.scientific.fill_data = false;
    options.read.max_array_bytes = 0;
    options.read.max_total_peaks = 0;
    options.read.max_total_array_bytes = 0;
    options.read.max_total_array_elements = 0;
    let mut log = Log::default();
    let report = run(&xml, &mut log, &options).unwrap();
    assert_eq!(
        report.delivered,
        MzMLCounts {
            spectra: 1,
            chromatograms: 1
        }
    );
    assert!(log.spectra[0].peaks.is_empty());
    assert_eq!(log.spectra[0].ms_level, 1);
    assert!(log.chromatograms[0].peaks.is_empty());
    assert!(!log.spectra[0].data_processing.is_empty());
    let loaded =
        mzml::read_with_load_options(Cursor::new(&xml), &options.load, &options.read).unwrap();
    assert!(loaded.spectra[0].peaks.is_empty());
    assert_eq!(loaded.settings.metadata["run-key"].as_str().unwrap(), "v");
    assert!(run(&xml, &mut Log::default(), &TransformOptions::default()).is_err());
}
#[test]
fn skipping_chromatograms_preserves_header_and_spectra_cpp017() {
    let mut options = TransformOptions::default();
    options.load.scientific.skip_chromatograms = true;
    let mut log = Log::default();
    let report = run(&grid(2, 2), &mut log, &options).unwrap();
    assert_eq!(
        report.expected,
        Some(MzMLCounts {
            spectra: 2,
            chromatograms: 0
        })
    );
    assert_eq!(
        report.delivered,
        MzMLCounts {
            spectra: 2,
            chromatograms: 0
        }
    );
    assert!(log.settings.is_some());
}
#[test]
fn tiny_administrative_and_parser_budgets_preserve_destination() {
    let xml = grid(1, 1);
    for which in 0..4 {
        let mut options = TransformOptions::default();
        match which {
            0 => options.max_bytes = 0,
            1 => options.max_work = 0,
            2 => options.read.max_param_bytes = 0,
            _ => options.read.max_xml_bytes = 8,
        };
        let mut log = Log::default();
        let mut old = MSExperiment::default();
        old.settings.comment = "owned".repeat(1024);
        let before = old.clone();
        assert!(into(&xml, &mut log, &mut old, &options).is_err());
        assert_eq!(old, before);
    }
}

#[test]
fn chromatogram_soft_stop_leaves_pending_spectrum_undelivered() {
    let mut options = TransformOptions::default();
    options.load.scientific.max_data_pool_size = 2;
    let mut log = Log {
        stop: Some(4),
        ..Default::default()
    };
    let mut destination = MSExperiment::default();
    let report = into(&grid(3, 2), &mut log, &mut destination, &options).unwrap();
    assert!(report.stopped);
    assert_eq!(
        report.delivered,
        MzMLCounts {
            spectra: 2,
            chromatograms: 2
        }
    );
    assert_eq!(&log.events[2..], ["s0", "s1", "c0", "c1"]);
    assert_eq!(destination.spectra.len(), 2);
    assert_eq!(destination.chromatograms.len(), 1);
    assert_eq!(destination.chromatograms[0].native_id, "c0");
}
#[test]
fn failure_in_final_reservation_preserves_prior_contents_after_all_callbacks() {
    fn destination(capacity: usize) -> MSExperiment {
        let mut spectra = Vec::with_capacity(capacity);
        spectra.push(MSSpectrum {
            peaks: vec![Peak1D::new(7.0, 8.0)],
            ..Default::default()
        });
        let mut value = MSExperiment {
            spectra,
            ..Default::default()
        };
        value.settings.comment = "prior".into();
        value
    }
    let xml = grid(8, 0);
    let mut options = TransformOptions::default();
    options.load.scientific.max_data_pool_size = 1;
    // Find the administrative allowance required when the caller already owns
    // sufficient final capacity. No scientific numeric oracle is derived here.
    let (mut low, mut high) = (0, 200_000);
    while low < high {
        let middle = (low + high) / 2;
        options.max_bytes = middle;
        if into(&xml, &mut Log::default(), &mut destination(9), &options).is_ok() {
            high = middle
        } else {
            low = middle + 1
        }
    }
    assert!(low < 200_000);
    options.max_bytes = low;
    let mut old = destination(1);
    let before = old.clone();
    let pointer = old.spectra[0].peaks.as_ptr();
    let mut log = Log::default();
    assert!(into(&xml, &mut log, &mut old, &options).is_err());
    assert_eq!(log.spectra.len(), 8);
    assert_eq!(old, before);
    assert_eq!(old.spectra[0].peaks.as_ptr(), pointer);
}
#[test]
fn source_seeded_metadata_assignments_and_existing_source_file_dedup() {
    let original = mzml::read(Cursor::new(SOURCE)).unwrap();
    let mut old = MSExperiment {
        settings: original.settings.clone(),
        ..Default::default()
    };
    let files = old.settings.source_files.len();
    let contacts = old.settings.contacts.len();
    old.settings.sample.name = "old".into();
    old.settings.document.identifier = "old".into();
    old.settings.comment = "unrelated".into();
    old.settings.document.loaded_file_path = "preserved.raw".into();
    into(
        SOURCE,
        &mut Log::default(),
        &mut old,
        &TransformOptions::default(),
    )
    .unwrap();
    assert_eq!(old.settings.source_files.len(), files);
    assert_eq!(old.settings.contacts.len(), 2 * contacts);
    assert_eq!(old.settings.sample, original.settings.sample);
    assert_eq!(
        old.settings.document.identifier,
        original.settings.document.identifier
    );
    assert_eq!(old.settings.comment, "unrelated");
    assert_eq!(old.settings.document.loaded_file_path, "preserved.raw");
    assert_eq!(
        old.settings.instrument_configurations,
        original.settings.instrument_configurations
    );
}
#[test]
fn raw_f64_selection_and_auxiliary_permutation_precede_callback() {
    use openms::kernel::NumericRange;
    let mut s = spec(0, 4.0, 1, &[3.0, 1.0, 2.0], &[30.0, 1.0 + 1e-9, 20.0]);
    let extra = array("1000786", &[300.0, 100.0, 200.0]).replace("value=\"\"", "value=\"extra\"");
    // Only the nonstandard-array CV carries its name; precision/compression
    // value attributes are insignificant in this independent XML fixture.
    s = s
        .replace(
            "<binaryDataArrayList count=\"2\">",
            "<binaryDataArrayList count=\"3\">",
        )
        .replace(
            "</binaryDataArrayList>",
            &format!("{extra}</binaryDataArrayList>"),
        );
    let mut options = TransformOptions::default();
    options
        .load
        .scientific
        .set_mz_range(NumericRange { min: 1.0, max: 3.0 });
    options.load.scientific.set_intensity_range(NumericRange {
        min: 1.0 + 5e-10,
        max: 30.0,
    });
    let mut log = Log::default();
    run(&doc(&[s], &[]), &mut log, &options).unwrap();
    assert_eq!(
        log.spectra[0]
            .peaks
            .iter()
            .map(|p| p.mz)
            .collect::<Vec<_>>(),
        [1.0, 2.0]
    );
    assert_eq!(log.spectra[0].peaks[0].intensity, 1.0);
    assert_eq!(log.spectra[0].float_data_arrays[0].data, [100.0, 200.0]);
}
#[test]
fn signed_raw_count_setup_and_header_stop_do_not_count_filtered_records() {
    let first = grid(0, 0).replace("spectrumList count=\"0\"", "spectrumList count=\"-2\"");
    let second = grid(0, 0);
    let mut calls = 0;
    let mut log = Log::default();
    let result = mzml::transform_from(
        || {
            calls += 1;
            Ok(Cursor::new(if calls == 1 { &first } else { &second }))
        },
        &mut log,
        &TransformOptions::default(),
    )
    .unwrap();
    assert_eq!(result.expected, Some(MzMLCounts::default()));
    let mut options = TransformOptions {
        skip_full_count: true,
        ..Default::default()
    };
    options.load.scientific.metadata_only = true;
    let xml = grid(1, 1);
    let cut = xml.find("<spectrumList").unwrap();
    let tail = &xml[cut..];
    let end = cut + tail.find('>').unwrap() + 1;
    let prefix = &xml[..end];
    let mut log = Log::default();
    assert_eq!(
        run(prefix, &mut log, &options).unwrap().expected,
        Some(MzMLCounts::default())
    );
}
#[test]
fn no_population_validates_descriptors_and_identities_without_decoding_payload() {
    let good = grid(1, 0).replace("<binary>", "<binary>!");
    let mut options = TransformOptions::default();
    options.load.scientific.fill_data = false;
    assert!(run(&good, &mut Log::default(), &options).is_ok());
    let missing_type = good.replacen(&cv("1000514", ""), "", 1);
    let missing_precision = good.replacen(&cv("1000523", ""), "", 1);
    let missing_compression = good.replacen(&cv("1000576", ""), "", 1);
    let integer_primary = good.replacen(&cv("1000523", ""), &cv("1000519", ""), 1);
    let wrong_coordinate = good.replacen(&cv("1000514", ""), &cv("1000595", ""), 1);
    let duplicate_primary = good.replacen(&cv("1000515", ""), &cv("1000514", ""), 1);
    let wrong_length = good.replacen(
        "<binaryDataArray encodedLength",
        "<binaryDataArray arrayLength=\"9\" encodedLength",
        1,
    );
    let auxiliary = array("1000786", &[3.0, 4.0]).replace(
        "MS:1000786\" name=\"fixture\" value=\"\"",
        "MS:1000786\" name=\"fixture\" value=\"dup\"",
    );
    let duplicate_aux = good
        .replace(
            "<binaryDataArrayList count=\"2\">",
            "<binaryDataArrayList count=\"4\">",
        )
        .replace(
            "</binaryDataArrayList>",
            &format!("{auxiliary}{auxiliary}</binaryDataArrayList>"),
        );
    let numpress_integer = good
        .replacen(&cv("1000523", ""), &cv("1000519", ""), 1)
        .replacen(&cv("1000576", ""), &cv("1002312", ""), 1);
    for (name, xml) in [
        ("type", missing_type),
        ("precision", missing_precision),
        ("compression", missing_compression),
        ("integer primary", integer_primary),
        ("wrong coordinate", wrong_coordinate),
        ("duplicate primary", duplicate_primary),
        ("length", wrong_length),
        ("duplicate auxiliary", duplicate_aux),
        ("Numpress type", numpress_integer),
    ] {
        let error = run(&xml, &mut Log::default(), &options).unwrap_err();
        assert!(!error.to_string().contains("base64"), "{name}: {error}");
    }
}
#[test]
fn setup_resolves_header_groups_and_typed_run_groups_before_settings_callback() {
    let groups = format!(
        "<referenceableParamGroupList count=\"2\"><referenceableParamGroup id=\"header\"><userParam name=\"header-key\" value=\"retained\" type=\"xsd:string\"/></referenceableParamGroup><referenceableParamGroup id=\"run-group\">{}<userParam name=\"typed\" value=\"42\" type=\"xsd:integer\"/></referenceableParamGroup></referenceableParamGroupList>",
        cv("1000858", "fraction-7")
    );
    let xml=grid(0,0).replace("<softwareList",&format!("{groups}<softwareList"))
        .replace("<instrumentConfiguration id=\"ic\"/>","<instrumentConfiguration id=\"ic\"><referenceableParamGroupRef ref=\"header\"/></instrumentConfiguration>")
        .replace("<spectrumList", "<referenceableParamGroupRef ref=\"run-group\"/><spectrumList");
    let mut log = Log::default();
    let mut destination = MSExperiment::default();
    into(
        &xml,
        &mut log,
        &mut destination,
        &TransformOptions::default(),
    )
    .unwrap();
    let settings = log.settings.unwrap();
    assert_eq!(settings, destination.settings);
    assert_eq!(
        settings.instrument.metadata["header-key"].as_str().unwrap(),
        "retained"
    );
    assert_eq!(settings.metadata["typed"].as_i64().unwrap(), 42);
    assert_eq!(settings.fraction_identifier, "fraction-7");
}
#[test]
fn unused_primary_and_auxiliary_metadata_are_not_materialized_without_data() {
    let mut spectrum = spec(0, 1.0, 1, &[1.0], &[2.0]);
    let extra = array("1000786", &[3.0]).replace(
        "MS:1000786\" name=\"fixture\" value=\"\"",
        "MS:1000786\" name=\"fixture\" value=\"aux\"",
    );
    spectrum = spectrum
        .replace(
            "<binaryDataArrayList count=\"2\">",
            "<binaryDataArrayList count=\"3\">",
        )
        .replace(
            "</binaryDataArrayList>",
            &format!("{extra}</binaryDataArrayList>"),
        )
        .replace(
            "<binary>",
            "<userParam name=\"unused\" value=\"3.5\" type=\"xsd:double\"/><binary>!",
        );
    let xml = doc(&[spectrum], &[]);
    let mut options = TransformOptions::default();
    options.load.scientific.fill_data = false;
    let mut log = Log::default();
    run(&xml, &mut log, &options).unwrap();
    let s = &log.spectra[0];
    assert!(s.metadata.is_empty());
    assert!(s.float_data_arrays.is_empty());
    assert!(s.integer_data_arrays.is_empty());
    assert!(s.string_data_arrays.is_empty());
    assert_eq!(s.rt, 1.0);
    assert!(!s.data_processing.is_empty());
}
#[test]
fn numpress_source_reference_bytes_are_decoded_before_callbacks() {
    // Existing independently generated/pinned Numpress source golden values.
    let values = [100.0, 200.0, 300.00005, 400.00010];
    let mz = "QWR64UAAAADo//8/0P//f1kSgA==";
    let intensity = "ZGaMXCFQkQ==";
    let a = |role: &str, mode: &str, text: &str| {
        format!(
            "<binaryDataArray encodedLength=\"{}\">{}{}<binary>{text}</binary></binaryDataArray>",
            text.len(),
            cv(role, ""),
            cv(mode, "")
        )
    };
    let spectrum = format!(
        "<spectrum id=\"s0\" defaultArrayLength=\"4\"><binaryDataArrayList count=\"2\">{}{}</binaryDataArrayList></spectrum>",
        a("1000514", "1002312", mz),
        a("1000515", "1002313", intensity)
    );
    let mut log = Log::default();
    run(
        &doc(&[spectrum], &[]),
        &mut log,
        &TransformOptions::default(),
    )
    .unwrap();
    for (p, expected) in log.spectra[0].peaks.iter().zip(values) {
        assert!((p.mz - expected).abs() < 1e-6);
    }
    assert_eq!(
        log.spectra[0]
            .peaks
            .iter()
            .map(|p| p.intensity)
            .collect::<Vec<_>>(),
        [100.0, 200.0, 300.0, 400.0]
    );
}
#[test]
fn paths_use_magic_compression_and_streams_can_supply_one_byte_buffers() {
    use std::io::{BufReader, Write};
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Directory(std::path::PathBuf);
    impl Drop for Directory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let directory = Directory(std::env::temp_dir().join(format!(
        "openms-consumer-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    )));
    std::fs::create_dir(&directory.0).unwrap();
    let xml = grid(2, 1);
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(xml.as_bytes()).unwrap();
    let mut bzip = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bzip.write_all(xml.as_bytes()).unwrap();
    for (index, bytes) in [
        xml.as_bytes().to_vec(),
        gzip.finish().unwrap(),
        bzip.finish().unwrap(),
    ]
    .into_iter()
    .enumerate()
    {
        let path = directory.0.join(format!("file-{index}.data"));
        std::fs::write(&path, bytes).unwrap();
        let mut log = Log::default();
        let report = mzml::transform(&path, &mut log).unwrap();
        assert_eq!(
            report.delivered,
            MzMLCounts {
                spectra: 2,
                chromatograms: 1
            }
        );
        let mut destination = MSExperiment::default();
        mzml::transform_into(&path, &mut Log::default(), &mut destination).unwrap();
        assert_eq!(destination.spectra.len(), 2);
        assert!(destination.settings.document.loaded_file_path.is_empty());
    }
    let mut log = Log::default();
    mzml::transform_from(
        || Ok(BufReader::with_capacity(1, Cursor::new(&xml))),
        &mut log,
        &TransformOptions::default(),
    )
    .unwrap();
    assert_eq!(log.spectra.len(), 2);
}
#[test]
fn callback_numeric_mutations_are_retained_without_post_validation() {
    struct Mutator(Log);
    impl MSDataConsumer for Mutator {
        fn set_expected_size(&mut self, s: usize, c: usize) -> Result<()> {
            self.0.set_expected_size(s, c)
        }
        fn set_experimental_settings(&mut self, s: &ExperimentalSettings) -> Result<()> {
            self.0.set_experimental_settings(s)
        }
        fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
            s.peaks[0].intensity = f32::NAN;
            s.integer_data_arrays
                .push(openms::kernel::DataArray::new("unaligned", vec![1]));
            Ok(ControlFlow::Continue(()))
        }
        fn consume_chromatogram(&mut self, c: &mut MSChromatogram) -> Result<ControlFlow<()>> {
            self.0.consume_chromatogram(c)
        }
    }
    let xml = grid(1, 0);
    let mut destination = MSExperiment::default();
    let mut consumer = Mutator(Log::default());
    mzml::transform_from_into(
        || Ok(Cursor::new(&xml)),
        &mut consumer,
        &mut destination,
        &TransformOptions::default(),
    )
    .unwrap();
    assert!(destination.spectra[0].peaks[0].intensity.is_nan());
    assert_eq!(destination.spectra[0].integer_data_arrays[0].data, [1]);
    assert!(destination.validate().is_err());
}

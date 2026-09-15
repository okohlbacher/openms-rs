// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "mzml")]
use base64::{Engine, engine::general_purpose::STANDARD};
use openms::format::{
    indexed_mzml::IndexedMzMLDecoder,
    mzml::{self, PeakWriteLimits},
    peak_options::{NumpressCompression as Mode, NumpressConfig, PeakFileOptions},
};
use openms::kernel::{DataArray, Precursor};
use openms::{ChromatogramPeak, Error, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};
use std::{
    io::{self, Cursor, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::{AtomicU64, Ordering},
};
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static N: AtomicU64 = AtomicU64::new(0);
        let p = std::env::temp_dir().join(format!(
            "openms-writing-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn path(&self, n: &str) -> PathBuf {
        self.0.join(n)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn experiment() -> MSExperiment {
    MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            peaks: vec![
                Peak1D::new(100.123456789, 1.25),
                Peak1D::new(200.25, 2.5),
                Peak1D::new(-0.0, -0.0),
            ],
            ..Default::default()
        }],
        ..Default::default()
    }
}
fn output(e: &MSExperiment, o: &PeakFileOptions) -> Vec<u8> {
    let mut v = Vec::new();
    let r = mzml::write_with_peak_options(&mut v, e, o).unwrap();
    assert_eq!(r.xml_bytes, v.len() as u64);
    assert_eq!(r.indexed, o.write_index);
    v
}
fn arrays(raw: &[u8]) -> Vec<(&str, &str)> {
    std::str::from_utf8(raw)
        .unwrap()
        .split("<binaryDataArray encodedLength=")
        .skip(1)
        .map(|part| {
            let b = part
                .split_once("<binary>")
                .unwrap()
                .1
                .split_once("</binary>")
                .unwrap()
                .0;
            let pre = part.split_once("<binary>").unwrap().0;
            (pre, b)
        })
        .collect()
}
fn independent(raw: &[u8]) -> usize {
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("tools/mzml_writing/check_output.py");
    let mut child = Command::new("python3")
        .arg(script)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Python3 independent oracle");
    child.stdin.take().unwrap().write_all(raw).unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    String::from_utf8(result.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap()
}
#[test]
fn legacy_defaults_and_explicit_plain_source_options_keep_exact_output() {
    let e = experiment();
    // Plain output: `write_with_options` and `write_index = false`.
    let mut old = Vec::new();
    mzml::write_with_options(&mut old, &e, &Default::default()).unwrap();
    let mut o = PeakFileOptions::default();
    o.write_index = false;
    assert_eq!(output(&e, &o), old);
    let indexed = output(&e, &PeakFileOptions::default());
    assert!(indexed.starts_with(b"<?xml"));
    independent(&indexed);
    // `mzml::write` is indexed by default, as source `MzMLFile::store` with
    // default `PeakFileOptions` (`write_index_ = true`; the C++ Release tools
    // write `indexedmzML`, benchmark smoke run 2026-09-14). Its single pass
    // produces exactly the prepared two-pass writer's bytes.
    let mut streamed = Vec::new();
    mzml::write(&mut streamed, &e).unwrap();
    assert_eq!(streamed, indexed);
    assert_eq!(
        mzml::read(Cursor::new(old)).unwrap(),
        mzml::read(Cursor::new(indexed)).unwrap()
    );
}
#[test]
fn all_source_precision_pairs_match_independent_struct_pack_bytes() {
    let e = experiment();
    for row in include_str!("data/mzml_writing/precision.tsv")
        .lines()
        .skip(1)
    {
        let cells: Vec<_> = row.split('\t').collect();
        if cells[2] != "0" {
            continue;
        }
        let mut o = PeakFileOptions::default();
        o.write_index = false;
        o.mz_32_bit = cells[0] == "1";
        o.intensity_32_bit = cells[1] == "1";
        let raw = output(&e, &o);
        let a = arrays(&raw);
        assert_eq!(a[0].1, cells[5]);
        assert_eq!(a[1].1, cells[6]);
        for (i, width) in [(0, cells[3]), (1, cells[4])] {
            assert!(a[i].0.contains(if width == "4" {
                "MS:1000521"
            } else {
                "MS:1000523"
            }));
        }
    }
}
#[test]
fn mass_numpress_request_forces_both_primary_f64_fallbacks() {
    let e = experiment();
    for row in include_str!("data/mzml_writing/precision.tsv")
        .lines()
        .skip(1)
    {
        let cells: Vec<_> = row.split('\t').collect();
        if cells[2] != "1" {
            continue;
        }
        let mut o = PeakFileOptions::default();
        o.write_index = false;
        o.mz_32_bit = cells[0] == "1";
        o.intensity_32_bit = cells[1] == "1";
        // A coarse fixed linear factor cannot meet the source relative-error check.
        let _ = o.set_numpress_configuration_mass_time(NumpressConfig {
            compression: Mode::Linear,
            estimate_fixed_point: false,
            fixed_point: 0.001,
            error_tolerance: 1e-12,
            ..Default::default()
        });
        let mut raw = Vec::new();
        let report = mzml::write_with_peak_options(&mut raw, &e, &o).unwrap();
        let a = arrays(&raw);
        assert_eq!(report.binary.fallback_arrays, 1);
        assert_eq!(a[0].1, cells[5]);
        assert_eq!(a[1].1, cells[6]);
        assert!(a.iter().all(|p| p.0.contains("MS:1000523")));
    }
}
#[test]
fn original_linear_codec_bytes_and_legacy_fallback_precision_survive() {
    let e = MSExperiment {
        spectra: vec![MSSpectrum {
            native_id: "scan=1".into(),
            peaks: [100., 200., 300.00005, 400.00010]
                .into_iter()
                .map(|v| Peak1D::new(v, 1.0))
                .collect(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let config = NumpressConfig {
        compression: Mode::Linear,
        ..Default::default()
    };
    let mut o = PeakFileOptions::default();
    o.mz_32_bit = true;
    let _ = o.set_numpress_configuration_mass_time(config);
    let raw = output(&e, &o);
    let a = arrays(&raw);
    assert_eq!(a[0].1, "QWR64UAAAADo//8/0P//f1kSgA==");
    assert!(a[1].0.contains("MS:1000523"));
    independent(&raw);
    let mut legacy = Vec::new();
    mzml::write_with_numpress(
        &mut legacy,
        &e,
        &mzml::NumpressWriteOptions {
            mass_time: config,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(arrays(&legacy)[1].0.contains("MS:1000521"));
    assert!(!String::from_utf8(legacy).unwrap().contains("indexedmzML"));
}
#[test]
fn selected_ion_discriminator_is_consumed_by_every_writer() {
    let mut e = experiment();
    let mut p = Precursor {
        mz: 250.,
        isolation_target_mz: Some(250.),
        ..Default::default()
    };
    p.cv_terms
        .metadata
        .insert("selected ion m/z".into(), 999.0.try_into().unwrap());
    p.cv_terms
        .metadata
        .insert("confidence".into(), "kept".into());
    e.spectra[0].precursors.push(p.clone());
    e.chromatograms.push(MSChromatogram {
        precursor: p,
        ..Default::default()
    });
    for mode in 0..3 {
        let mut raw = Vec::new();
        match mode {
            0 => mzml::write(&mut raw, &e).unwrap(),
            1 => {
                mzml::write_with_numpress(&mut raw, &e, &Default::default()).unwrap();
            }
            _ => {
                mzml::write_with_peak_options(&mut raw, &e, &Default::default()).unwrap();
            }
        }
        let text = std::str::from_utf8(&raw).unwrap();
        assert!(!text.contains("<userParam name=\"selected ion m/z\""));
        assert!(text.contains("name=\"selected ion m/z\" value=\"999\""));
        let read = mzml::read(Cursor::new(raw)).unwrap();
        for p in [
            &read.spectra[0].precursors[0],
            &read.chromatograms[0].precursor,
        ] {
            assert_eq!(p.mz, 999.);
            assert_eq!(p.isolation_target_mz, Some(250.));
            assert_eq!(
                p.cv_terms
                    .metadata
                    .get("confidence")
                    .unwrap()
                    .as_str()
                    .unwrap(),
                "kept"
            );
        }
    }
}
#[test]
fn tpp_suppresses_only_isolation_and_forces_zero_charge_on_both_record_kinds() {
    let mut e = experiment();
    e.spectra[0].precursors.push(Precursor {
        mz: 500.,
        isolation_target_mz: Some(499.),
        isolation_window_lower_offset: 2.,
        charge: 0,
        ..Default::default()
    });
    e.chromatograms.push(MSChromatogram::default());
    let mut o = PeakFileOptions::default();
    o.force_tpp_compatibility = true;
    let raw = output(&e, &o);
    let text = std::str::from_utf8(&raw).unwrap();
    // Product retains its independent isolationWindow; precursor blocks do not.
    for block in text.split("<precursor>").skip(1) {
        let p = block.split_once("</precursor>").unwrap().0;
        assert!(!p.contains("isolationWindow"));
        assert!(p.contains("name=\"charge state\" value=\"0\""));
    }
    assert_eq!(text.matches("<precursor>").count(), 2);
    independent(&raw);
}
#[test]
fn bad_selected_ion_type_and_units_fail_before_any_writer_output() {
    let mut e = experiment();
    e.spectra[0].precursors.push(Precursor::default());
    let values = [
        "999".into(),
        openms::metadata::MetaValue::try_from(999.0)
            .unwrap()
            .with_unit(openms::metadata::Unit::new("MS:1000040", "m/z", "MS").unwrap())
            .unwrap(),
    ];
    for value in values {
        e.spectra[0].precursors[0]
            .cv_terms
            .metadata
            .insert("selected ion m/z".into(), value);
        for mode in 0..3 {
            let mut raw = b"old".to_vec();
            let r = match mode {
                0 => mzml::write(&mut raw, &e),
                1 => mzml::write_with_numpress(&mut raw, &e, &Default::default()).map(|_| ()),
                _ => mzml::write_with_peak_options(&mut raw, &e, &Default::default()).map(|_| ()),
            };
            assert!(r.is_err());
            assert_eq!(raw, b"old");
        }
    }
}

#[test]
fn index_rows_match_every_utf8_escaped_record_and_existing_decoder() {
    let directory = Directory::new();
    for (spectra, chroms) in [(true, false), (false, true), (true, true)] {
        let mut e = MSExperiment::default();
        if spectra {
            e.spectra = vec![
                MSSpectrum {
                    native_id: "scan=λ&\"1".into(),
                    ..Default::default()
                },
                MSSpectrum::default(),
            ];
        }
        if chroms {
            e.chromatograms = vec![
                MSChromatogram {
                    native_id: "κ & \"two\"".into(),
                    ..Default::default()
                },
                MSChromatogram::default(),
            ];
        }
        let bytes = output(&e, &Default::default());
        independent(&bytes);
        let path = directory.path("index.mzML");
        std::fs::write(&path, &bytes).unwrap();
        assert!(mzml::has_index(&path).unwrap());
        let decoder = IndexedMzMLDecoder::default();
        let offset = decoder.find_index_list_offset(&path).unwrap().unwrap();
        assert!(bytes[offset as usize..].starts_with(b"<indexList "));
        let index = decoder.parse_offsets(&path, offset).unwrap();
        assert_eq!(index.spectra.len(), e.spectra.len());
        assert_eq!(index.chromatograms.len(), e.chromatograms.len());
        for (id, offset) in index.spectra.iter().chain(&index.chromatograms) {
            assert!(
                std::str::from_utf8(&bytes[*offset as usize..])
                    .unwrap()
                    .starts_with(if id.starts_with("scan=") || id.starts_with("index=") {
                        "<spectrum "
                    } else {
                        "<chromatogram "
                    })
            );
        }
    }
}
struct ShortWriter {
    bytes: Vec<u8>,
    chunk: usize,
    fail_after: Option<usize>,
}
impl Write for ShortWriter {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self.fail_after == Some(self.bytes.len()) {
            return Err(io::Error::other("deliberate write failure"));
        }
        let n = b.len().min(self.chunk).min(
            self.fail_after
                .map_or(usize::MAX, |end| end - self.bytes.len()),
        );
        self.bytes.extend_from_slice(&b[..n]);
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
#[test]
fn real_sha1_covers_block_boundaries_and_only_successful_partial_writes() {
    let mut e = experiment();
    let mut found = std::collections::BTreeSet::new();
    for n in 0..64 {
        e.spectra[0]
            .metadata
            .insert("padding".into(), "a".repeat(n).into());
        let bytes = output(&e, &Default::default());
        let start = std::str::from_utf8(&bytes)
            .unwrap()
            .find("<fileChecksum>")
            .unwrap()
            + 14;
        let residue = start % 64;
        if [55, 56, 63, 0].contains(&residue) {
            assert_eq!(independent(&bytes), residue);
            found.insert(residue);
            for chunk in [1, 7, 63, 64] {
                let mut short = ShortWriter {
                    bytes: Vec::new(),
                    chunk,
                    fail_after: None,
                };
                mzml::write_with_peak_options(&mut short, &e, &Default::default()).unwrap();
                assert_eq!(short.bytes, bytes);
            }
        }
    }
    assert_eq!(found, [0, 55, 56, 63].into());
    let expected = output(&e, &Default::default());
    let end = expected.len() / 2;
    let mut writer = ShortWriter {
        bytes: Vec::new(),
        chunk: 7,
        fail_after: Some(end),
    };
    assert!(mzml::write_with_peak_options(&mut writer, &e, &Default::default()).is_err());
    assert_eq!(writer.bytes, expected[..end]);
}
#[test]
fn indexed_empty_rejects_before_path_and_external_output_cpp050() {
    let e = MSExperiment::default();
    let mut v = b"existing".to_vec();
    assert!(matches!(
        mzml::write_with_peak_options(&mut v, &e, &Default::default()),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(v, b"existing");
    struct NoPath;
    impl AsRef<Path> for NoPath {
        fn as_ref(&self) -> &Path {
            panic!("path requested before indexed-empty check")
        }
    }
    assert!(mzml::store_with_peak_options(NoPath, &e, &Default::default()).is_err());
    let mut plain = PeakFileOptions::default();
    plain.write_index = false;
    let mut old = Vec::new();
    mzml::write(&mut old, &e).unwrap();
    assert_eq!(output(&e, &plain), old);
}
#[test]
fn source_options_use_chromatogram_time_precision_and_auxiliary_configs() {
    let mut e = experiment();
    e.chromatograms.push(MSChromatogram {
        native_id: "c".into(),
        peaks: e.spectra[0]
            .peaks
            .iter()
            .map(|p| ChromatogramPeak::new(p.mz, p.intensity))
            .collect(),
        ..Default::default()
    });
    e.spectra[0]
        .float_data_arrays
        .push(DataArray::new("mean charge array", vec![1., 2., 3.]));
    e.spectra[0]
        .integer_data_arrays
        .push(DataArray::new("charge array", vec![1, 2, 3]));
    e.spectra[0].string_data_arrays.push(DataArray::new(
        "names",
        vec!["a".into(), "b".into(), "c".into()],
    ));
    let mut o = PeakFileOptions::default();
    o.mz_32_bit = true;
    o.zlib_compression = true;
    o.set_numpress_configuration_float_data_array(NumpressConfig {
        compression: Mode::Linear,
        ..Default::default()
    });
    let mut raw = Vec::new();
    let report = mzml::write_with_peak_options(&mut raw, &e, &o).unwrap();
    assert_eq!(report.binary.fallback_arrays, 1);
    let a = arrays(&raw);
    assert!(a[0].0.contains("MS:1000521"));
    assert!(a[5].0.contains("MS:1000521"));
    assert!(a[2].0.contains("MS:1000521"));
    assert!(a[3].0.contains("MS:1000519"));
    assert!(a[4].0.contains("MS:1001479"));
    assert!(a.iter().all(|x| x.0.contains("MS:1000574")));
    independent(&raw);
    let read = mzml::read(Cursor::new(raw)).unwrap();
    for (p, c) in read.spectra[0]
        .peaks
        .iter()
        .zip(&read.chromatograms[0].peaks)
    {
        assert_eq!(p.mz, c.rt);
    }
    assert_eq!(
        read.spectra[0].float_data_arrays,
        e.spectra[0].float_data_arrays
    );
}
#[test]
fn read_only_options_never_filter_sort_or_clone_ms_level_selection_on_store() {
    let e = experiment();
    let expected = output(&e, &Default::default());
    let mut o = PeakFileOptions::default();
    o.metadata_only = true;
    o.fill_data = false;
    o.skip_chromatograms = true;
    o.skip_xml_checks = true;
    o.always_append_data = true;
    o.sort_spectra_by_mz = false;
    o.sort_chromatograms_by_rt = false;
    o.force_mq_compatibility = true;
    o.write_supplemental_data = false;
    o.precursor_mz_selected_ion = false;
    o.max_data_pool_size = 0;
    o.set_rt_range(openms::kernel::NumericRange {
        min: f64::NAN,
        max: f64::NAN,
    });
    o.set_mz_range(openms::kernel::NumericRange {
        min: 1e9,
        max: 1e10,
    });
    for i in 0..10000 {
        o.add_ms_level(i).unwrap();
    }
    let pointer = o.ms_levels().as_ptr();
    assert_eq!(output(&e, &o), expected);
    assert_eq!(pointer, o.ms_levels().as_ptr());
}
#[test]
fn requested_f32_overflow_and_cumulative_limits_leave_output_and_path_unchanged() {
    let directory = Directory::new();
    let destination = directory.path("kept.mzML");
    std::fs::write(&destination, b"old file").unwrap();
    let mut e = experiment();
    e.spectra[0].peaks[0].mz = f64::MAX;
    let mut o = PeakFileOptions::default();
    o.mz_32_bit = true;
    let mut bytes = b"old".to_vec();
    assert!(mzml::write_with_peak_options(&mut bytes, &e, &o).is_err());
    assert_eq!(bytes, b"old");
    assert!(mzml::store_with_peak_options(&destination, &e, &o).is_err());
    assert_eq!(std::fs::read(&destination).unwrap(), b"old file");
    let e = experiment();
    let o = PeakFileOptions::default();
    let expected = output(&e, &o);
    let exact = PeakWriteLimits {
        max_xml_bytes: expected.len() as u64,
        ..Default::default()
    };
    let mut exact_bytes = Vec::new();
    mzml::write_with_peak_options_and_limits(&mut exact_bytes, &e, &o, &exact).unwrap();
    assert_eq!(exact_bytes, expected);
    for limits in [
        PeakWriteLimits {
            max_xml_bytes: expected.len() as u64 - 1,
            ..exact
        },
        PeakWriteLimits {
            max_work: 0,
            ..exact
        },
        PeakWriteLimits {
            max_bytes: 1,
            ..exact
        },
        PeakWriteLimits {
            max_xml_bytes: u64::MAX,
            ..exact
        },
    ] {
        let mut bytes = b"old".to_vec();
        assert!(mzml::write_with_peak_options_and_limits(&mut bytes, &e, &o, &limits).is_err());
        assert_eq!(bytes, b"old");
    }
}
#[test]
fn compressed_path_indexes_address_decoded_xml_and_publish_atomically() {
    use std::io::Read;
    let directory = Directory::new();
    let e = experiment();
    let expected = output(&e, &Default::default());
    for suffix in ["mzML", "mzML.gz", "mzML.bz2"] {
        let path = directory.path(&format!("file.{suffix}"));
        let report = mzml::store_with_peak_options(&path, &e, &Default::default()).unwrap();
        assert_eq!(report.xml_bytes, expected.len() as u64);
        let bytes = std::fs::read(&path).unwrap();
        let mut decoded = Vec::new();
        if suffix.ends_with(".gz") {
            flate2::read::GzDecoder::new(bytes.as_slice())
                .read_to_end(&mut decoded)
                .unwrap();
        } else if suffix.ends_with(".bz2") {
            bzip2::read::BzDecoder::new(bytes.as_slice())
                .read_to_end(&mut decoded)
                .unwrap();
        } else {
            decoded = bytes;
        }
        assert_eq!(decoded, expected);
        independent(&decoded);
        let mut expected_peaks = e.spectra[0].peaks.clone();
        expected_peaks.sort_by(|a, b| a.mz.total_cmp(&b.mz));
        assert_eq!(mzml::load(&path).unwrap().spectra[0].peaks, expected_peaks);
    }
}
#[test]
fn indexed_output_passes_actual_pinned_xsd_for_both_record_kinds_and_tpp() {
    if Command::new("xmllint").arg("--version").output().is_err() {
        eprintln!("xmllint unavailable; indexed XSD test not executed");
        return;
    }
    let mut e = experiment();
    e.chromatograms.push(MSChromatogram {
        native_id: "trace".into(),
        ..Default::default()
    });
    e.spectra[0].precursors.push(Precursor::default());
    for tpp in [false, true] {
        let mut o = PeakFileOptions::default();
        o.force_tpp_compatibility = tpp;
        let raw = output(&e, &o);
        let mut child = Command::new("xmllint")
            .args(["--nonet", "--noout", "--schema"])
            .arg(
                Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/data/mzml_writing/mzML_idx_1_10.xsd"),
            )
            .arg("-")
            .stdin(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&raw).unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        independent(&raw);
    }
}

#[test]
fn markup_and_index_allowances_accumulate_across_records_before_external_output() {
    let mut single = experiment();
    single.spectra[0]
        .metadata
        .insert("large & key".into(), "<&\"λ".repeat(4096).into());
    let mut repeated = single.clone();
    for i in 1..8 {
        let mut spectrum = single.spectra[0].clone();
        spectrum.native_id = format!("scan={}", i + 1);
        repeated.spectra.push(spectrum);
    }
    mzml::write_with_peak_options(std::io::sink(), &repeated, &Default::default()).unwrap();
    for byte_limit in [false, true] {
        let mut low = 0;
        let mut high = if byte_limit {
            PeakWriteLimits::default().max_bytes
        } else {
            PeakWriteLimits::default().max_work
        };
        // Find the public allowance accepted for one complete record, independently
        // of implementation constants; the same allowance must not reset per record.
        while low < high {
            let middle = low + (high - low) / 2;
            let mut limits = PeakWriteLimits::default();
            if byte_limit {
                limits.max_bytes = middle;
            } else {
                limits.max_work = middle;
            }
            if mzml::write_with_peak_options_and_limits(
                std::io::sink(),
                &single,
                &Default::default(),
                &limits,
            )
            .is_ok()
            {
                high = middle;
            } else {
                low = middle + 1;
            }
        }
        let mut limits = PeakWriteLimits::default();
        if byte_limit {
            limits.max_bytes = high;
        } else {
            limits.max_work = high;
        }
        mzml::write_with_peak_options_and_limits(
            std::io::sink(),
            &single,
            &Default::default(),
            &limits,
        )
        .unwrap();
        let mut untouched = b"existing output".to_vec();
        assert!(
            mzml::write_with_peak_options_and_limits(
                &mut untouched,
                &repeated,
                &Default::default(),
                &limits
            )
            .is_err()
        );
        assert_eq!(untouched, b"existing output");
    }
}

#[test]
fn three_source_codec_options_reproduce_original_bytes_with_and_without_zlib() {
    use std::io::Read;
    let literal = [100.0, 200.0, 300.00005, 400.00010];
    let mut spectrum = MSSpectrum::from_peaks(literal.map(|v| Peak1D::new(v, v as f32)).to_vec());
    spectrum.native_id = "scan=1".into();
    spectrum.float_data_arrays.push(DataArray::new(
        "signal to noise array",
        literal.map(|v| v as f32).to_vec(),
    ));
    let e = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    for zlib in [false, true] {
        let mut options = PeakFileOptions::default();
        options.zlib_compression = zlib;
        let _ = options.set_numpress_configuration_mass_time(NumpressConfig {
            compression: Mode::Linear,
            ..Default::default()
        });
        options.set_numpress_configuration_intensity(NumpressConfig {
            compression: Mode::Pic,
            ..Default::default()
        });
        options.set_numpress_configuration_float_data_array(NumpressConfig {
            compression: Mode::Slof,
            ..Default::default()
        });
        let mut raw = Vec::new();
        let report = mzml::write_with_peak_options(&mut raw, &e, &options).unwrap();
        assert_eq!(report.binary.encoded_arrays, 3);
        assert_eq!(report.binary.fallback_arrays, 0);
        for ((metadata, text), (plain, combined, original)) in arrays(&raw).into_iter().zip([
            ("MS:1002312", "MS:1002746", "QWR64UAAAADo//8/0P//f1kSgA=="),
            ("MS:1002313", "MS:1002747", "ZGaMXCFQkQ=="),
            ("MS:1002314", "MS:1002748", "QMVagAAAAAAZxX3ivPP8/w=="),
        ]) {
            assert!(metadata.contains(if zlib { combined } else { plain }));
            assert!(metadata.contains("MS:1000523"));
            let mut actual = STANDARD.decode(text).unwrap();
            if zlib {
                let mut decoded = Vec::new();
                flate2::read::ZlibDecoder::new(actual.as_slice())
                    .read_to_end(&mut decoded)
                    .unwrap();
                actual = decoded;
            }
            assert_eq!(actual, STANDARD.decode(original).unwrap());
        }
        independent(&raw);
        assert_eq!(
            mzml::read(Cursor::new(raw)).unwrap().spectra[0].peaks.len(),
            4
        );
    }
}

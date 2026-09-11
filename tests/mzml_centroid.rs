#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::{
    Error, MSExperiment, MSSpectrum, Peak1D,
    format::mzml::{self, CentroidInfoLimits, LoadOptions, ReadOptions, SpecInfo},
    kernel::{NumericRange, SpectrumType},
    metadata::{DataProcessing, ProcessingAction},
};
use std::{
    collections::BTreeMap,
    io::Write,
    path::{Path, PathBuf},
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
};
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "openms-centroid-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn spectrum(level: u32, kind: SpectrumType, points: usize) -> MSSpectrum {
    MSSpectrum {
        ms_level: level,
        spectrum_type: kind,
        peaks: (0..points)
            .map(|i| Peak1D::new(i as f64 + 1.0, 1.0))
            .collect(),
        ..Default::default()
    }
}
fn xml(mut spectra: Vec<MSSpectrum>) -> Vec<u8> {
    for (i, s) in spectra.iter_mut().enumerate() {
        s.native_id = format!("scan={i}");
    }
    let mut bytes = Vec::new();
    mzml::write(
        &mut bytes,
        &MSExperiment {
            spectra,
            ..Default::default()
        },
    )
    .unwrap();
    bytes
}
fn inspect(
    path: &Path,
    n: usize,
    options: &LoadOptions,
) -> openms::Result<BTreeMap<u32, SpecInfo>> {
    mzml::centroid_info_with_options(
        path,
        n,
        options,
        &ReadOptions::default(),
        CentroidInfoLimits::default(),
    )
}
fn counts(c: usize, p: usize, u: usize) -> SpecInfo {
    SpecInfo {
        count_centroided: c,
        count_profile: p,
        count_unknown: u,
    }
}
#[test]
fn global_recognized_quota_counts_unknowns_and_includes_last_stop_record() {
    let directory = Directory::new();
    let path = directory.write(
        "mixed.mzML",
        &xml(vec![
            spectrum(1, SpectrumType::Unknown, 0),
            spectrum(1, SpectrumType::Centroid, 0),
            spectrum(2, SpectrumType::Unknown, 4),
            spectrum(2, SpectrumType::Profile, 0),
            spectrum(3, SpectrumType::Centroid, 0),
        ]),
    );
    assert_eq!(
        inspect(&path, 1, &LoadOptions::default()).unwrap(),
        [(1, counts(1, 0, 1))].into()
    );
    assert_eq!(
        inspect(&path, 2, &LoadOptions::default()).unwrap(),
        [(1, counts(1, 0, 1)), (2, counts(0, 1, 1))].into()
    );
    let all: BTreeMap<_, _> = [
        (1, counts(1, 0, 1)),
        (2, counts(0, 1, 1)),
        (3, counts(1, 0, 0)),
    ]
    .into();
    assert_eq!(inspect(&path, 10, &LoadOptions::default()).unwrap(), all);
    assert_eq!(mzml::centroid_info(&path).unwrap(), all);
    let eleven = directory.write(
        "eleven.mzML",
        &xml((0..11)
            .map(|_| spectrum(1, SpectrumType::Centroid, 0))
            .collect()),
    );
    assert_eq!(
        mzml::centroid_info(&eleven).unwrap(),
        [(1, counts(10, 0, 0))].into()
    );
}
#[test]
fn stored_type_history_and_actual_estimator_precedence_survive_transport() {
    let picked = Arc::new(DataProcessing {
        actions: [ProcessingAction::PeakPicking].into(),
        ..Default::default()
    });
    let mut stored = spectrum(1, SpectrumType::Profile, 0);
    stored.data_processing.push(Arc::clone(&picked));
    let mut inferred = spectrum(2, SpectrumType::Unknown, 0);
    inferred.data_processing.push(picked);
    let mut data = spectrum(3, SpectrumType::Unknown, 0);
    data.peaks = [0.0, 2.0, 3.0, 10.0, 3.0, 2.0, 0.0]
        .into_iter()
        .enumerate()
        .map(|(i, y)| Peak1D::new(i as f64 * 0.1, y))
        .collect();
    let directory = Directory::new();
    let path = directory.write("precedence.mzML", &xml(vec![stored, inferred, data]));
    assert_eq!(
        inspect(&path, 3, &LoadOptions::default()).unwrap(),
        [
            (1, counts(0, 1, 0)),
            (2, counts(1, 0, 0)),
            (3, counts(0, 1, 0))
        ]
        .into()
    );
}
#[test]
fn filters_run_before_quota_and_metadata_only_returns_no_counts() {
    let mut specs = vec![
        spectrum(1, SpectrumType::Centroid, 0),
        spectrum(2, SpectrumType::Unknown, 0),
        spectrum(2, SpectrumType::Profile, 0),
        spectrum(2, SpectrumType::Centroid, 0),
    ];
    for (i, s) in specs.iter_mut().enumerate() {
        s.rt = i as f64;
    }
    let directory = Directory::new();
    let path = directory.write("filtered.mzML", &xml(specs));
    let mut options = LoadOptions::default();
    options.scientific.add_ms_level(2).unwrap();
    options
        .scientific
        .set_rt_range(NumericRange { min: 1.0, max: 3.0 });
    assert_eq!(
        inspect(&path, 1, &options).unwrap(),
        [(2, counts(0, 1, 1))].into()
    );
    options.scientific.metadata_only = true;
    assert!(inspect(&path, 1, &options).unwrap().is_empty());
    options.scientific.metadata_only = false;
    options.skip_spectra = true;
    assert!(inspect(&path, 1, &options).unwrap().is_empty());
}
#[test]
fn forced_population_preserves_all_caller_options_after_success_and_errors_cpp018() {
    let directory = Directory::new();
    let path = directory.write(
        "data.mzML",
        &xml(vec![spectrum(1, SpectrumType::Unknown, 6)]),
    );
    let mut options = LoadOptions::default();
    options.scientific.fill_data = false;
    options.scientific.max_data_pool_size = 1;
    options.scientific.skip_chromatograms = true;
    options.scientific.add_ms_level(1).unwrap();
    let before = options.scientific.clone();
    let pointer = options.scientific.ms_levels().as_ptr();
    assert_eq!(
        inspect(&path, 1, &options).unwrap(),
        [(1, counts(1, 0, 0))].into()
    );
    assert!(inspect(&directory.0.join("missing.mzML"), 1, &options).is_err());
    assert!(
        mzml::centroid_info_with_options(
            &path,
            1,
            &options,
            &ReadOptions::default(),
            CentroidInfoLimits {
                max_points: 5,
                ..Default::default()
            }
        )
        .is_err()
    );
    let invalid = directory.write("bad.mzML", b"not XML");
    assert!(inspect(&invalid, 1, &options).is_err());
    assert_eq!(options.scientific, before);
    assert_eq!(options.scientific.ms_levels().as_ptr(), pointer);
}
struct NoPath;
impl AsRef<Path> for NoPath {
    fn as_ref(&self) -> &Path {
        panic!("path touched before validation")
    }
}
#[test]
fn zero_quota_and_option_clone_limits_fail_before_path_access_cpp048() {
    assert!(matches!(
        mzml::centroid_info_with_options(
            NoPath,
            0,
            &LoadOptions::default(),
            &ReadOptions::default(),
            CentroidInfoLimits::default()
        ),
        Err(Error::InvalidValue(_))
    ));
    let mut options = LoadOptions::default();
    options.scientific.add_ms_level(1).unwrap();
    assert!(matches!(
        mzml::centroid_info_with_options(
            NoPath,
            1,
            &options,
            &ReadOptions::default(),
            CentroidInfoLimits {
                max_bytes: 3,
                ..Default::default()
            }
        ),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(options.scientific.ms_levels(), [1]);
}
#[test]
fn source_pool_predecode_can_fail_beyond_the_requested_quota() {
    let mut bytes = String::from_utf8(xml(vec![
        spectrum(1, SpectrumType::Centroid, 1),
        spectrum(1, SpectrumType::Centroid, 1),
    ]))
    .unwrap();
    let index = bytes.match_indices("<binary>").nth(2).unwrap().0 + "<binary>".len();
    bytes.insert(index, '!');
    let directory = Directory::new();
    let path = directory.write("bad-second.mzML", bytes.as_bytes());
    let mut options = LoadOptions::default();
    options.scientific.max_data_pool_size = 1;
    assert_eq!(
        inspect(&path, 1, &options).unwrap(),
        [(1, counts(1, 0, 0))].into()
    );
    options.scientific.max_data_pool_size = 2;
    assert!(inspect(&path, 1, &options).is_err());
}
#[test]
fn one_pass_soft_stop_does_not_touch_the_unread_xml_tail() {
    let bytes = String::from_utf8(xml(vec![spectrum(1, SpectrumType::Centroid, 0)]))
        .unwrap()
        .replace("</spectrumList>", "</spectrumList><broken");
    let directory = Directory::new();
    let path = directory.write("tail.mzML", bytes.as_bytes());
    let mut options = LoadOptions::default();
    options.scientific.max_data_pool_size = 1;
    assert_eq!(
        inspect(&path, 1, &options).unwrap(),
        [(1, counts(1, 0, 0))].into()
    );
    // A setup count pass would read through this sole list into the bad tail.
    assert!(mzml::load_size(&path).is_err());
}
#[test]
fn compressed_inputs_and_unsorted_option_follow_the_existing_consumer() {
    let mut s = spectrum(1, SpectrumType::Unknown, 0);
    s.peaks = [0.4, 0.1, 0.2, 0.3, 0.2, 0.5, 0.0]
        .into_iter()
        .zip([0.0, 2.0, 3.0, 10.0, 3.0, 2.0, 0.0])
        .map(|(x, y)| Peak1D::new(x, y))
        .collect();
    let bytes = xml(vec![s]);
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&bytes).unwrap();
    let mut bz = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    bz.write_all(&bytes).unwrap();
    let directory = Directory::new();
    let mut options = LoadOptions::default();
    options.scientific.sort_spectra_by_mz = false;
    for (i, data) in [bytes, gz.finish().unwrap(), bz.finish().unwrap()]
        .into_iter()
        .enumerate()
    {
        let path = directory.write(&format!("compressed-{i}.data"), &data);
        assert_eq!(
            inspect(&path, 1, &options).unwrap(),
            [(1, counts(0, 1, 0))].into()
        );
    }
}
#[test]
fn estimator_history_and_output_map_limits_accumulate_across_spectra() {
    let directory = Directory::new();
    let one = directory.write(
        "one.mzML",
        &xml(vec![spectrum(1, SpectrumType::Unknown, 6)]),
    );
    let two = directory.write(
        "two.mzML",
        &xml(vec![
            spectrum(1, SpectrumType::Unknown, 6),
            spectrum(1, SpectrumType::Unknown, 6),
        ]),
    );
    let options = LoadOptions::default();
    let read = ReadOptions::default();
    // Isolate a classification-only work threshold with one successful result;
    // replaying the same allowance over two records must fail cumulatively.
    let (mut low, mut high) = (0, 1000);
    while low < high {
        let middle = (low + high) / 2;
        let limit = CentroidInfoLimits {
            max_work: middle,
            ..Default::default()
        };
        if mzml::centroid_info_with_options(&one, 10, &options, &read, limit).is_ok() {
            high = middle
        } else {
            low = middle + 1
        }
    }
    assert!(low < 1000);
    assert!(
        mzml::centroid_info_with_options(
            &two,
            10,
            &options,
            &read,
            CentroidInfoLimits {
                max_work: low,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(
        mzml::centroid_info_with_options(
            &one,
            10,
            &options,
            &read,
            CentroidInfoLimits {
                max_bytes: 6 * 16,
                ..Default::default()
            }
        )
        .is_err()
    );
    let mut history = spectrum(1, SpectrumType::Unknown, 0);
    history.data_processing.push(Arc::new(DataProcessing {
        actions: [ProcessingAction::PeakPicking].into(),
        ..Default::default()
    }));
    let path = directory.write("history.mzML", &xml(vec![history]));
    assert!(
        mzml::centroid_info_with_options(
            &path,
            1,
            &options,
            &read,
            CentroidInfoLimits {
                max_work: 13,
                ..Default::default()
            }
        )
        .is_err()
    );
}

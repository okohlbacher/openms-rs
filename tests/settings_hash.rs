// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of the `std::hash` and `operator<<` sections of
//! `ChromatogramSettings_test.cpp` and `SpectrumSettings_test.cpp` (core SDK
//! bc9cc12), cited by source line. Literals are tier 3, transcribed from the
//! class tests.
//!
//! The settings records hold floating values and are therefore not `Eq`, so
//! they cannot key a `HashSet` directly; the class tests' `unordered_set` and
//! `unordered_map` checks run through [`Buckets`], which is exactly what a
//! hash container does: bucket by hash, then compare with `==`.

use std::collections::HashMap;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash};

use openms::kernel::{Precursor, SpectrumType};
use openms::metadata::{
    Acquisition, ChromatogramSettings, ChromatogramType, DataProcessing, ProcessingAction, Product,
    ScanMode, ScanWindow, Software, SpectrumSettings,
};

fn hash_of<T: Hash>(value: &T) -> u64 {
    BuildHasherDefault::<DefaultHasher>::default().hash_one(value)
}

/// A hash map keyed by a `Hash + PartialEq` type.
struct Buckets<K, V> {
    buckets: HashMap<u64, Vec<(K, V)>>,
}

impl<K: Hash + PartialEq + Clone, V: Copy> Buckets<K, V> {
    fn new() -> Self {
        Self {
            buckets: HashMap::new(),
        }
    }
    fn insert(&mut self, key: &K, value: V) {
        let bucket = self.buckets.entry(hash_of(key)).or_default();
        match bucket.iter_mut().find(|(k, _)| k == key) {
            Some(slot) => slot.1 = value,
            None => bucket.push((key.clone(), value)),
        }
    }
    fn get(&self, key: &K) -> Option<V> {
        self.buckets
            .get(&hash_of(key))
            .and_then(|bucket| bucket.iter().find(|(k, _)| k == key))
            .map(|(_, v)| *v)
    }
    fn len(&self) -> usize {
        self.buckets.values().map(Vec::len).sum()
    }
}

fn chromatogram(native_id: &str) -> ChromatogramSettings {
    ChromatogramSettings {
        native_id: native_id.into(),
        ..ChromatogramSettings::default()
    }
}

fn populated_chromatogram() -> ChromatogramSettings {
    let mut cs = ChromatogramSettings {
        native_id: "native_id_1".into(),
        comment: "test comment".into(),
        chromatogram_type: ChromatogramType::SelectedReactionMonitoring,
        ..ChromatogramSettings::default()
    };
    cs.precursor.mz = 500.5;
    cs.product.mz = 200.2;
    cs.acquisition_info.method_of_combination = "sum".into();
    cs.instrument_settings.scan_windows.push(ScanWindow {
        begin: 100.0,
        end: 1000.0,
        ..ScanWindow::default()
    });
    cs
}

// ChromatogramSettings_test.cpp L427 [EXTRA] std::hash<ChromatogramSettings>
#[test]
fn chromatogram_section_extra_hash() {
    let cs1 = ChromatogramSettings::default();
    let cs2 = ChromatogramSettings::default();
    assert!(cs1 == cs2);
    assert_eq!(hash_of(&cs1), hash_of(&cs2));

    let cs1 = populated_chromatogram();
    let cs2 = populated_chromatogram();
    assert!(cs1 == cs2);
    assert_eq!(hash_of(&cs1), hash_of(&cs2));

    let cs3 = chromatogram("different_id");
    assert!(cs1 != cs3);
    assert_ne!(hash_of(&cs1), hash_of(&cs3));

    let mut set = Buckets::new();
    let (s1, s2, s3) = (
        chromatogram("id1"),
        chromatogram("id2"),
        chromatogram("id1"),
    );
    set.insert(&s1, ());
    set.insert(&s2, ());
    assert_eq!(set.len(), 2);
    set.insert(&s3, ());
    assert_eq!(set.len(), 2);
    assert!(set.get(&s1).is_some());
    assert!(set.get(&s2).is_some());
    assert!(set.get(&s3).is_some());

    let mut map = Buckets::new();
    let (k1, k2, k3) = (
        chromatogram("key1"),
        chromatogram("key2"),
        chromatogram("key1"),
    );
    map.insert(&k1, 1);
    map.insert(&k2, 2);
    assert_eq!(map.len(), 2);
    assert_eq!(map.get(&k1), Some(1));
    assert_eq!(map.get(&k2), Some(2));
    assert_eq!(map.get(&k3), Some(1));
    map.insert(&k3, 3);
    assert_eq!(map.get(&k1), Some(3));

    assert_eq!(hash_of(&cs1), hash_of(&cs1));
}

/// `operator<<` (ChromatogramSettings.cpp:149-154, SpectrumSettings.cpp:185-190)
/// writes two banner lines and nothing of the record; derived from the source,
/// which has no class-test section for it.
#[test]
fn stream_output_is_the_two_banner_lines() {
    assert_eq!(
        populated_chromatogram().to_string(),
        "-- CHROMATOGRAMSETTINGS BEGIN --\n-- CHROMATOGRAMSETTINGS END --\n"
    );
    assert_eq!(
        SpectrumSettings::default().to_string(),
        "-- SPECTRUMSETTINGS BEGIN --\n-- SPECTRUMSETTINGS END --\n"
    );
}

fn populated_spectrum() -> SpectrumSettings {
    let mut s1 = SpectrumSettings {
        native_id: "scan=1".into(),
        comment: "test comment".into(),
        spectrum_type: SpectrumType::Centroid,
        ..SpectrumSettings::default()
    };
    s1.instrument_settings.scan_mode = ScanMode::MsnSpectrum;
    s1.instrument_settings.zoom_scan = true;
    s1.instrument_settings.scan_windows.push(ScanWindow {
        begin: 100.0,
        end: 2000.0,
        ..ScanWindow::default()
    });
    s1.acquisition_info.method_of_combination = "sum".into();
    s1.acquisition_info.acquisitions.push(Acquisition {
        identifier: "acq_1".into(),
        ..Acquisition::default()
    });
    s1.source_file.name = "test.mzML".into();
    s1.source_file.path = "/data/".into();
    s1.source_file.file_type = "mzML".into();
    let precursor = Precursor {
        mz: 500.5,
        intensity: 1000.0,
        charge: 2,
        activation_energy: 35.0,
        isolation_window_lower_offset: 1.5,
        isolation_window_upper_offset: 1.5,
        ..Precursor::default()
    };
    s1.precursors.push(precursor.into());
    s1.products.push(Product {
        mz: 250.25,
        ..Product::default()
    });
    let mut processing = DataProcessing {
        software: Software::new("TestTool", "1.0"),
        ..DataProcessing::default()
    };
    processing.actions.insert(ProcessingAction::PeakPicking);
    s1.data_processing.push(processing);
    s1.metadata.insert("key1".into(), "value1".into());
    s1.metadata.insert("key2".into(), 42i32.into());
    s1
}

// SpectrumSettings_test.cpp L435 [EXTRA] std::hash<SpectrumSettings>
#[test]
fn spectrum_section_extra_hash() {
    let s1 = populated_spectrum();
    let s2 = s1.clone();
    assert!(s1 == s2);
    assert_eq!(hash_of(&s1), hash_of(&s2));

    let empty1 = SpectrumSettings::default();
    let empty2 = SpectrumSettings::default();
    assert_eq!(hash_of(&empty1), hash_of(&empty2));

    let diff = SpectrumSettings {
        native_id: "scan=999".into(),
        comment: "different comment".into(),
        spectrum_type: SpectrumType::Profile,
        ..SpectrumSettings::default()
    };
    assert_ne!(hash_of(&s1), 0);

    let mut set = Buckets::new();
    set.insert(&s1, ());
    set.insert(&empty1, ());
    set.insert(&diff, ());
    assert_eq!(set.len(), 3);
    set.insert(&s2, ());
    assert_eq!(set.len(), 3);

    let mut map = Buckets::new();
    map.insert(&s1, "first");
    map.insert(&empty1, "empty");
    map.insert(&diff, "different");
    assert_eq!(map.len(), 3);
    assert_eq!(map.get(&s1), Some("first"));
    assert_eq!(map.get(&s2), Some("first"));

    let mut modified = s1.clone();
    modified.comment = "modified comment".into();
    assert!(s1 != modified);
    // The class test only computes both hashes; they differ here too.
    assert_ne!(hash_of(&s1), hash_of(&modified));
}

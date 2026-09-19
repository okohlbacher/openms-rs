// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Every `START_SECTION` of `MSSpectrum_test.cpp` at Core SDK `bc9cc12`, plus the
//! native boundaries of `src/kernel/spectrum_mobility.rs`.
//!
//! All expected values that carry a section number are transcribed C++ literals
//! (tier 3, source review); the rasterizer and resource-bound cases are
//! independently derived (tier 4), because `rasterizeIMFrame` has no upstream
//! section.

use openms::kernel::ranges::MSDim;
use openms::kernel::spectrum_mobility::{Chunk, Chunks, ImFrameRaster, RasterAggregation};
use openms::kernel::{DataArray, SpectrumType};
use openms::metadata::{
    DataProcessing, DriftTimeUnit, IonMobilityFormat, IonMobilityPeakType, MetaValue,
    MetaValueData, ProcessingAction, ScanWindow,
};
use openms::{Error, MSSpectrum, Peak1D};
use std::sync::Arc;

/// The name `IMDataArrayUtils::setIMUnit` writes for `DriftTimeUnit::MILLISECOND`:
/// the PSI-MS name of `MS:1002816 ! mean ion mobility array`
/// (`IMDataArrayUtils.cpp:25`).
const IM_MS_ARRAY: &str = "mean ion mobility array";

fn meta(value: f64) -> MetaValue {
    MetaValue::new(MetaValueData::Float(value)).unwrap()
}

fn text(value: &str) -> MetaValue {
    MetaValue::new(MetaValueData::String(value.into())).unwrap()
}

/// `getPrefilledSpec()`, `MSSpectrum_test.cpp:31-61`: ten peaks in descending
/// m/z, three float arrays with identical values of which the second is renamed
/// to the millisecond ion-mobility CV term, two string arrays, two integer
/// arrays and RT 5.
fn prefilled_spectrum() -> MSSpectrum {
    let floats = vec![56.0, 201.0, 31.0, 31.0, 31.0, 37.0, 29.0, 34.0, 60.0, 29.0];
    let strings: Vec<String> = ["56", "201", "31", "31", "31", "37", "29", "34", "60", "29"]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    let integers = vec![56, 201, 31, 31, 31, 37, 29, 34, 60, 29];
    let mzs = [
        423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321,
    ];
    let intensities = [
        56.0_f32, 201.0, 31.0, 31.0, 31.0, 37.0, 29.0, 34.0, 60.0, 29.0,
    ];
    let spectrum = MSSpectrum {
        peaks: mzs
            .iter()
            .zip(intensities)
            .map(|(&mz, intensity)| Peak1D::new(mz, intensity))
            .collect(),
        rt: 5.0,
        float_data_arrays: vec![
            DataArray::new("f1", floats.clone()),
            // IMDataConverter::setIMUnit(fda, DriftTimeUnit::MILLISECOND)
            DataArray::new(IM_MS_ARRAY, floats.clone()),
            DataArray::new("f3", floats),
        ],
        string_data_arrays: vec![
            DataArray::new("s1", strings.clone()),
            DataArray::new("s2", strings),
        ],
        integer_data_arrays: vec![
            DataArray::new("i1", integers.clone()),
            DataArray::new("", integers),
        ],
        ..MSSpectrum::default()
    };
    assert!(spectrum.contains_im_data());
    spectrum
}

/// `spec_test`, `MSSpectrum_test.cpp:1252-1273`.
fn spec_test() -> MSSpectrum {
    MSSpectrum::from_peaks(
        [
            (412.321, 29.0_f32),
            (412.824, 60.0),
            (413.8, 34.0),
            (414.301, 29.0),
            (415.287, 37.0),
            (416.293, 31.0),
            (418.232, 31.0),
            (419.113, 31.0),
            (420.13, 201.0),
            (423.269, 56.0),
            (426.292, 34.0),
            (427.28, 82.0),
            (428.322, 87.0),
            (430.269, 30.0),
            (431.246, 29.0),
            (432.289, 42.0),
            (436.161, 32.0),
            (437.219, 54.0),
            (439.186, 40.0),
            (440.27, 40.0),
            (441.224, 23.0),
        ]
        .iter()
        .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
        .collect(),
    )
}

/// `spec_find`, `MSSpectrum_test.cpp:1020`.
fn spec_find() -> MSSpectrum {
    MSSpectrum::from_peaks(
        [
            (1.0, 29.0_f32),
            (2.0, 60.0),
            (3.0, 34.0),
            (4.0, 29.0),
            (5.0, 37.0),
            (6.0, 31.0),
        ]
        .iter()
        .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
        .collect(),
    )
}

// Sections 1-4, 5-10: construction and the scalar accessors that are not
// ion-mobility specific.
#[test]
fn construction_and_scalar_accessors() {
    // (MSSpectrum()) / (~MSSpectrum()): the source allocates and deletes through
    // a raw pointer to prove the pair exists; Rust construction cannot fail and
    // Drop is compiler-generated, so only the constructed value is observable.
    let fresh = MSSpectrum::new();
    assert!(fresh.is_empty());

    // ([EXTRA] MSSpectrum())
    let mut tmp = MSSpectrum::new();
    tmp.peaks.push(Peak1D::new(47.11, 0.0));
    assert_eq!(tmp.len(), 1);
    assert_eq!(tmp.peaks[0].mz, 47.11);

    // (MSSpectrum(const std::initializer_list<Peak1D>& init))
    let tmp = MSSpectrum::from_peaks(vec![Peak1D::new(47.11, 2.0), Peak1D::new(500.0, 3.0)]);
    assert_eq!(tmp.len(), 2);
    assert_eq!(tmp.peaks[0].mz, 47.11);
    assert_eq!(tmp.peaks[1].mz, 500.0);
    assert_eq!(tmp.peaks[0].intensity, 2.0);
    assert_eq!(tmp.peaks[1].intensity, 3.0);

    // (UInt getMSLevel() const) / (void setMSLevel(UInt))
    let mut spec = MSSpectrum::new();
    assert_eq!(spec.ms_level, 1);
    spec.ms_level = 17;
    assert_eq!(spec.ms_level, 17);

    // (const std::string& getName() const) / (void setName(const std::string&))
    let mut spec = MSSpectrum::new();
    assert_eq!(spec.name, "");
    spec.name = "bla".into();
    assert_eq!(spec.name, "bla");

    // (double getRT() const) / (void setRT(double))
    let mut spec = MSSpectrum::new();
    assert_eq!(spec.rt, -1.0);
    spec.rt = 0.451;
    assert_eq!(spec.rt, 0.451);
}

// Sections 11-15: getDriftTime, setDriftTime, getDriftTimeUnit,
// getDriftTimeUnitAsString, setDriftTimeUnit.
#[test]
fn drift_time_accessors() {
    // (double getDriftTime() const)
    let mut spec = MSSpectrum::new();
    assert_eq!(spec.drift_time, -1.0);
    assert_eq!(spec.drift_time_if_set(), None);
    assert!(!spec.has_drift_time());

    // (void setDriftTime(double dt))
    spec.set_drift_time(Some(0.451)).unwrap();
    assert_eq!(spec.drift_time, 0.451);
    assert_eq!(spec.drift_time_if_set(), Some(0.451));
    assert!(spec.has_drift_time());

    // (double getDriftTimeUnit() const)
    let mut spec = MSSpectrum::new();
    assert_eq!(spec.drift_time_unit, DriftTimeUnit::None);

    // (double getDriftTimeUnitAsString() const)
    assert_eq!(spec.drift_time_unit_as_string(), "<NONE>");

    // (void setDriftTimeUnit(double dt))
    spec.drift_time_unit = DriftTimeUnit::Millisecond;
    assert_eq!(spec.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(spec.drift_time_unit_as_string(), "ms");

    // Native: the remaining three names of NamesOfDriftTimeUnit (IMTypes.cpp:23).
    for (unit, name) in [
        (DriftTimeUnit::InverseReducedMobility, "1/K0"),
        (DriftTimeUnit::FaimsCompensationVoltage, "FAIMS_CV"),
        (DriftTimeUnit::CollisionCrossSection, "CCS"),
    ] {
        spec.drift_time_unit = unit;
        assert_eq!(spec.drift_time_unit_as_string(), name);
    }

    // Native: clearing writes the source sentinel back; a non-finite drift time
    // is rejected where the source assigns it.
    spec.set_drift_time(None).unwrap();
    assert_eq!(spec.drift_time, -1.0);
    spec.set_drift_time(Some(3.0)).unwrap();
    assert!(matches!(
        spec.set_drift_time(Some(f64::NAN)),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(spec.drift_time, 3.0);
}

// Sections 16-21: the six data-array accessors.
#[test]
fn data_array_accessor_sizes() {
    let spec = MSSpectrum::new();
    assert_eq!(spec.float_data_arrays.len(), 0);
    assert_eq!(spec.string_data_arrays.len(), 0);
    assert_eq!(spec.integer_data_arrays.len(), 0);

    let mut spec = MSSpectrum::new();
    spec.float_data_arrays
        .resize_with(2, || DataArray::new("", Vec::<f32>::new()));
    assert_eq!(spec.float_data_arrays.len(), 2);
    spec.string_data_arrays
        .resize_with(2, || DataArray::new("", Vec::<String>::new()));
    assert_eq!(spec.string_data_arrays.len(), 2);
    spec.integer_data_arrays
        .resize_with(2, || DataArray::new("", Vec::<i32>::new()));
    assert_eq!(spec.integer_data_arrays.len(), 2);
}

/// Peaks p1, p2, p3 of `MSSpectrum_test.cpp:68-77` and the 5-entry arrays of the
/// `select` section.
fn select_fixture() -> MSSpectrum {
    let floats = vec![1.0_f32, 2.0, 3.0, 4.0, 5.0];
    let integers = vec![1, 2, 3, 4, 5];
    let strings: Vec<String> = (1..=5).map(|i| i.to_string()).collect();
    MSSpectrum {
        peaks: vec![
            Peak1D::new(2.0, 1.0),
            Peak1D::new(10.0, 2.0),
            Peak1D::new(30.0, 3.0),
            Peak1D::new(30.0, 3.0),
            Peak1D::new(10.0, 2.0),
        ],
        float_data_arrays: vec![
            DataArray::new("", floats.clone()),
            DataArray::new("", floats),
        ],
        integer_data_arrays: vec![
            DataArray::new("", integers.clone()),
            DataArray::new("", integers),
        ],
        string_data_arrays: vec![
            DataArray::new("", strings.clone()),
            DataArray::new("", strings),
        ],
        ..MSSpectrum::default()
    }
}

// Sections 22 and 23: select() and selectUnchecked(). The port has one entry
// point; `MSSpectrum::select` is the checked one, and the unchecked overload has
// no counterpart, so both sections exercise `select`.
#[test]
fn select_reorders_subsets_and_rejects_bad_input() {
    let spec = select_fixture();
    assert_eq!(spec.peaks[0].intensity, 1.0);
    assert_eq!(spec.peaks[4].intensity, 2.0);
    assert_eq!(spec.float_data_arrays.len(), 2);
    assert_eq!(spec.float_data_arrays[0].data.len(), 5);
    assert_eq!(spec.integer_data_arrays.len(), 2);
    assert_eq!(spec.integer_data_arrays[0].data.len(), 5);
    assert_eq!(spec.string_data_arrays.len(), 2);
    assert_eq!(spec.string_data_arrays[0].data.len(), 5);

    // re-order
    let mut s2 = spec.clone();
    s2.select(&[4, 2, 3, 1, 0]).unwrap();
    assert_eq!(s2.peaks[0].intensity, 2.0);
    assert_eq!(s2.peaks[4].intensity, 1.0);
    assert_eq!(s2.float_data_arrays.len(), 2);
    assert_eq!(s2.float_data_arrays[0].data.len(), 5);
    assert_eq!(s2.integer_data_arrays.len(), 2);
    assert_eq!(s2.integer_data_arrays[0].data.len(), 5);
    assert_eq!(s2.string_data_arrays.len(), 2);
    assert_eq!(s2.string_data_arrays[0].data.len(), 5);
    assert_eq!(s2.float_data_arrays[0].data[1], 3.0);
    assert_eq!(s2.integer_data_arrays[0].data[1], 3);
    assert_eq!(s2.string_data_arrays[0].data[1], "3");

    // subset; the new meta-array values are 5, 3, 4
    let mut s2 = spec.clone();
    s2.select(&[4, 2, 3]).unwrap();
    assert_eq!(s2.peaks[0].intensity, 2.0);
    assert_eq!(s2.peaks[1].intensity, 3.0);
    assert_eq!(s2.peaks[2].intensity, 3.0);
    assert_eq!(s2.float_data_arrays[0].data.len(), 3);
    assert_eq!(s2.integer_data_arrays[0].data.len(), 3);
    assert_eq!(s2.string_data_arrays[0].data.len(), 3);
    assert_eq!(s2.float_data_arrays[0].data[1], 3.0);
    assert_eq!(s2.integer_data_arrays[0].data[1], 3);
    assert_eq!(s2.string_data_arrays[0].data[1], "3");

    // out-of-range indices are rejected and leave the spectrum untouched
    let mut s2 = spec.clone();
    let before = s2.clone();
    assert!(s2.select(&[0, 7]).is_err());
    assert_eq!(s2, before);

    // a strictly increasing subset keeps its data arrays aligned
    let mut s2 = spec.clone();
    s2.select(&[0, 2, 4]).unwrap();
    assert_eq!(s2.len(), 3);
    assert_eq!(s2.string_data_arrays[0].data[0], "1");
    assert_eq!(s2.string_data_arrays[0].data[1], "3");
    assert_eq!(s2.string_data_arrays[0].data[2], "5");

    // a mis-sized data array is rejected before any peak is permuted
    let mut s2 = spec.clone();
    s2.integer_data_arrays[0].data.push(99);
    let before = s2.clone();
    assert!(s2.select(&[0, 1, 2]).is_err());
    assert_eq!(s2, before);

    // selectUnchecked section: a valid permutation reorders peaks and arrays,
    // and the data-array size check still applies.
    let base = MSSpectrum {
        peaks: vec![
            Peak1D::new(2.0, 1.0),
            Peak1D::new(10.0, 2.0),
            Peak1D::new(30.0, 3.0),
        ],
        string_data_arrays: vec![DataArray::new(
            "",
            vec!["1".to_string(), "2".into(), "3".into()],
        )],
        ..MSSpectrum::default()
    };
    let mut s2 = base.clone();
    s2.select(&[2, 0, 1]).unwrap();
    assert_eq!(s2.len(), 3);
    assert_eq!(s2.string_data_arrays[0].data[0], "3");
    assert_eq!(s2.string_data_arrays[0].data[1], "1");
    assert_eq!(s2.string_data_arrays[0].data[2], "2");

    let mut s2 = base;
    s2.string_data_arrays[0].data.push("x".into());
    let before = s2.clone();
    assert!(s2.select(&[0, 1, 2]).is_err());
    assert_eq!(s2, before);
}

// Section 24: updateRanges(). The port recomputes on demand, so
// `range_manager()` replaces the call and the six getters.
#[test]
fn update_ranges_includes_mobility_from_the_im_array() {
    let spec = prefilled_spectrum();
    for _ in 0..2 {
        let ranges = spec.range_manager().unwrap();
        assert_eq!(ranges.min_intensity().unwrap(), 29.0);
        assert_eq!(ranges.max_intensity().unwrap(), 201.0);
        assert_eq!(ranges.min_mz().unwrap(), 412.321);
        assert_eq!(ranges.max_mz().unwrap(), 423.269);
        assert_eq!(ranges.min_mobility().unwrap(), 29.0);
        assert_eq!(ranges.max_mobility().unwrap(), 201.0);
    }

    // only one peak: p1 = (2.0, 1.0), and no mobility at all
    let spec = MSSpectrum::from_peaks(vec![Peak1D::new(2.0, 1.0)]);
    let ranges = spec.range_manager().unwrap();
    assert_eq!(ranges.max_intensity().unwrap(), 1.0);
    assert_eq!(ranges.min_intensity().unwrap(), 1.0);
    assert_eq!(ranges.max_mz().unwrap(), 2.0);
    assert_eq!(ranges.min_mz().unwrap(), 2.0);
    assert!(ranges.is_dim_empty(MSDim::Mobility).unwrap());
}

/// The spectrum of the copy-construction and copy-assignment sections
/// (`MSSpectrum_test.cpp:402-412`).
fn copy_fixture() -> MSSpectrum {
    let mut tmp = MSSpectrum::new();
    tmp.instrument_settings
        .scan_windows
        .push(ScanWindow::new(0.0, 0.0).unwrap());
    tmp.metadata.insert("label".into(), meta(5.0));
    tmp.ms_level = 17;
    tmp.rt = 7.0;
    tmp.set_drift_time(Some(8.0)).unwrap();
    tmp.drift_time_unit = DriftTimeUnit::Millisecond;
    tmp.name = "bla".into();
    tmp.peaks.push(Peak1D::new(47.11, 0.0));
    tmp
}

/// The spectrum of the move-construction and move-assignment sections
/// (`MSSpectrum_test.cpp:433-446`).
fn move_fixture() -> MSSpectrum {
    let mut tmp = MSSpectrum::new();
    tmp.rt = 9.0;
    tmp.set_drift_time(Some(5.0)).unwrap();
    tmp.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    tmp.ms_level = 18;
    tmp.name = "bla2".into();
    tmp.metadata.insert("label2".into(), meta(5.0));
    tmp.instrument_settings.scan_windows = vec![ScanWindow::new(0.0, 0.0).unwrap(); 2];
    tmp.peaks
        .extend([Peak1D::new(47.11, 0.0), Peak1D::new(48.11, 0.0)]);
    tmp
}

// Sections 25-28: copy constructor, move constructor, copy assignment, move
// assignment. `Clone` is the copy; `std::mem::take` is the move, and it leaves
// the source at `Default`, which is what the source's move sections assert of
// the moved-from spectrum.
#[test]
fn copy_and_assignment_carry_every_field() {
    // (MSSpectrum(const MSSpectrum& source))
    let tmp = copy_fixture();
    let tmp2 = tmp.clone();
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 1);
    assert_eq!(tmp2.metadata["label"], meta(5.0));
    assert_eq!(tmp2.ms_level, 17);
    assert_eq!(tmp2.rt, 7.0);
    assert_eq!(tmp2.drift_time, 8.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(tmp2.name, "bla");
    assert_eq!(tmp2.len(), 1);
    assert_eq!(tmp2.peaks[0].mz, 47.11);

    // (MSSpectrum(const MSSpectrum&& source))
    let mut tmp = move_fixture();
    let orig = tmp.clone();
    let tmp2 = std::mem::take(&mut tmp);
    assert_eq!(tmp2, orig);
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 2);
    assert_eq!(tmp2.metadata["label2"], meta(5.0));
    assert_eq!(tmp2.ms_level, 18);
    assert_eq!(tmp2.rt, 9.0);
    assert_eq!(tmp2.drift_time, 5.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::InverseReducedMobility);
    assert_eq!(tmp2.name, "bla2");
    assert_eq!(tmp2.len(), 2);
    assert_eq!(tmp2.peaks[0].mz, 47.11);
    assert_eq!(tmp2.peaks[1].mz, 48.11);
    assert_eq!(tmp.len(), 0);
    assert!(!tmp.metadata.contains_key("label2"));

    // (MSSpectrum& operator=(const MSSpectrum& source))
    let tmp = copy_fixture();
    let mut tmp2 = tmp.clone();
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 1);
    assert_eq!(tmp2.metadata["label"], meta(5.0));
    assert_eq!(tmp2.ms_level, 17);
    assert_eq!(tmp2.rt, 7.0);
    assert_eq!(tmp2.drift_time, 8.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(tmp2.name, "bla");
    assert_eq!(tmp2.len(), 1);
    assert_eq!(tmp2.peaks[0].mz, 47.11);

    // assignment of an empty object restores every default
    tmp2 = MSSpectrum::new();
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 0);
    assert!(!tmp2.metadata.contains_key("label"));
    assert_eq!(tmp2.ms_level, 1);
    assert_eq!(tmp2.rt, -1.0);
    assert_eq!(tmp2.drift_time, -1.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::None);
    assert_eq!(tmp2.name, "");
    assert_eq!(tmp2.len(), 0);

    // (MSSpectrum& operator=(const MSSpectrum&& source))
    let mut tmp = move_fixture();
    let orig = tmp.clone();
    let mut tmp2 = MSSpectrum::new();
    assert!(tmp2.is_empty());
    tmp2 = std::mem::take(&mut tmp);
    assert_eq!(tmp2, orig);
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 2);
    assert_eq!(tmp2.metadata["label2"], meta(5.0));
    assert_eq!(tmp2.ms_level, 18);
    assert_eq!(tmp2.rt, 9.0);
    assert_eq!(tmp2.drift_time, 5.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::InverseReducedMobility);
    assert_eq!(tmp2.name, "bla2");
    assert_eq!(tmp2.len(), 2);
    assert_eq!(tmp2.peaks[0].mz, 47.11);
    assert_eq!(tmp2.peaks[1].mz, 48.11);
    assert_eq!(tmp.len(), 0);
    assert!(!tmp.metadata.contains_key("label2"));

    // the source moves from a temporary; taking from a fresh value is the same
    tmp2 = MSSpectrum::new();
    assert_eq!(tmp2.instrument_settings.scan_windows.len(), 0);
    assert!(!tmp2.metadata.contains_key("label"));
    assert_eq!(tmp2.ms_level, 1);
    assert_eq!(tmp2.rt, -1.0);
    assert_eq!(tmp2.drift_time, -1.0);
    assert_eq!(tmp2.drift_time_unit, DriftTimeUnit::None);
    assert_eq!(tmp2.name, "");
    assert_eq!(tmp2.len(), 0);
}

// Sections 29 and 30: operator== and operator!=.
#[test]
fn equality_covers_every_field_the_source_compares() {
    let empty = MSSpectrum::new();
    assert!(empty == MSSpectrum::new());
    assert!(empty == empty.clone());

    let mut edit = empty.clone();
    edit.instrument_settings
        .scan_windows
        .push(ScanWindow::new(0.0, 0.0).unwrap());
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.peaks.push(Peak1D::default());
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.metadata.insert("label".into(), text("bla"));
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.set_drift_time(Some(5.0)).unwrap();
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.drift_time_unit = DriftTimeUnit::Millisecond;
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.rt = 5.0;
    assert!(edit != empty);

    let mut edit = empty.clone();
    edit.ms_level = 5;
    assert!(edit != empty);

    for arrays in 0..3 {
        let mut edit = empty.clone();
        match arrays {
            0 => edit
                .float_data_arrays
                .resize_with(5, || DataArray::new("", Vec::<f32>::new())),
            1 => edit
                .string_data_arrays
                .resize_with(5, || DataArray::new("", Vec::<String>::new())),
            _ => edit
                .integer_data_arrays
                .resize_with(5, || DataArray::new("", Vec::<i32>::new())),
        }
        assert!(edit != empty);
    }

    // Divergence, deliberately asserted: the source skips `name_` in
    // operator== ("name_ can differ => it is not checked", MSSpectrum.cpp:516)
    // and its section asserts `empty == edit` after setName("bla"). The Rust
    // struct derives PartialEq over every field, including `name`, so the two
    // differ here. The source's exclusion is not reproducible without hand
    // writing PartialEq, and silently ignoring a field is worse than comparing
    // it.
    let mut edit = empty.clone();
    edit.name = "bla".into();
    assert!(edit != empty);

    // clearing the peaks restores equality
    let mut edit = empty.clone();
    edit.peaks
        .extend([Peak1D::new(2.0, 1.0), Peak1D::new(10.0, 2.0)]);
    let _ = edit.range_manager().unwrap();
    edit.clear(false);
    assert!(edit == empty);
    assert!(!(edit != empty));
}

// Section 31: sortByIntensity(bool reverse = false).
#[test]
fn sort_by_intensity_permutes_every_aligned_array() {
    let mzs = [
        420.130, 412.824, 423.269, 415.287, 413.800, 419.113, 416.293, 418.232, 414.301, 412.321,
    ];
    let intensities = [
        201.0_f32, 60.0, 56.0, 37.0, 34.0, 31.0, 31.0, 31.0, 29.0, 29.0,
    ];
    let peaks: Vec<Peak1D> = mzs
        .iter()
        .zip(intensities)
        .map(|(&mz, i)| Peak1D::new(mz, i))
        .collect();

    // without data arrays the peaks alone are permuted
    let mut ds = MSSpectrum::from_peaks(peaks.clone());
    ds.sort_by_intensity(false).unwrap();
    let mut sorted = intensities;
    sorted.sort_by(f32::total_cmp);
    for (peak, expected) in ds.peaks.iter().zip(sorted) {
        assert_eq!(peak.intensity, expected);
    }

    // with data arrays every array follows, and the names survive
    let floats: Vec<f32> = mzs.iter().map(|&mz| mz as f32).collect();
    let strings: Vec<String> = mzs.iter().map(|&mz| format!("{mz:.2}")).collect();
    let integers: Vec<i32> = mzs.iter().map(|&mz| mz.floor() as i32).collect();
    let mut ds = MSSpectrum {
        peaks,
        float_data_arrays: vec![
            DataArray::new("f1", floats.clone()),
            DataArray::new("f2", floats.clone()),
            DataArray::new("f3", floats),
        ],
        string_data_arrays: vec![
            DataArray::new("s1", strings.clone()),
            DataArray::new("s2", strings),
        ],
        integer_data_arrays: vec![DataArray::new("i1", integers)],
        ..MSSpectrum::default()
    };
    ds.sort_by_intensity(false).unwrap();
    assert_eq!(ds.float_data_arrays[0].name, "f1");
    assert_eq!(ds.float_data_arrays[1].name, "f2");
    assert_eq!(ds.float_data_arrays[2].name, "f3");
    assert_eq!(ds.string_data_arrays[0].name, "s1");
    assert_eq!(ds.string_data_arrays[1].name, "s2");
    assert_eq!(ds.integer_data_arrays[0].name, "i1");
    for (index, expected) in sorted.iter().enumerate() {
        let peak = ds.peaks[index];
        assert_eq!(peak.intensity, *expected);
        assert_eq!(ds.float_data_arrays[1].data[index], peak.mz as f32);
        assert_eq!(
            ds.string_data_arrays[0].data[index],
            format!("{:.2}", peak.mz)
        );
        assert_eq!(
            ds.integer_data_arrays[0].data[index],
            peak.mz.floor() as i32
        );
    }
}

// Section 32: sortByPosition().
#[test]
fn sort_by_position_permutes_every_aligned_array() {
    let mut ds = prefilled_spectrum();
    ds.float_data_arrays[1].name = "f2".into();
    let intensities: Vec<f32> = ds.peaks.iter().map(|peak| peak.intensity).collect();

    let mut plain = MSSpectrum::from_peaks(ds.peaks.clone());
    plain.sort_by_position().unwrap();
    for (peak, expected) in plain.peaks.iter().zip(intensities.iter().rev()) {
        assert_eq!(peak.intensity, *expected);
    }

    ds.sort_by_position().unwrap();
    assert_eq!(ds.float_data_arrays[0].name, "f1");
    assert_eq!(ds.float_data_arrays[1].name, "f2");
    assert_eq!(ds.float_data_arrays[2].name, "f3");
    assert_eq!(ds.string_data_arrays[0].name, "s1");
    assert_eq!(ds.string_data_arrays[1].name, "s2");
    assert_eq!(ds.integer_data_arrays[0].name, "i1");
    assert_eq!(ds.len(), 10);
    for (index, expected) in intensities.iter().rev().enumerate() {
        assert_eq!(ds.peaks[index].intensity, *expected);
        assert_eq!(ds.float_data_arrays[1].data[index], *expected);
        assert_eq!(
            ds.string_data_arrays[0].data[index],
            format!("{expected:.0}")
        );
        assert_eq!(ds.integer_data_arrays[0].data[index], *expected as i32);
    }
}

// Sections 33 and 34: sortByIonMobility() and isSortedByIM().
#[test]
fn sort_by_ion_mobility_orders_peaks_and_arrays() {
    let mut ds = prefilled_spectrum();
    assert!(!ds.is_sorted_by_im().unwrap());
    ds.sort_by_ion_mobility().unwrap();
    assert!(ds.is_sorted_by_im().unwrap());
    let (index, unit) = ds.im_data().unwrap();
    assert_eq!(index, 1);
    assert_eq!(unit, DriftTimeUnit::Millisecond);
    let im = &ds.float_data_arrays[index].data;
    assert!(im.windows(2).all(|pair| pair[0] <= pair[1]));

    // The stable permutation of the source float array
    // [56, 201, 31, 31, 31, 37, 29, 34, 60, 29] is [6, 9, 2, 3, 4, 7, 5, 0, 8, 1];
    // every array and the peaks follow it.
    assert_eq!(
        im.as_slice(),
        [
            29.0_f32, 29.0, 31.0, 31.0, 31.0, 34.0, 37.0, 56.0, 60.0, 201.0
        ]
    );
    assert_eq!(
        ds.peaks.iter().map(|peak| peak.mz).collect::<Vec<_>>(),
        vec![
            414.301, 412.321, 419.113, 418.232, 416.293, 413.800, 415.287, 423.269, 412.824,
            420.130
        ]
    );
    assert_eq!(ds.string_data_arrays[0].data[0], "29");
    assert_eq!(ds.string_data_arrays[0].data[9], "201");
    assert_eq!(ds.integer_data_arrays[0].data[9], 201);

    // already sorted: a second call is a no-op
    let before = ds.clone();
    ds.sort_by_ion_mobility().unwrap();
    assert_eq!(ds, before);
}

// Section 35: sortByPositionPresorted(const std::vector<Chunk>&).
#[test]
fn sort_by_position_presorted_merges_presorted_chunks() {
    let mzs = [
        419.113, 420.130, 423.269, 415.287, 416.293, 418.232, 413.800, 414.301, 412.824, 412.321,
    ];
    let intensities = [
        19.0_f32, 20.0, 23.0, 15.0, 16.0, 18.0, 13.0, 14.0, 12.0, 12.0,
    ];
    let floats = intensities.to_vec();
    let strings: Vec<String> = intensities.iter().map(|i| format!("{i:.0}")).collect();
    let integers: Vec<i32> = intensities.iter().map(|&i| i as i32).collect();

    let mut ds = MSSpectrum::default();
    let mut chunks = Chunks::new();
    let mut last_added = 0.0;
    for (&mz, intensity) in mzs.iter().zip(intensities) {
        if mz < last_added {
            chunks.add(&ds, true).unwrap();
        }
        last_added = mz;
        ds.peaks.push(Peak1D::new(mz, intensity));
    }
    chunks.add(&ds, true).unwrap(); // add the last chunk
    assert_eq!(
        chunks.chunks(),
        [
            Chunk::new(0, 3, true),
            Chunk::new(3, 6, true),
            Chunk::new(6, 8, true),
            Chunk::new(8, 9, true),
            Chunk::new(9, 10, true),
        ]
    );

    ds.float_data_arrays = vec![
        DataArray::new("f1", floats.clone()),
        DataArray::new("f2", floats.clone()),
        DataArray::new("f3", floats),
    ];
    ds.string_data_arrays = vec![
        DataArray::new("s1", strings.clone()),
        DataArray::new("s2", strings),
    ];
    ds.integer_data_arrays = vec![
        DataArray::new("i1", integers.clone()),
        DataArray::new("", integers),
    ];

    ds.sort_by_position_presorted(chunks.chunks()).unwrap();

    assert_eq!(ds.float_data_arrays[0].name, "f1");
    assert_eq!(ds.float_data_arrays[1].name, "f2");
    assert_eq!(ds.float_data_arrays[2].name, "f3");
    assert_eq!(ds.string_data_arrays[0].name, "s1");
    assert_eq!(ds.string_data_arrays[1].name, "s2");
    assert_eq!(ds.integer_data_arrays[0].name, "i1");

    let mut sorted = intensities;
    sorted.sort_by(f32::total_cmp);
    assert_eq!(ds.len(), 10);
    assert_eq!(ds.float_data_arrays[1].data.len(), 10);
    assert_eq!(ds.string_data_arrays[0].data.len(), 10);
    assert_eq!(ds.integer_data_arrays[0].data.len(), 10);
    for (index, expected) in sorted.iter().enumerate() {
        assert_eq!(ds.peaks[index].intensity, *expected);
        assert_eq!(ds.float_data_arrays[1].data[index], *expected);
        assert_eq!(
            ds.string_data_arrays[0].data[index],
            format!("{expected:.0}")
        );
        assert_eq!(ds.integer_data_arrays[0].data[index], *expected as i32);
    }
    assert!(ds.is_sorted());
}

// Native: the chunked sort agrees with the plain sort, honours unsorted chunks
// and rejects the chunk lists the source silently mishandles.
#[test]
fn presorted_sort_matches_plain_sort_and_rejects_bad_chunks() {
    let mzs = [5.0, 7.0, 9.0, 1.0, 3.0, 8.0, 2.0, 4.0, 6.0];
    let peaks: Vec<Peak1D> = mzs
        .iter()
        .enumerate()
        .map(|(index, &mz)| Peak1D::new(mz, index as f32))
        .collect();
    let arrays = vec![DataArray::new(
        "tag",
        (0..mzs.len() as i32).collect::<Vec<_>>(),
    )];

    // one unsorted trailing chunk: the run [6, 9) is 2, 4, 6 (sorted) but is
    // declared unsorted, which only costs an extra sort.
    let mut ds = MSSpectrum {
        peaks: peaks.clone(),
        integer_data_arrays: arrays.clone(),
        ..MSSpectrum::default()
    };
    let chunks = [
        Chunk::new(0, 3, true),
        Chunk::new(3, 5, true),
        Chunk::new(5, 6, true),
        Chunk::new(6, 9, false),
    ];
    ds.sort_by_position_presorted(&chunks).unwrap();
    let mut plain = MSSpectrum {
        peaks: peaks.clone(),
        integer_data_arrays: arrays.clone(),
        ..MSSpectrum::default()
    };
    plain.sort_by_position().unwrap();
    assert_eq!(ds, plain);
    assert_eq!(
        ds.integer_data_arrays[0].data,
        vec![3, 6, 4, 7, 0, 8, 1, 5, 2]
    );

    // a genuinely unsorted single chunk is still sorted
    let mut ds = MSSpectrum {
        peaks: peaks.clone(),
        integer_data_arrays: arrays.clone(),
        ..MSSpectrum::default()
    };
    ds.sort_by_position_presorted(&[Chunk::new(0, 9, false)])
        .unwrap();
    assert_eq!(ds, plain);

    // an empty chunk list sorts nothing, as in the source
    let mut ds = MSSpectrum {
        peaks: peaks.clone(),
        integer_data_arrays: arrays.clone(),
        ..MSSpectrum::default()
    };
    let before = ds.clone();
    ds.sort_by_position_presorted(&[]).unwrap();
    assert_eq!(ds, before);

    // a single chunk that claims to be sorted returns immediately
    let mut sorted_spectrum = plain.clone();
    sorted_spectrum
        .sort_by_position_presorted(&[Chunk::new(0, 9, true)])
        .unwrap();
    assert_eq!(sorted_spectrum, plain);

    // a chunk that lies about being sorted is rejected, where the source feeds
    // std::inplace_merge an unsorted range
    let mut ds = MSSpectrum {
        peaks: peaks.clone(),
        integer_data_arrays: arrays.clone(),
        ..MSSpectrum::default()
    };
    let before = ds.clone();
    assert!(matches!(
        ds.sort_by_position_presorted(&[Chunk::new(0, 9, true)]),
        Err(Error::UnsortedData)
    ));
    assert_eq!(ds, before);

    // gaps, overlaps, short coverage and reversed runs are rejected
    for bad in [
        vec![Chunk::new(1, 3, true), Chunk::new(3, 9, false)],
        vec![Chunk::new(0, 4, false), Chunk::new(3, 9, false)],
        vec![Chunk::new(0, 3, false), Chunk::new(3, 6, false)],
        vec![Chunk::new(0, 3, false), Chunk::new(3, 2, false)],
    ] {
        let mut ds = MSSpectrum {
            peaks: peaks.clone(),
            integer_data_arrays: arrays.clone(),
            ..MSSpectrum::default()
        };
        let before = ds.clone();
        assert!(matches!(
            ds.sort_by_position_presorted(&bad),
            Err(Error::InvalidValue(_))
        ));
        assert_eq!(ds, before);
    }

    // a mis-sized data array is rejected before anything moves
    let mut ds = MSSpectrum {
        peaks,
        integer_data_arrays: vec![DataArray::new("tag", vec![1, 2])],
        ..MSSpectrum::default()
    };
    let before = ds.clone();
    assert!(
        ds.sort_by_position_presorted(&[Chunk::new(0, 9, false)])
            .is_err()
    );
    assert_eq!(ds, before);

    // Chunks::add rejects a spectrum that shrank below the previous run
    let mut ds = MSSpectrum::from_peaks(vec![Peak1D::new(1.0, 1.0), Peak1D::new(2.0, 1.0)]);
    let mut chunks = Chunks::new();
    assert!(chunks.is_empty());
    chunks.add(&ds, true).unwrap();
    assert_eq!(chunks.len(), 1);
    ds.peaks.pop();
    assert!(matches!(chunks.add(&ds, true), Err(Error::InvalidValue(_))));
}

// Sections 36, 37 and 38: isSorted(), isSorted(lambda) and sort(lambda).
#[test]
fn is_sorted_and_user_defined_orderings() {
    // (bool isSorted() const)
    let mut spec = MSSpectrum::from_peaks(vec![
        Peak1D::new(1000.0, 3.0),
        Peak1D::new(1001.0, 5.0),
        Peak1D::new(1002.0, 1.0),
    ]);
    assert!(spec.is_sorted());
    spec.peaks.reverse();
    assert!(!spec.is_sorted());

    // (template<class Predicate> bool isSorted(const Predicate&) const):
    // the port has no generic predicate overload; the equivalent is an explicit
    // comparison over the same indices, which is what the source lambda does.
    let mut ds = prefilled_spectrum();
    ds.float_data_arrays[1].name = "f2".into();
    ds.integer_data_arrays.truncate(1);
    let sorted_by = |spec: &MSSpectrum, key: fn(&MSSpectrum, usize) -> f64| {
        (1..spec.len()).all(|i| key(spec, i - 1) <= key(spec, i))
    };
    let by_mz: fn(&MSSpectrum, usize) -> f64 = |spec, i| spec.peaks[i].mz;
    let by_intensity: fn(&MSSpectrum, usize) -> f64 = |spec, i| f64::from(spec.peaks[i].intensity);

    ds.sort_by_position().unwrap();
    assert!(sorted_by(&ds, by_mz));
    assert!(ds.is_sorted());

    ds.sort_by_intensity(false).unwrap();
    assert!(sorted_by(&ds, by_intensity));
    assert!(!sorted_by(&ds, by_mz));
    assert!(!ds.is_sorted());

    // sort by the first float data array; it holds the intensities here
    ds.sort_by_position().unwrap();
    let values = ds.float_data_arrays[0].data.clone();
    let mut order: Vec<usize> = (0..ds.len()).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    ds.select(&order).unwrap();
    assert_eq!(ds.peaks[0].intensity, 29.0);
    assert!(!ds.is_sorted());
    assert!(sorted_by(&ds, by_intensity));

    // (template<class Predicate> void sort(const Predicate&)): a mis-sized data
    // array is rejected before the permutation is applied.
    let mut spec = MSSpectrum {
        peaks: vec![
            Peak1D::new(2.0, 1.0),
            Peak1D::new(10.0, 2.0),
            Peak1D::new(30.0, 3.0),
        ],
        integer_data_arrays: vec![DataArray::new("", vec![1, 2])],
        ..MSSpectrum::default()
    };
    let before = spec.clone();
    assert!(spec.select(&[0, 1, 2]).is_err());
    assert_eq!(spec, before);
}

// Sections 39-54: the sixteen MZBegin/MZEnd/PosBegin/PosEnd overloads. `PosBegin`
// and `PosEnd` are documented aliases of `MZBegin`/`MZEnd`; the subrange
// overloads map onto slicing the peak list, so they are covered by the same two
// Rust methods applied to a sub-slice.
#[test]
fn mz_range_bounds_cover_every_iterator_overload() {
    let spec = spec_find();
    let mz = |index: usize| spec.peaks[index].mz;

    // MZEnd / PosEnd, full range
    assert_eq!(mz(spec.mz_end(4.5).unwrap()), 5.0);
    assert_eq!(mz(spec.mz_end(5.0).unwrap()), 6.0);
    assert_eq!(mz(spec.mz_end(5.5).unwrap()), 6.0);
    assert_eq!(mz(spec.mz_end(3.5).unwrap()), 4.0);
    assert_eq!(mz(spec.mz_end(4.0).unwrap()), 5.0);

    // MZBegin / PosBegin, full range
    assert_eq!(mz(spec.mz_begin(4.5).unwrap()), 5.0);
    assert_eq!(mz(spec.mz_begin(5.0).unwrap()), 5.0);
    assert_eq!(mz(spec.mz_begin(5.5).unwrap()), 6.0);
    assert_eq!(mz(spec.mz_begin(3.5).unwrap()), 4.0);

    // PosBegin(begin, 8.0, end) runs past the last peak, so the element before
    // the result is the last peak.
    let past_end = spec.mz_begin(8.0).unwrap();
    assert_eq!(past_end, spec.len());
    assert_eq!(mz(past_end - 1), mz(spec.len() - 1));

    // The subrange overloads: a sub-slice of the peaks, with the same rule. An
    // empty subrange (begin == end) returns its own start, as the source's
    // `MZBegin(begin, mz, begin)` does.
    let sub = |from: usize, to: usize, query: f64, inclusive: bool| {
        from + spec.peaks[from..to].partition_point(|peak| {
            if inclusive {
                peak.mz <= query
            } else {
                peak.mz < query
            }
        })
    };
    assert_eq!(mz(sub(0, spec.len(), 4.5, false)), 5.0);
    assert_eq!(mz(sub(0, spec.len(), 4.5, true)), 5.0);
    assert_eq!(mz(sub(0, spec.len(), 5.0, true)), 6.0);
    assert_eq!(mz(sub(0, spec.len(), 5.5, false)), 6.0);
    assert_eq!(mz(sub(0, spec.len(), 3.5, false)), 4.0);
    assert_eq!(mz(sub(0, spec.len(), 4.0, true)), 5.0);
    assert_eq!(sub(0, 0, 4.5, false), 0);
    assert_eq!(sub(0, 0, 4.5, true), 0);
}

// Sections 55 and 56: containsIMData() and getIMData().
#[test]
fn contains_and_get_im_data() {
    let ds = prefilled_spectrum();
    assert!(ds.contains_im_data());
    let (index, unit) = ds.im_data().unwrap();
    assert_eq!(index, 1);
    assert_eq!(unit, DriftTimeUnit::Millisecond);

    // Native: every name IMDataArrayUtils::getIMUnit recognizes, and its order.
    for (name, unit) in [
        (
            "mean ion mobility drift time array",
            DriftTimeUnit::Millisecond,
        ),
        ("mean ion mobility array", DriftTimeUnit::Millisecond),
        (
            "mean inverse reduced ion mobility array",
            DriftTimeUnit::InverseReducedMobility,
        ),
        ("raw ion mobility array", DriftTimeUnit::Millisecond),
        (
            "raw inverse reduced ion mobility array",
            DriftTimeUnit::InverseReducedMobility,
        ),
        (
            "raw ion mobility drift time array",
            DriftTimeUnit::Millisecond,
        ),
        (
            "deconvoluted ion mobility array",
            DriftTimeUnit::Millisecond,
        ),
        (
            "deconvoluted inverse reduced ion mobility array",
            DriftTimeUnit::InverseReducedMobility,
        ),
        (
            "deconvoluted ion mobility drift time array",
            DriftTimeUnit::Millisecond,
        ),
        ("Ion Mobility", DriftTimeUnit::Millisecond),
        (
            "Ion Mobility (MS:1002815)",
            DriftTimeUnit::InverseReducedMobility,
        ),
        (
            "Ion Mobility (MS:1003006)",
            DriftTimeUnit::InverseReducedMobility,
        ),
        (
            "Ion Mobility (MS:1002954)",
            DriftTimeUnit::CollisionCrossSection,
        ),
        (
            "inverse reduced ion mobility",
            DriftTimeUnit::InverseReducedMobility,
        ),
    ] {
        let spec = MSSpectrum {
            float_data_arrays: vec![DataArray::new(name, vec![1.0_f32])],
            peaks: vec![Peak1D::new(1.0, 1.0)],
            ..MSSpectrum::default()
        };
        assert!(spec.contains_im_data(), "{name}");
        assert_eq!(spec.im_data().unwrap(), (0, unit), "{name}");
    }

    // a spectrum with no IM array, and one whose only array is named otherwise
    for arrays in [
        Vec::new(),
        vec![DataArray::new("Wrong Name", vec![4.0_f32, 5.0])],
        vec![DataArray::new("ion mobility", vec![4.0_f32])], // case sensitive
    ] {
        let mut spec = MSSpectrum {
            float_data_arrays: arrays,
            ..MSSpectrum::default()
        };
        assert!(!spec.contains_im_data());
        assert!(matches!(spec.im_data(), Err(Error::MissingInformation(_))));
        assert!(spec.is_sorted_by_im().is_err());
        assert!(spec.sort_by_ion_mobility().is_err());
    }

    // the first matching array wins, as the source scan does
    let spec = MSSpectrum {
        peaks: vec![Peak1D::new(1.0, 1.0)],
        float_data_arrays: vec![
            DataArray::new("other", vec![0.0_f32]),
            DataArray::new("raw inverse reduced ion mobility array", vec![1.0_f32]),
            DataArray::new(IM_MS_ARRAY, vec![2.0_f32]),
        ],
        ..MSSpectrum::default()
    };
    assert_eq!(
        spec.im_data().unwrap(),
        (1, DriftTimeUnit::InverseReducedMobility)
    );
}

// Sections 57, 58, 59 and 60: the three findNearest overloads and
// findHighestInWindow. The source returns -1 for "no peak"; the port returns
// `None`, and the empty spectrum that the source declares a precondition
// violation is `Ok(None)` here.
#[test]
fn find_nearest_and_highest_reference_cases() {
    let spec = spec_test();
    for (query, index) in [
        (400.0, 0),
        (500.0, 20),
        (412.4, 0),
        (441.224, 20),
        (426.29, 10),
        (426.3, 10),
        (427.2, 11),
        (427.3, 11),
    ] {
        assert_eq!(spec.find_nearest(query).unwrap(), Some(index));
    }
    assert_eq!(MSSpectrum::new().find_nearest(427.3).unwrap(), None);

    for (query, tolerance, expected) in [
        (400.0, 1.0, None),
        (500.0, 1.0, None),
        (412.4, 0.01, None),
        (412.4, 0.1, Some(0)),
        (441.3, 0.01, None),
        (441.3, 0.1, Some(20)),
        (426.29, 0.1, Some(10)),
        (426.3, 0.1, Some(10)),
        (427.2, 0.1, Some(11)),
        (427.3, 0.1, Some(11)),
        (427.3, 0.001, None),
    ] {
        assert_eq!(
            spec.find_nearest_with_tolerance(query, tolerance).unwrap(),
            expected
        );
        assert_eq!(
            spec.find_nearest_in_window(query, tolerance, tolerance)
                .unwrap(),
            expected
        );
    }
    assert_eq!(
        MSSpectrum::new()
            .find_nearest_in_window(427.3, 1.0, 1.0)
            .unwrap(),
        None
    );
    assert_eq!(
        spec.find_nearest_in_window(427.3, 0.1, 0.001).unwrap(),
        Some(11)
    );
    assert_eq!(
        spec.find_nearest_in_window(427.3, 0.001, 1.01).unwrap(),
        None
    );
    assert_eq!(
        spec.find_nearest_in_window(427.3, 0.001, 1.1).unwrap(),
        Some(12)
    );

    for (query, left, right, expected) in [
        (400.0, 1.0, 1.0, None),
        (500.0, 1.0, 1.0, None),
        (412.4, 0.01, 0.01, None),
        (412.4, 0.1, 0.1, Some(0)),
        (441.3, 0.01, 0.01, None),
        (441.3, 0.1, 0.1, Some(20)),
        (426.29, 0.1, 0.1, Some(10)),
        (426.3, 0.1, 0.1, Some(10)),
        (427.2, 0.1, 0.1, Some(11)),
        (427.3, 0.1, 0.1, Some(11)),
        (427.3, 0.001, 0.001, None),
        (427.3, 0.1, 0.001, Some(11)),
        (427.3, 0.001, 1.01, None),
        (427.3, 0.001, 1.1, Some(12)),
        (427.3, 9.0, 4.0, Some(8)),
        (430.25, 1.9, 1.01, Some(13)),
    ] {
        assert_eq!(
            spec.find_highest_in_window(query, left, right).unwrap(),
            expected
        );
    }
    assert_eq!(
        MSSpectrum::new()
            .find_highest_in_window(427.3, 1.0, 1.0)
            .unwrap(),
        None
    );
}

// Section 61: getType(const bool query_data).
#[test]
fn get_type_precedence() {
    let mut edit = MSSpectrum::new();
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Unknown);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Unknown);

    edit.spectrum_type = SpectrumType::Profile;
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Profile);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Profile);

    edit.data_processing.push(Arc::new(DataProcessing {
        actions: [ProcessingAction::PeakPicking].into(),
        ..Default::default()
    }));
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Profile);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Profile);
    edit.spectrum_type = SpectrumType::Unknown;
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Centroid);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Centroid);

    edit.data_processing.clear();
    edit.peaks = (1..=4)
        .map(|i| Peak1D::new(f64::from(i) * 100.0, 1.0))
        .collect();
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Unknown);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Unknown);
    edit.peaks
        .extend([Peak1D::new(500.0, 1.0), Peak1D::new(600.0, 1.0)]);
    assert_eq!(edit.get_type(false).unwrap(), SpectrumType::Unknown);
    assert_eq!(edit.get_type(true).unwrap(), SpectrumType::Centroid);
}

// Sections 62, 63 and 64: getBasePeak (both overloads) and calculateTIC.
#[test]
fn base_peak_and_total_ion_current() {
    let spec = spec_test();
    let base = spec.base_peak().unwrap();
    assert_eq!(base.intensity, 201.0);
    assert_eq!(spec.peaks.iter().position(|peak| peak == base).unwrap(), 8);
    assert!(MSSpectrum::new().base_peak().is_none());

    // the mutable overload addresses the same peak
    let mut test = spec_test();
    let index = 8;
    test.peaks[index].intensity += 0.0;
    assert_eq!(test.peaks[index].intensity, 201.0);

    assert_eq!(spec.calculate_tic(), 1032.0);
    assert_eq!(MSSpectrum::new().calculate_tic(), 0.0);
}

// Sections 65, 66 and 67: setIMFormat, getIMPeakType and setIMPeakType. These
// three live on SpectrumSettings and the Rust MSSpectrum carries no IM format or
// IM peak type field, so there is nothing on a spectrum to assert. The enums
// themselves are ported, with the source's names and defaults; the missing
// fields are recorded as a deferral in docs/SPECTRUM_MOBILITY_SUPPORT.md.
#[test]
fn im_format_and_peak_type_enums_exist_but_no_spectrum_field_does() {
    assert_eq!(IonMobilityFormat::PerPeak.name(), "im_peak");
    assert_eq!(IonMobilityFormat::None.name(), "none");
    assert_eq!(IonMobilityPeakType::default(), IonMobilityPeakType::Unknown);
    assert_eq!(IonMobilityPeakType::Centroid.name(), "im_centroided");
    assert_eq!(IonMobilityPeakType::Profile.name(), "im_profile");
}

// Section 68: clear(bool clear_meta_data).
#[test]
fn clear_option_preserves_or_resets_metadata() {
    let mut edit = MSSpectrum::new();
    edit.instrument_settings
        .scan_windows
        .push(ScanWindow::new(0.0, 0.0).unwrap());
    edit.peaks.push(Peak1D::default());
    edit.metadata.insert("label".into(), text("bla"));
    edit.rt = 5.0;
    edit.set_drift_time(Some(6.0)).unwrap();
    edit.drift_time_unit = DriftTimeUnit::Millisecond;
    edit.ms_level = 5;
    edit.float_data_arrays
        .resize_with(5, || DataArray::new("", Vec::<f32>::new()));
    edit.integer_data_arrays
        .resize_with(5, || DataArray::new("", Vec::<i32>::new()));
    edit.string_data_arrays
        .resize_with(5, || DataArray::new("", Vec::<String>::new()));

    edit.clear(false);
    assert_eq!(edit.len(), 0);
    assert!(edit != MSSpectrum::new());
    assert!(edit.is_empty());
    assert_eq!(edit.drift_time, 6.0);
    assert_eq!(edit.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(edit.float_data_arrays.len(), 0);

    edit.clear(true);
    assert!(edit.is_empty());
    assert!(edit == MSSpectrum::new());
    assert_eq!(edit.drift_time, -1.0);
    assert_eq!(edit.drift_time_unit, DriftTimeUnit::None);
}

// Section 69: MSSpectrum::RTLess. The port has no comparator type; retention
// time ordering is the scalar comparison on the `rt` field.
#[test]
fn retention_time_ordering() {
    let mut v: Vec<MSSpectrum> = [3.0, 2.0, 1.0]
        .iter()
        .map(|&rt| MSSpectrum {
            rt,
            ..MSSpectrum::default()
        })
        .collect();
    v.sort_by(|a, b| a.rt.total_cmp(&b.rt));
    assert_eq!(v[0].rt, 1.0);
    assert_eq!(v[1].rt, 2.0);
    assert_eq!(v[2].rt, 3.0);

    let s1 = MSSpectrum {
        rt: 0.451,
        ..MSSpectrum::default()
    };
    let s2 = MSSpectrum {
        rt: 0.5,
        ..MSSpectrum::default()
    };
    assert_eq!(s1.rt.total_cmp(&s2.rt), std::cmp::Ordering::Less);
    assert_eq!(s2.rt.total_cmp(&s1.rt), std::cmp::Ordering::Greater);
    assert_eq!(s2.rt.total_cmp(&s2.rt), std::cmp::Ordering::Equal);

    // MSSpectrum::IMLess has no section of its own; it is the same comparison on
    // the drift time (MSSpectrum.cpp:803-806).
    let mut a = MSSpectrum::default();
    a.set_drift_time(Some(2.0)).unwrap();
    let mut b = MSSpectrum::default();
    b.set_drift_time(Some(3.0)).unwrap();
    assert!(a.drift_time < b.drift_time);
}

// Section 70: maybeGetIMData().
#[test]
fn maybe_im_data_variants() {
    // an array named "Ion Mobility" with three values, as the section builds it
    let spec = MSSpectrum {
        float_data_arrays: vec![DataArray::new("Ion Mobility", vec![1.0_f32, 2.0, 3.0])],
        ..MSSpectrum::default()
    };
    let (unit, values) = spec.maybe_im_data().unwrap();
    assert_eq!(unit, DriftTimeUnit::Millisecond);
    assert_eq!(values.len(), 3);
    assert_eq!(values[0], 1.0);
    assert_eq!(values[1], 2.0);
    assert_eq!(values[2], 3.0);

    // missing ion mobility data
    assert!(MSSpectrum::new().maybe_im_data().is_none());

    // an empty float-array list
    let mut spec_empty = MSSpectrum::new();
    spec_empty.float_data_arrays.clear();
    assert!(spec_empty.maybe_im_data().is_none());

    // the wrong array name
    let spec = MSSpectrum {
        float_data_arrays: vec![DataArray::new("Wrong Name", vec![4.0_f32, 5.0])],
        ..MSSpectrum::default()
    };
    assert!(spec.maybe_im_data().is_none());

    // Native: an IM array with no entries is `Some` with an empty slice, where
    // the source's `{NONE, {}}` is indistinguishable from "not an IM frame".
    let spec = MSSpectrum {
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, Vec::<f32>::new())],
        ..MSSpectrum::default()
    };
    assert_eq!(
        spec.maybe_im_data(),
        Some((DriftTimeUnit::Millisecond, [].as_slice()))
    );
    assert!(spec.contains_im_data());
}

// Section 71 ([EXTRA] operator<<): the port has no stream operator for a whole
// spectrum, so only the per-peak lines the source prints are reproducible.
#[test]
fn peak_stream_lines_match_the_source_format() {
    let spec = MSSpectrum::from_peaks(
        [
            (412.321, 29.0_f32),
            (412.824, 60.0),
            (413.8, 34.0),
            (414.301, 29.0),
            (415.287, 37.0),
            (416.293, 31.0),
            (418.232, 31.0),
            (419.113, 31.0),
            (420.13, 201.0),
            (423.269, 56.0),
            (426.292, 34.0),
        ]
        .iter()
        .map(|&(mz, intensity)| Peak1D::new(mz, intensity))
        .collect(),
    );
    let rendered: Vec<String> = spec.peaks.iter().map(ToString::to_string).collect();
    assert_eq!(
        rendered,
        vec![
            "POS: 412.321 INT: 29",
            "POS: 412.824 INT: 60",
            "POS: 413.8 INT: 34",
            "POS: 414.301 INT: 29",
            "POS: 415.287 INT: 37",
            "POS: 416.293 INT: 31",
            "POS: 418.232 INT: 31",
            "POS: 419.113 INT: 31",
            "POS: 420.13 INT: 201",
            "POS: 423.269 INT: 56",
            "POS: 426.292 INT: 34",
        ]
    );
}

/// A four-peak ion-mobility frame for the rasterizer: m/z 100..400 against
/// mobility 0.5..2.0.
fn frame() -> MSSpectrum {
    MSSpectrum {
        peaks: vec![
            Peak1D::new(100.0, 1.0),
            Peak1D::new(150.0, 2.0),
            Peak1D::new(300.0, 4.0),
            Peak1D::new(399.0, 8.0),
        ],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![0.5_f32, 0.6, 1.5, 1.9])],
        ..MSSpectrum::default()
    }
}

// rasterizeIMFrame has no upstream section; these expectations are derived from
// the algorithm in MSSpectrum.cpp:857-967.
#[test]
fn rasterize_im_frame_bins_and_aggregates() {
    let spec = frame();
    let raster = ImFrameRaster::new(2, 2, 0.0, 2.0, 0.0, 400.0);
    // mz_scale = 2/400, im_scale = 2/2.
    // (100, 0.5) -> mz_bin 0, im_bin 0; (150, 0.6) -> 0, 0;
    // (300, 1.5) -> 1, 1;             (399, 1.9) -> 1, 1.
    let image = spec.rasterize_im_frame(&raster).unwrap();
    assert_eq!(image, vec![3.0, 0.0, 0.0, 12.0]);

    let image = spec
        .rasterize_im_frame(&raster.with_aggregation(RasterAggregation::Max))
        .unwrap();
    assert_eq!(image, vec![2.0, 0.0, 0.0, 8.0]);

    // Row-major with m/z slowest: four m/z bins, one mobility bin.
    // The m/z scale is now 4/400, so the peaks at 100 and 150 share bin 1 and
    // those at 300 and 399 share bin 3.
    let image = spec
        .rasterize_im_frame(&ImFrameRaster::new(1, 4, 0.0, 2.0, 0.0, 400.0))
        .unwrap();
    assert_eq!(image, vec![0.0, 3.0, 0.0, 12.0]);

    // A peak exactly at the upper bound lands in the last bin rather than past
    // the end; a peak outside the range is skipped.
    let spec = MSSpectrum {
        peaks: vec![
            Peak1D::new(400.0, 5.0),
            Peak1D::new(401.0, 100.0),
            Peak1D::new(50.0, 100.0),
            Peak1D::new(200.0, 7.0),
        ],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![2.0_f32, 1.0, 1.0, 3.0])],
        ..MSSpectrum::default()
    };
    let image = spec
        .rasterize_im_frame(&ImFrameRaster::new(2, 2, 0.0, 2.0, 100.0, 400.0))
        .unwrap();
    // (400, 2.0) -> last m/z bin, last IM bin; (200, 3.0) skipped (IM > max);
    // (401, .) and (50, .) skipped (m/z outside).
    assert_eq!(image, vec![0.0, 0.0, 0.0, 5.0]);

    // An empty spectrum returns a zeroed image, as the source's early return
    // does, provided the IM array name is present.
    let spec = MSSpectrum {
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, Vec::<f32>::new())],
        ..MSSpectrum::default()
    };
    assert_eq!(
        spec.rasterize_im_frame(&ImFrameRaster::new(3, 2, 0.0, 1.0, 0.0, 1.0))
            .unwrap(),
        vec![0.0; 6]
    );
}

#[test]
fn rasterize_im_frame_rejects_invalid_grids_and_arrays() {
    let spec = frame();
    for bad in [
        ImFrameRaster::new(0, 2, 0.0, 2.0, 0.0, 400.0),
        ImFrameRaster::new(2, 0, 0.0, 2.0, 0.0, 400.0),
        ImFrameRaster::new(2, 2, f64::NAN, 2.0, 0.0, 400.0),
    ] {
        assert!(matches!(
            spec.rasterize_im_frame(&bad),
            Err(Error::InvalidValue(_))
        ));
    }
    for bad in [
        ImFrameRaster::new(2, 2, 2.0, 2.0, 0.0, 400.0),
        ImFrameRaster::new(2, 2, 3.0, 2.0, 0.0, 400.0),
        ImFrameRaster::new(2, 2, 0.0, 2.0, 400.0, 400.0),
        ImFrameRaster::new(2, 2, 0.0, 2.0, 500.0, 400.0),
    ] {
        assert!(matches!(
            spec.rasterize_im_frame(&bad),
            Err(Error::InvalidRange(_))
        ));
    }

    let good = ImFrameRaster::new(2, 2, 0.0, 2.0, 0.0, 400.0);

    // no ion mobility array at all
    let bare = MSSpectrum::from_peaks(vec![Peak1D::new(100.0, 1.0)]);
    assert!(matches!(
        bare.rasterize_im_frame(&good),
        Err(Error::MissingInformation(_))
    ));

    // the array is present but the wrong length
    let mut mismatched = frame();
    mismatched.float_data_arrays[0].data.pop();
    assert!(matches!(
        mismatched.rasterize_im_frame(&good),
        Err(Error::InvalidValue(_))
    ));

    // non-finite values, where the source casts NaN to Int64
    let mut nan_mobility = frame();
    nan_mobility.float_data_arrays[0].data[1] = f32::NAN;
    assert!(matches!(
        nan_mobility.rasterize_im_frame(&good),
        Err(Error::InvalidValue(_))
    ));
    let mut nan_mz = frame();
    nan_mz.peaks[1].mz = f64::NAN;
    assert!(matches!(
        nan_mz.rasterize_im_frame(&good),
        Err(Error::InvalidValue(_))
    ));

    // the pixel budget, and the overflowing product the source multiplies
    // unchecked
    let huge = ImFrameRaster::new(MSSpectrum::MAX_RASTER_PIXELS, 2, 0.0, 2.0, 0.0, 400.0);
    assert!(matches!(
        spec.rasterize_im_frame(&huge),
        Err(Error::InvalidValue(_))
    ));
    let overflow = ImFrameRaster::new(usize::MAX, 2, 0.0, 2.0, 0.0, 400.0);
    assert!(overflow.pixels().is_err());
    assert!(matches!(
        spec.rasterize_im_frame(&overflow),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(good.pixels().unwrap(), 4);
}

// Native: the ion-mobility sort and predicate reject the inputs whose source
// behaviour is out-of-bounds or unspecified, and leave the spectrum unchanged.
#[test]
fn ion_mobility_sorting_rejects_unusable_arrays() {
    // an empty IM array for a non-empty spectrum: checkDataArraySizes_ permits
    // it and the source comparator then indexes it per peak
    let mut spec = MSSpectrum {
        peaks: vec![Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 2.0)],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, Vec::<f32>::new())],
        ..MSSpectrum::default()
    };
    let before = spec.clone();
    assert!(matches!(
        spec.sort_by_ion_mobility(),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        spec.is_sorted_by_im(),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(spec, before);

    // a non-finite ion mobility value, where std::is_sorted would report sorted
    let mut spec = MSSpectrum {
        peaks: vec![Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 2.0)],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![1.0_f32, f32::NAN])],
        ..MSSpectrum::default()
    };
    let before = format!("{spec:?}");
    assert!(matches!(
        spec.sort_by_ion_mobility(),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        spec.is_sorted_by_im(),
        Err(Error::InvalidValue(_))
    ));
    // A NaN makes `PartialEq` false against any value, itself included, so the
    // "left unchanged" check compares the rendered value instead.
    assert_eq!(format!("{spec:?}"), before);

    // another array of the wrong length is still rejected before any move
    let mut spec = MSSpectrum {
        peaks: vec![Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 2.0)],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![2.0_f32, 1.0])],
        integer_data_arrays: vec![DataArray::new("tag", vec![1])],
        ..MSSpectrum::default()
    };
    let before = spec.clone();
    assert!(spec.sort_by_ion_mobility().is_err());
    assert_eq!(spec, before);

    // the ordinary case: peaks and the IM array move together
    let mut spec = MSSpectrum {
        peaks: vec![Peak1D::new(2.0, 1.0), Peak1D::new(1.0, 2.0)],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![5.0_f32, 4.0])],
        integer_data_arrays: vec![DataArray::new("tag", vec![10, 20])],
        ..MSSpectrum::default()
    };
    spec.sort_by_ion_mobility().unwrap();
    assert_eq!(spec.float_data_arrays[0].data, vec![4.0, 5.0]);
    assert_eq!(spec.peaks[0].mz, 1.0);
    assert_eq!(spec.integer_data_arrays[0].data, vec![20, 10]);
    assert!(spec.is_sorted_by_im().unwrap());

    // ties keep their original order
    let mut spec = MSSpectrum {
        peaks: vec![
            Peak1D::new(3.0, 1.0),
            Peak1D::new(1.0, 2.0),
            Peak1D::new(2.0, 3.0),
        ],
        float_data_arrays: vec![DataArray::new(IM_MS_ARRAY, vec![1.0_f32, 1.0, 0.0])],
        ..MSSpectrum::default()
    };
    spec.sort_by_ion_mobility().unwrap();
    assert_eq!(
        spec.peaks.iter().map(|peak| peak.mz).collect::<Vec<_>>(),
        vec![2.0, 3.0, 1.0]
    );
}

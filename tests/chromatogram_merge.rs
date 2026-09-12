//! Section-by-section port of `MSChromatogram_test.cpp` (43 sections) and
//! `Mobilogram_test.cpp` (48 sections) at Core SDK
//! `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
//!
//! Every expected value is transcribed from those two class tests, which makes
//! this tier-3 evidence (source review) throughout. No C++ was built or run. The
//! per-section audit, including which Rust function covers which section, is in
//! `docs/CHROMATOGRAM_MERGE_SUPPORT.md` and `docs/MOBILOGRAM_SUPPORT.md`.

use openms::Error;
use openms::constants::user_param::MERGED_CHROMATOGRAM_MZS;
use openms::kernel::chromatogram_merge::{
    ChromatogramMergeLimits, ChromatogramMergeOptions, ChromatogramSortLimits, MergedDataArrays,
    chromatogram_mz_less, merge_rt_key,
};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MobilityPeak1D, Mobilogram, NumericRange,
};
use openms::metadata::{DriftTimeUnit, MetaValue, MetaValueData, Product, ScanWindow};

fn peak(rt: f64, intensity: f32) -> ChromatogramPeak {
    ChromatogramPeak::new(rt, intensity)
}

fn mobility_peak(mobility: f64, intensity: f32) -> MobilityPeak1D {
    MobilityPeak1D::new(mobility, intensity)
}

fn text_array(name: &str, values: &[&str]) -> DataArray<String> {
    DataArray::new(name, values.iter().map(|v| (*v).to_owned()).collect())
}

// MSChromatogram_test.cpp:47-66 -- the shared p1..p5 fixture peaks.
const P1: ChromatogramPeak = ChromatogramPeak::new(2.0, 1.0);
const P2: ChromatogramPeak = ChromatogramPeak::new(10.0, 2.0);
const P3: ChromatogramPeak = ChromatogramPeak::new(30.0, 3.0);
const P4: ChromatogramPeak = ChromatogramPeak::new(30.0, 6.0);
const P5: ChromatogramPeak = ChromatogramPeak::new(30.0001, 3.0);

// MSChromatogram_test.cpp:504-529 -- the 21-point search fixture.
fn search_fixture() -> MSChromatogram {
    MSChromatogram::from_peaks(vec![
        peak(412.321, 29.0),
        peak(412.824, 60.0),
        peak(413.8, 34.0),
        peak(414.301, 29.0),
        peak(415.287, 37.0),
        peak(416.293, 31.0),
        peak(418.232, 31.0),
        peak(419.113, 31.0),
        peak(420.13, 201.0),
        peak(423.269, 56.0),
        peak(426.292, 34.0),
        peak(427.28, 82.0),
        peak(428.322, 87.0),
        peak(430.269, 30.0),
        peak(431.246, 29.0),
        peak(432.289, 42.0),
        peak(436.161, 32.0),
        peak(437.219, 54.0),
        peak(439.186, 40.0),
        peak(440.27, 40.0),
        peak(441.224, 23.0),
    ])
}

// MSChromatogram_test.cpp:549-567 and 792-799 -- seven points at RT 1..7.
fn ladder() -> MSChromatogram {
    MSChromatogram::from_peaks((1..=7).map(|rt| peak(f64::from(rt), 0.0)).collect())
}

// MSChromatogram_test.cpp:354-371 -- RT order 3,1,2 with mirroring arrays and no
// float array, the shape that used to defeat the fast-path guard.
fn annotated() -> MSChromatogram {
    MSChromatogram {
        peaks: vec![peak(3.0, 10.0), peak(1.0, 30.0), peak(2.0, 20.0)],
        string_data_arrays: vec![text_array("s1", &["rt3", "rt1", "rt2"])],
        integer_data_arrays: vec![DataArray::new("i1", vec![3, 1, 2])],
        ..MSChromatogram::default()
    }
}

// MSChromatogram_test.cpp:450-466 -- four points with tied RTs and tied
// intensities, carrying their original identity A..D / 1..4.
fn tied() -> MSChromatogram {
    MSChromatogram {
        peaks: vec![
            peak(1.0, 10.0),
            peak(2.0, 20.0),
            peak(1.0, 20.0),
            peak(2.0, 10.0),
        ],
        string_data_arrays: vec![text_array("s", &["A", "B", "C", "D"])],
        integer_data_arrays: vec![DataArray::new("i", vec![1, 2, 3, 4])],
        ..MSChromatogram::default()
    }
}

fn identities(chromatogram: &MSChromatogram) -> (Vec<&str>, Vec<i32>) {
    (
        chromatogram.string_data_arrays[0]
            .data
            .iter()
            .map(String::as_str)
            .collect(),
        chromatogram.integer_data_arrays[0].data.clone(),
    )
}

// --------------------------------------------------------------------------
// MSChromatogram_test.cpp
// --------------------------------------------------------------------------

/// Sections `MSChromatogram()` (34), `virtual ~MSChromatogram()` (41),
/// `getName() const` (69) and `setName` (78).
#[test]
fn source_construction_destruction_and_name() {
    let mut chromatogram = MSChromatogram::new();
    assert_eq!(chromatogram, MSChromatogram::default());
    // The C++ section only checks that `new MSChromatogram()` is not null and
    // that `delete` runs; ownership is static here, so the counterpart is that a
    // value constructs and drops without any explicit step.
    assert_eq!(chromatogram.name, "");
    chromatogram.name = "my_fancy_name".into();
    assert_eq!(chromatogram.name, "my_fancy_name");
    drop(chromatogram);
}

/// Sections `const FloatDataArrays& getFloatDataArrays() const` (84) through
/// `IntegerDataArrays& getIntegerDataArrays()` (121): the six accessors read
/// empty and the mutable ones can be resized in place.
#[test]
fn source_data_array_accessors() {
    let chromatogram = MSChromatogram::new();
    assert_eq!(chromatogram.float_data_arrays.len(), 0);
    assert_eq!(chromatogram.string_data_arrays.len(), 0);
    assert_eq!(chromatogram.integer_data_arrays.len(), 0);

    let mut chromatogram = MSChromatogram::new();
    chromatogram
        .float_data_arrays
        .resize(2, DataArray::default());
    assert_eq!(chromatogram.float_data_arrays.len(), 2);
    let mut chromatogram = MSChromatogram::new();
    chromatogram
        .string_data_arrays
        .resize(2, DataArray::default());
    assert_eq!(chromatogram.string_data_arrays.len(), 2);
    let mut chromatogram = MSChromatogram::new();
    chromatogram
        .integer_data_arrays
        .resize(2, DataArray::default());
    assert_eq!(chromatogram.integer_data_arrays.len(), 2);
}

/// Section `void sortByIntensity(bool reverse=false)` (129).
#[test]
fn source_sort_by_intensity_permutes_every_array() {
    let rts = [
        420.130, 412.824, 423.269, 415.287, 413.800, 419.113, 416.293, 418.232, 414.301, 412.321,
    ];
    let intensities = [201.0, 60.0, 56.0, 37.0, 34.0, 31.0, 31.0, 31.0, 29.0, 29.0];
    let strings = [
        "420.13", "412.82", "423.27", "415.29", "413.80", "419.11", "416.29", "418.23", "414.30",
        "412.32",
    ];
    let ints = [420, 412, 423, 415, 413, 419, 416, 418, 414, 412];

    // First half of the section: no data arrays, intensities come out ascending.
    let mut plain = MSChromatogram::from_peaks(
        (0..rts.len())
            .map(|i| peak(rts[i], intensities[i]))
            .collect(),
    );
    plain.sort_by_intensity(false).unwrap();
    let mut ascending = intensities;
    ascending.sort_by(f32::total_cmp);
    for (point, expected) in plain.peaks.iter().zip(ascending) {
        assert_eq!(point.intensity, expected);
    }

    // Second half: three float, two string and one integer array, all named and
    // all mirroring their point, must follow the same permutation.
    let mut chromatogram = MSChromatogram::from_peaks(
        (0..rts.len())
            .map(|i| peak(rts[i], intensities[i]))
            .collect(),
    );
    for name in ["f1", "f2", "f3"] {
        chromatogram.float_data_arrays.push(DataArray::new(
            name,
            rts.iter().map(|&r| r as f32).collect(),
        ));
    }
    for name in ["s1", "s2"] {
        chromatogram
            .string_data_arrays
            .push(text_array(name, &strings));
    }
    chromatogram
        .integer_data_arrays
        .push(DataArray::new("i1", ints.to_vec()));
    chromatogram.sort_by_intensity(false).unwrap();

    for (index, name) in ["f1", "f2", "f3"].iter().enumerate() {
        assert_eq!(chromatogram.float_data_arrays[index].name, *name);
    }
    assert_eq!(chromatogram.string_data_arrays[0].name, "s1");
    assert_eq!(chromatogram.string_data_arrays[1].name, "s2");
    assert_eq!(chromatogram.integer_data_arrays[0].name, "i1");

    // Stable ascending intensity moves original indices 8,9,5,6,7,4,3,2,1,0.
    let expected_rts = [
        414.301, 412.321, 419.113, 416.293, 418.232, 413.800, 415.287, 423.269, 412.824, 420.130,
    ];
    let expected_strings = [
        "414.30", "412.32", "419.11", "416.29", "418.23", "413.80", "415.29", "423.27", "412.82",
        "420.13",
    ];
    let expected_ints = [414, 412, 419, 416, 418, 413, 415, 423, 412, 420];
    for i in 0..rts.len() {
        assert_eq!(chromatogram.peaks[i].intensity, ascending[i]);
        assert_eq!(chromatogram.peaks[i].rt, expected_rts[i]);
        // The C++ section compares the float array against the point's own RT
        // with TOLERANCE_ABSOLUTE(0.0001); f32 rounding at 420 is ~3e-5.
        let stored = f64::from(chromatogram.float_data_arrays[1].data[i]);
        assert!((stored - chromatogram.peaks[i].rt).abs() < 1e-4);
        assert_eq!(
            chromatogram.string_data_arrays[0].data[i],
            expected_strings[i]
        );
        assert_eq!(
            chromatogram.integer_data_arrays[0].data[i],
            expected_ints[i]
        );
    }
}

/// Section `void sortByPosition()` (229).
#[test]
fn source_sort_by_position_permutes_every_array() {
    let rts = [
        423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321,
    ];
    let intensities = [
        56.0f32, 201.0, 31.0, 31.0, 31.0, 37.0, 29.0, 34.0, 60.0, 29.0,
    ];
    let strings = ["56", "201", "31", "31", "31", "37", "29", "34", "60", "29"];
    let ints = [56, 201, 31, 31, 31, 37, 29, 34, 60, 29];

    let mut plain = MSChromatogram::from_peaks(
        (0..rts.len())
            .map(|i| peak(rts[i], intensities[i]))
            .collect(),
    );
    plain.sort_by_position().unwrap();
    // The input is strictly descending in RT, so sorting reverses it.
    for (point, expected) in plain.peaks.iter().zip(intensities.iter().rev()) {
        assert_eq!(point.intensity, *expected);
    }

    let mut chromatogram = MSChromatogram::from_peaks(
        (0..rts.len())
            .map(|i| peak(rts[i], intensities[i]))
            .collect(),
    );
    for name in ["f1", "f2", "f3"] {
        chromatogram
            .float_data_arrays
            .push(DataArray::new(name, intensities.to_vec()));
    }
    for name in ["s1", "s2"] {
        chromatogram
            .string_data_arrays
            .push(text_array(name, &strings));
    }
    // The C++ section installs two integer arrays and names only the first.
    for _ in 0..2 {
        chromatogram
            .integer_data_arrays
            .push(DataArray::new("", ints.to_vec()));
    }
    chromatogram.integer_data_arrays[0].name = "i1".into();
    chromatogram.sort_by_position().unwrap();

    for (index, name) in ["f1", "f2", "f3"].iter().enumerate() {
        assert_eq!(chromatogram.float_data_arrays[index].name, *name);
    }
    assert_eq!(chromatogram.string_data_arrays[0].name, "s1");
    assert_eq!(chromatogram.string_data_arrays[1].name, "s2");
    assert_eq!(chromatogram.integer_data_arrays[0].name, "i1");
    for (i, expected) in intensities.iter().rev().enumerate() {
        assert_eq!(chromatogram.peaks[i].intensity, *expected);
        assert_eq!(chromatogram.float_data_arrays[1].data[i], *expected);
        assert_eq!(
            chromatogram.string_data_arrays[0].data[i],
            format!("{expected:.0}")
        );
        assert_eq!(
            chromatogram.integer_data_arrays[0].data[i],
            *expected as i32
        );
    }
}

/// Section `bool isSorted() const` (320).
#[test]
fn source_is_sorted() {
    let mut chromatogram = MSChromatogram::from_peaks(vec![
        peak(1000.0, 1.0),
        peak(1001.0, 1.0),
        peak(1002.0, 1.0),
    ]);
    assert!(chromatogram.is_sorted());
    chromatogram.peaks.reverse();
    assert!(!chromatogram.is_sorted());
}

/// Section `[EXTRA] sorting reorders string and integer data arrays even when no
/// float data array is present` (345), the 47-assertion regression section.
#[test]
fn source_sorting_and_selection_regression() {
    // sortByPosition: RT 1,2,3 and the arrays follow.
    let mut chromatogram = annotated();
    assert!(chromatogram.float_data_arrays.is_empty());
    chromatogram.sort_by_position().unwrap();
    assert_eq!(chromatogram.peaks[0].rt, 1.0);
    assert_eq!(chromatogram.peaks[1].rt, 2.0);
    assert_eq!(chromatogram.peaks[2].rt, 3.0);
    assert_eq!(
        chromatogram.string_data_arrays[0].data,
        ["rt1", "rt2", "rt3"]
    );
    assert_eq!(chromatogram.integer_data_arrays[0].data, [1, 2, 3]);
    assert_eq!(chromatogram.string_data_arrays[0].name, "s1");
    assert_eq!(chromatogram.integer_data_arrays[0].name, "i1");

    // sortByIntensity ascending: 10,20,30 -> RT 3,2,1.
    let mut chromatogram = annotated();
    chromatogram.sort_by_intensity(false).unwrap();
    assert_eq!(chromatogram.peaks[0].intensity, 10.0);
    assert_eq!(chromatogram.peaks[1].intensity, 20.0);
    assert_eq!(chromatogram.peaks[2].intensity, 30.0);
    assert_eq!(
        chromatogram.string_data_arrays[0].data,
        ["rt3", "rt2", "rt1"]
    );
    assert_eq!(chromatogram.integer_data_arrays[0].data, [3, 2, 1]);

    // sortByIntensity reverse: 30,20,10 -> RT 1,2,3.
    let mut chromatogram = annotated();
    chromatogram.sort_by_intensity(true).unwrap();
    assert_eq!(chromatogram.peaks[0].intensity, 30.0);
    assert_eq!(chromatogram.peaks[2].intensity, 10.0);
    assert_eq!(chromatogram.string_data_arrays[0].data[0], "rt1");
    assert_eq!(chromatogram.string_data_arrays[0].data[2], "rt3");
    assert_eq!(chromatogram.integer_data_arrays[0].data[0], 1);
    assert_eq!(chromatogram.integer_data_arrays[0].data[2], 3);

    // A float array present must be permuted too.
    let mut chromatogram = annotated();
    chromatogram
        .float_data_arrays
        .push(DataArray::new("f", vec![3.5, 1.5, 2.5]));
    chromatogram.sort_by_intensity(true).unwrap();
    assert_eq!(chromatogram.float_data_arrays[0].data, [1.5, 2.5, 3.5]);

    // A mis-sized array is rejected and nothing changes.
    let mut chromatogram = annotated();
    chromatogram.integer_data_arrays[0].data.push(99);
    let before = chromatogram.clone();
    assert!(chromatogram.sort_by_position().is_err());
    assert_eq!(chromatogram, before);
    assert_eq!(chromatogram.len(), 3);
    assert_eq!(chromatogram.peaks[0].rt, 3.0);
    assert_eq!(chromatogram.string_data_arrays[0].data[0], "rt3");

    // Out-of-range selection indices are rejected, not undefined.
    let mut chromatogram = annotated();
    let before = chromatogram.clone();
    assert!(chromatogram.select(&[0, 7]).is_err());
    assert_eq!(chromatogram, before);
    assert_eq!(chromatogram.len(), 3);

    // A valid unique reorder keeps the arrays aligned.
    let mut chromatogram = annotated();
    chromatogram.select(&[2, 0]).unwrap();
    assert_eq!(chromatogram.len(), 2);
    assert_eq!(chromatogram.string_data_arrays[0].data, ["rt2", "rt3"]);

    // Stability across tied keys, for all three sorts.
    let mut chromatogram = tied();
    chromatogram.sort_by_position().unwrap();
    assert_eq!(chromatogram.peaks[0].rt, 1.0);
    assert_eq!(chromatogram.peaks[1].rt, 1.0);
    assert_eq!(chromatogram.peaks[2].rt, 2.0);
    assert_eq!(chromatogram.peaks[3].rt, 2.0);
    assert_eq!(
        identities(&chromatogram),
        (vec!["A", "C", "B", "D"], vec![1, 3, 2, 4])
    );

    let mut chromatogram = tied();
    chromatogram.sort_by_intensity(false).unwrap();
    assert_eq!(
        identities(&chromatogram),
        (vec!["A", "D", "B", "C"], vec![1, 4, 2, 3])
    );

    let mut chromatogram = tied();
    chromatogram.sort_by_intensity(true).unwrap();
    assert_eq!(
        identities(&chromatogram),
        (vec!["B", "C", "A", "D"], vec![2, 3, 1, 4])
    );

    // sort(lambda) validates array sizes before the predicate is ever invoked,
    // so an undersized array is an error rather than an out-of-bounds read.
    let mut chromatogram = annotated();
    chromatogram.integer_data_arrays[0].data.pop();
    let before = chromatogram.clone();
    let mut calls = 0usize;
    let rejected = chromatogram.sort_by(|c, i, j| {
        calls += 1;
        c.integer_data_arrays[0].data[i] < c.integer_data_arrays[0].data[j]
    });
    assert!(rejected.is_err());
    assert_eq!(calls, 0);
    assert_eq!(chromatogram, before);
    assert_eq!(chromatogram.len(), 3);
}

/// The rest of source `template<class Predicate> void sort(const Predicate&)`
/// (`MSChromatogram.h:240`), whose ordering contract the class test exercises
/// only through the rejection case above. Modelled on the equivalent
/// `Mobilogram_test.cpp:614-646` section, which sorts by a data array.
#[test]
fn predicate_sort_follows_a_data_array_and_is_bounded() {
    let mut chromatogram = MSChromatogram {
        peaks: vec![peak(1.0, 10.0), peak(2.0, 20.0), peak(3.0, 30.0)],
        integer_data_arrays: vec![DataArray::new("rank", vec![2, 0, 1])],
        string_data_arrays: vec![text_array("tag", &["a", "b", "c"])],
        ..MSChromatogram::default()
    };
    chromatogram
        .sort_by(|c, i, j| c.integer_data_arrays[0].data[i] < c.integer_data_arrays[0].data[j])
        .unwrap();
    assert_eq!(chromatogram.integer_data_arrays[0].data, [0, 1, 2]);
    assert_eq!(chromatogram.peaks[0].rt, 2.0);
    assert_eq!(chromatogram.peaks[1].rt, 3.0);
    assert_eq!(chromatogram.peaks[2].rt, 1.0);
    assert_eq!(chromatogram.string_data_arrays[0].data, ["b", "c", "a"]);

    for limits in [
        ChromatogramSortLimits {
            max_peaks: 2,
            ..Default::default()
        },
        ChromatogramSortLimits {
            max_bytes: 0,
            ..Default::default()
        },
    ] {
        let mut candidate = chromatogram.clone();
        let before = candidate.clone();
        let mut calls = 0usize;
        assert!(
            candidate
                .sort_by_with_limits(
                    |_, _, _| {
                        calls += 1;
                        false
                    },
                    limits
                )
                .is_err()
        );
        assert_eq!(calls, 0);
        assert_eq!(candidate, before);
    }
}

/// Section `Size findNearest(CoordinateType rt) const` (504).
#[test]
fn source_find_nearest() {
    let chromatogram = search_fixture();
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
        assert_eq!(chromatogram.find_nearest(query).unwrap(), Some(index));
    }
    // The source throws Exception::Precondition on an empty chromatogram; this
    // port reports absence instead.
    assert_eq!(MSChromatogram::new().find_nearest(427.3).unwrap(), None);
}

/// Sections `Iterator RTBegin(CoordinateType rt)` (549) through
/// `ConstIterator RTEnd(ConstIterator, CoordinateType, ConstIterator) const`
/// (762). C++ has eight overloads for mutable/const and whole/subrange; Rust has
/// one shared pair per bound, so each overload's literals are checked against it.
#[test]
fn source_rt_bounds_and_subranges() {
    let chromatogram = ladder();
    let rt_at = |index: usize| chromatogram.peaks[index].rt;

    // RTBegin(rt), both the mutable and the const overload.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_begin(4.5).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_begin(5.0).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_begin(5.5).unwrap()), 6.0);
    }
    // RTEnd(rt), both overloads.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_end(4.5).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_end(5.0).unwrap()), 6.0);
        assert_eq!(rt_at(chromatogram.rt_end(5.5).unwrap()), 6.0);
    }
    // RTBegin(begin, rt, end), both overloads. An empty subrange returns its
    // own start, as lower_bound(begin, .., begin) returns begin.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_begin_in(4.5, 0..7).unwrap()), 5.0);
        assert_eq!(chromatogram.rt_begin_in(4.5, 0..0).unwrap(), 0);
    }
    // RTEnd(begin, rt, end), both overloads.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_end_in(3.5, 0..7).unwrap()), 4.0);
        assert_eq!(rt_at(chromatogram.rt_end_in(5.0, 0..7).unwrap()), 6.0);
        assert_eq!(chromatogram.rt_end_in(4.5, 0..0).unwrap(), 0);
    }

    // Native checks the source leaves undefined.
    assert!(chromatogram.rt_begin_in(4.5, 4..99).is_err());
    let inverted = std::ops::Range { start: 5, end: 2 };
    assert!(chromatogram.rt_end_in(4.5, inverted).is_err());
    assert!(chromatogram.rt_begin(f64::NAN).is_err());
    let unsorted = MSChromatogram::from_peaks(vec![peak(2.0, 0.0), peak(1.0, 0.0)]);
    assert!(matches!(unsorted.rt_begin(1.5), Err(Error::UnsortedData)));
    assert!(matches!(
        unsorted.rt_begin_in(1.5, 0..2),
        Err(Error::UnsortedData)
    ));
}

/// Sections `Iterator PosBegin(CoordinateType rt)` (801) through
/// `ConstIterator PosEnd(ConstIterator, CoordinateType, ConstIterator) const`
/// (891). The source declares these as aliases of RTBegin()/RTEnd(); this checks
/// each section's literals and that the alias really is the same answer.
#[test]
fn source_pos_bounds_are_rt_bounds() {
    let chromatogram = ladder();
    let rt_at = |index: usize| chromatogram.peaks[index].rt;
    let last = chromatogram.len() - 1;

    // PosBegin(rt) and its const overload.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_begin(4.5).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_begin(5.0).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_begin(5.5).unwrap()), 6.0);
    }
    // PosBegin(begin, rt, end), mutable overload.
    assert_eq!(rt_at(chromatogram.rt_begin_in(4.5, 0..7).unwrap()), 5.0);
    assert_eq!(rt_at(chromatogram.rt_begin_in(5.5, 0..7).unwrap()), 6.0);
    assert_eq!(chromatogram.rt_begin_in(4.5, 0..0).unwrap(), 0);
    // A query past the last point lands on the past-the-end index.
    assert_eq!(chromatogram.rt_begin_in(8.0, 0..7).unwrap(), 7);
    assert_eq!(
        rt_at(chromatogram.rt_begin_in(8.0, 0..7).unwrap() - 1),
        rt_at(last)
    );
    // PosBegin(begin, rt, end), const overload.
    assert_eq!(rt_at(chromatogram.rt_begin_in(3.5, 0..7).unwrap()), 4.0);
    assert_eq!(rt_at(chromatogram.rt_begin_in(4.5, 0..7).unwrap()), 5.0);

    // PosEnd(rt) and its const overload.
    for _ in 0..2 {
        assert_eq!(rt_at(chromatogram.rt_end(4.5).unwrap()), 5.0);
        assert_eq!(rt_at(chromatogram.rt_end(5.0).unwrap()), 6.0);
        assert_eq!(rt_at(chromatogram.rt_end(5.5).unwrap()), 6.0);
    }
    // PosEnd(begin, rt, end), mutable overload.
    assert_eq!(rt_at(chromatogram.rt_end_in(3.5, 0..7).unwrap()), 4.0);
    assert_eq!(rt_at(chromatogram.rt_end_in(4.0, 0..7).unwrap()), 5.0);
    assert_eq!(chromatogram.rt_end_in(4.5, 0..0).unwrap(), 0);
    // PosEnd(begin, rt, end), const overload.
    assert_eq!(rt_at(chromatogram.rt_end_in(4.5, 0..7).unwrap()), 5.0);
    assert_eq!(rt_at(chromatogram.rt_end_in(5.0, 0..7).unwrap()), 6.0);

    // The alias identity itself, over every query the sections use.
    for query in [3.5, 4.0, 4.5, 5.0, 5.5, 8.0] {
        assert_eq!(
            chromatogram.rt_begin_in(query, 0..7).unwrap(),
            chromatogram.rt_begin(query).unwrap()
        );
        assert_eq!(
            chromatogram.rt_end_in(query, 0..7).unwrap(),
            chromatogram.rt_end(query).unwrap()
        );
    }
}

// MSChromatogram_test.cpp:908-921 -- the record used by all four copy/move
// sections: one scan window, a "label" meta value, product m/z, a name and a
// single point.
fn record() -> MSChromatogram {
    let mut chromatogram = MSChromatogram {
        product: Product {
            mz: 7.0,
            ..Product::default()
        },
        name: "bla".into(),
        peaks: vec![peak(47.11, 0.0)],
        ..MSChromatogram::default()
    };
    chromatogram
        .instrument_settings
        .scan_windows
        .push(ScanWindow::default());
    chromatogram.metadata.insert(
        "label".into(),
        MetaValue::new(MetaValueData::Float(5.0)).unwrap(),
    );
    chromatogram
}

/// Sections `MSChromatogram(const MSChromatogram&)` (908),
/// `MSChromatogram(const MSChromatogram&&)` (934),
/// `operator=(const MSChromatogram&)` (972) and
/// `operator=(const MSChromatogram&&)` (1005).
#[test]
fn source_copy_move_and_assignment() {
    // Copy construction.
    let original = record();
    let copy = original.clone();
    assert_eq!(copy.instrument_settings.scan_windows.len(), 1);
    assert_eq!(copy.metadata["label"].as_f64().unwrap(), 5.0);
    assert_eq!(copy.product.mz, 7.0);
    assert_eq!(copy.name, "bla");
    assert_eq!(copy.len(), 1);
    assert_eq!(copy.peaks[0].rt, 47.11);

    // Move construction. The C++ section additionally asserts the move
    // constructor is noexcept and that the moved-from value is empty; Rust moves
    // never panic and statically forbid touching the source afterwards, so those
    // two have no runtime counterpart.
    let original = record();
    let reference = original.clone();
    let moved = original;
    assert_eq!(moved, reference);
    assert_eq!(moved.instrument_settings.scan_windows.len(), 1);
    assert_eq!(moved.metadata["label"].as_f64().unwrap(), 5.0);
    assert_eq!(moved.product.mz, 7.0);
    assert_eq!(moved.name, "bla");
    assert_eq!(moved.len(), 1);
    assert_eq!(moved.peaks[0].rt, 47.11);

    // Copy assignment, then assignment of a default value. The C++ record for
    // the assignment sections has no scan window.
    let mut source = record();
    source.instrument_settings.scan_windows.clear();
    let mut target = MSChromatogram::new();
    assert_eq!(target.len(), 0);
    target = source.clone();
    assert_eq!(target.metadata["label"].as_f64().unwrap(), 5.0);
    assert_eq!(target.product.mz, 7.0);
    assert_eq!(target.name, "bla");
    assert_eq!(target.len(), 1);
    assert_eq!(target.peaks[0].rt, 47.11);
    target = MSChromatogram::new();
    assert_eq!(target.instrument_settings.scan_windows.len(), 0);
    assert!(!target.metadata.contains_key("label"));
    assert_eq!(target.product.mz, 0.0);
    assert_eq!(target.name, "");
    assert_eq!(target.len(), 0);

    // Move assignment.
    let mut target = MSChromatogram::new();
    assert_eq!(target.len(), 0);
    target = source.clone();
    assert_eq!(target, source);
    assert_eq!(target.metadata["label"].as_f64().unwrap(), 5.0);
    assert_eq!(target.product.mz, 7.0);
    assert_eq!(target.name, "bla");
    assert_eq!(target.len(), 1);
    assert_eq!(target.peaks[0].rt, 47.11);
    target = MSChromatogram::default();
    assert_eq!(target.instrument_settings.scan_windows.len(), 0);
    assert!(!target.metadata.contains_key("label"));
    assert_eq!(target.product.mz, 0.0);
    assert_eq!(target.name, "");
    assert_eq!(target.len(), 0);
}

/// Sections `bool operator==(const MSChromatogram&) const` (1057) and
/// `bool operator!=(const MSChromatogram&) const` (1103).
#[test]
fn source_equality_ignores_only_the_name() {
    let empty = MSChromatogram::new();
    let edit = MSChromatogram::new();
    assert!(edit.source_equal(&empty));
    assert!(edit == empty);

    let mut edit = empty.clone();
    edit.instrument_settings
        .scan_windows
        .resize(1, ScanWindow::default());
    assert!(!edit.source_equal(&empty));

    let mut edit = empty.clone();
    edit.peaks.resize(1, ChromatogramPeak::default());
    assert!(!edit.source_equal(&empty));

    let mut edit = empty.clone();
    edit.metadata.insert("label".into(), "bla".into());
    assert!(!empty.source_equal(&edit));
    edit.product = Product {
        mz: 5.0,
        ..Product::default()
    };
    assert!(!empty.source_equal(&edit));

    let mut edit = empty.clone();
    edit.float_data_arrays.resize(5, DataArray::default());
    assert!(!empty.source_equal(&edit));

    let mut edit = empty.clone();
    edit.string_data_arrays.resize(5, DataArray::default());
    assert!(!empty.source_equal(&edit));

    let mut edit = empty.clone();
    edit.integer_data_arrays.resize(5, DataArray::default());
    assert!(!empty.source_equal(&edit));

    // The name is not part of the source predicate. The derived PartialEq does
    // compare it, which is the one deliberate difference.
    let mut edit = empty.clone();
    edit.name = "bla".into();
    assert!(empty.source_equal(&edit));
    assert!(empty != edit);

    // Points plus a range computation, then clear(false), is equal to empty
    // again: this port has no range cache to leave behind.
    let mut edit = empty.clone();
    edit.peaks.push(P1);
    edit.peaks.push(P2);
    edit.range_manager().unwrap();
    edit.clear(false);
    assert!(empty.source_equal(&edit));
    assert!(empty == edit);
}

/// Section `virtual void updateRanges()` (1152). The source caches the ranges in
/// the inherited RangeManager; this port computes them on demand, so "call it
/// twice to check the initialization" becomes "two calls agree".
#[test]
fn source_update_ranges() {
    let mut chromatogram = MSChromatogram::from_peaks(vec![P1, P2, P1]);
    let first = chromatogram.range_manager().unwrap();
    let second = chromatogram.range_manager().unwrap();
    assert_eq!(first, second);
    assert_eq!(second.max_intensity().unwrap(), 2.0);
    assert_eq!(second.min_intensity().unwrap(), 1.0);
    assert_eq!(second.max_rt().unwrap(), 10.0);
    assert_eq!(second.min_rt().unwrap(), 2.0);
    assert_eq!(
        chromatogram.ranges().unwrap().rt,
        Some(NumericRange {
            min: 2.0,
            max: 10.0
        })
    );

    chromatogram.clear(true);
    chromatogram.peaks.push(P1);
    let ranges = chromatogram.range_manager().unwrap();
    assert_eq!(ranges.max_intensity().unwrap(), 1.0);
    assert_eq!(ranges.min_intensity().unwrap(), 1.0);
    assert_eq!(ranges.max_rt().unwrap(), 2.0);
    assert_eq!(ranges.min_rt().unwrap(), 2.0);
}

/// Section `void clear(bool clear_meta_data)` (1179).
#[test]
fn source_clear() {
    let mut edit = MSChromatogram::new();
    edit.instrument_settings
        .scan_windows
        .resize(1, ScanWindow::default());
    edit.peaks.resize(1, ChromatogramPeak::default());
    edit.metadata.insert("label".into(), "bla".into());
    edit.product.mz = 5.0;
    edit.float_data_arrays.resize(5, DataArray::default());
    edit.integer_data_arrays.resize(5, DataArray::default());
    edit.string_data_arrays.resize(5, DataArray::default());

    edit.clear(false);
    assert_eq!(edit.len(), 0);
    assert!(edit.is_empty());
    assert!(!edit.source_equal(&MSChromatogram::new()));
    // The arrays are parallel to the points, so dropping the points drops them.
    assert!(edit.float_data_arrays.is_empty());
    assert!(edit.integer_data_arrays.is_empty());
    assert!(edit.string_data_arrays.is_empty());

    edit.clear(true);
    assert_eq!(edit.len(), 0);
    assert!(edit.is_empty());
    assert!(edit.source_equal(&MSChromatogram::new()));
    assert_eq!(edit, MSChromatogram::new());
}

/// Sections `double getMZ() const` (1207) and
/// `[MSChromatogram::MZLess] bool operator()(...) const` (1218).
#[test]
fn source_product_mz_and_mz_less() {
    let mut chromatogram = MSChromatogram::new();
    assert_eq!(chromatogram.product.mz, 0.0);
    chromatogram.product = Product {
        mz: 0.1,
        ..Product::default()
    };
    assert_eq!(chromatogram.product.mz, 0.1);

    let mut a = MSChromatogram::new();
    a.product.mz = 1000.0;
    let mut b = MSChromatogram::new();
    b.product.mz = 1000.1;
    assert!(chromatogram_mz_less(&a, &b));
    assert!(!chromatogram_mz_less(&b, &a));
    assert!(!chromatogram_mz_less(&a, &a));
}

/// Section `void mergePeaks(MSChromatogram& other)` (1237).
#[test]
fn source_merge_peaks() {
    let mut a = MSChromatogram::from_peaks(vec![P1, P3]);
    let b = MSChromatogram::from_peaks(vec![P2, P5]);
    let mut expected = MSChromatogram::from_peaks(vec![P1, P2, P4]);
    a.sort_by_position().unwrap();
    expected.sort_by_position().unwrap();

    a.merge_peaks(&b, true).unwrap();
    expected.metadata.insert(
        MERGED_CHROMATOGRAM_MZS.to_owned(),
        MetaValue::new(MetaValueData::FloatList(vec![b.product.mz])).unwrap(),
    );
    assert!(a.source_equal(&expected));
    assert_eq!(a, expected);
    // RT 30.0 and 30.0001 share millisecond bucket 30000, so their intensities
    // 3 and 3 are summed into the single point (30.0, 6.0) that `p4` describes.
    assert_eq!(a.peaks, vec![P1, P2, P4]);
    assert_eq!(a.merged_chromatogram_mzs().unwrap(), [0.0]);
}

/// Section `std::ostream& operator<<(std::ostream&, const MSChromatogram&)`
/// (1257).
#[test]
fn source_stream_layout() {
    let mut chromatogram = MSChromatogram::new();
    chromatogram.product.mz = 1000.0;
    chromatogram.peaks.push(peak(47.11, 0.0));
    let rendered = chromatogram.to_string();
    assert!(rendered.contains("MSCHROMATOGRAM BEGIN"));
    assert!(rendered.contains("47.11"));
    // The full source layout: the settings operator prints only its own two
    // delimiters and ignores the settings themselves.
    assert_eq!(
        rendered,
        "-- MSCHROMATOGRAM BEGIN --\n\
         -- CHROMATOGRAMSETTINGS BEGIN --\n\
         -- CHROMATOGRAMSETTINGS END --\n\
         POS: 47.11 INT: 0\n\
         -- MSCHROMATOGRAM END --\n"
    );
    assert_eq!(
        format!("{chromatogram:.2}"),
        "-- MSCHROMATOGRAM BEGIN --\n\
         -- CHROMATOGRAMSETTINGS BEGIN --\n\
         -- CHROMATOGRAMSETTINGS END --\n\
         POS: 47.11 INT: 0.00\n\
         -- MSCHROMATOGRAM END --\n"
    );
}

// --------------------------------------------------------------------------
// Native merge behaviour beyond the class test's single section
// --------------------------------------------------------------------------

/// The millisecond bucketing that decides which points are summed, and the
/// resulting non-distance behaviour recorded at
/// [`merge_rt_key`](openms::kernel::chromatogram_merge::merge_rt_key).
#[test]
fn merge_key_is_a_millisecond_bucket_not_a_distance() {
    assert_eq!(merge_rt_key(0.0), 0.0);
    assert_eq!(merge_rt_key(0.0005), 1.0); // round() breaks halves away from zero
    assert_eq!(merge_rt_key(0.00149), 1.0);
    assert_eq!(merge_rt_key(0.00151), 2.0);
    assert_eq!(merge_rt_key(-0.0005), -1.0);

    // 0.00099 s apart and merged ...
    let mut near = MSChromatogram::from_peaks(vec![peak(0.0005, 1.0)]);
    near.merge_peaks(&MSChromatogram::from_peaks(vec![peak(0.00149, 2.0)]), false)
        .unwrap();
    assert_eq!(near.peaks, vec![peak(0.0005, 3.0)]);

    // ... while 0.00002 s apart stay separate.
    let mut far = MSChromatogram::from_peaks(vec![peak(0.00149, 1.0)]);
    far.merge_peaks(&MSChromatogram::from_peaks(vec![peak(0.00151, 2.0)]), false)
        .unwrap();
    assert_eq!(far.peaks, vec![peak(0.00149, 1.0), peak(0.00151, 2.0)]);
}

/// Draining either side, the destination's retention time winning a tie, and the
/// destination's own product m/z being left alone.
#[test]
fn merge_drains_both_sides_and_keeps_the_destination_product() {
    let mut chromatogram = MSChromatogram::from_peaks(vec![peak(1.0, 1.0), peak(2.0, 1.0)]);
    chromatogram.product.mz = 123.5;
    let other = MSChromatogram::from_peaks(vec![peak(3.0, 5.0), peak(4.0, 5.0)]);
    chromatogram.merge_peaks(&other, false).unwrap();
    assert_eq!(
        chromatogram.peaks,
        vec![
            peak(1.0, 1.0),
            peak(2.0, 1.0),
            peak(3.0, 5.0),
            peak(4.0, 5.0)
        ]
    );
    assert_eq!(chromatogram.product.mz, 123.5);

    // The other way round: `other` is drained first.
    let mut chromatogram = MSChromatogram::from_peaks(vec![peak(9.0, 1.0)]);
    chromatogram
        .merge_peaks(&MSChromatogram::from_peaks(vec![peak(1.0, 2.0)]), false)
        .unwrap();
    assert_eq!(chromatogram.peaks, vec![peak(1.0, 2.0), peak(9.0, 1.0)]);

    // A tie keeps the destination's retention time, not the other's.
    let mut chromatogram = MSChromatogram::from_peaks(vec![peak(5.0004, 1.0)]);
    chromatogram
        .merge_peaks(&MSChromatogram::from_peaks(vec![peak(5.0, 2.0)]), false)
        .unwrap();
    assert_eq!(chromatogram.peaks, vec![peak(5.0004, 3.0)]);

    // Merging with an empty chromatogram changes nothing; merging into an empty
    // one copies.
    let mut chromatogram = MSChromatogram::from_peaks(vec![peak(1.0, 1.0)]);
    chromatogram
        .merge_peaks(&MSChromatogram::new(), false)
        .unwrap();
    assert_eq!(chromatogram.peaks, vec![peak(1.0, 1.0)]);
    let mut chromatogram = MSChromatogram::new();
    chromatogram
        .merge_peaks(&MSChromatogram::from_peaks(vec![peak(1.0, 1.0)]), false)
        .unwrap();
    assert_eq!(chromatogram.peaks, vec![peak(1.0, 1.0)]);
}

/// `add_meta` appends to an existing list, and a non-list value is an error
/// where the source's `DataValue::toDoubleList` throws.
#[test]
fn merge_metadata_list_accumulates_and_rejects_a_wrong_type() {
    let mut chromatogram = MSChromatogram::new();
    assert!(chromatogram.merged_chromatogram_mzs().unwrap().is_empty());
    let mut other = MSChromatogram::from_peaks(vec![peak(1.0, 1.0)]);
    other.product.mz = 100.5;
    chromatogram.merge_peaks(&other, true).unwrap();
    other.product.mz = 200.25;
    chromatogram.merge_peaks(&other, true).unwrap();
    assert_eq!(
        chromatogram.merged_chromatogram_mzs().unwrap(),
        [100.5, 200.25]
    );

    let mut wrong = MSChromatogram::new();
    wrong
        .metadata
        .insert(MERGED_CHROMATOGRAM_MZS.to_owned(), "not a list".into());
    assert!(wrong.merged_chromatogram_mzs().is_err());
    let before = wrong.clone();
    assert!(wrong.merge_peaks(&other, true).is_err());
    assert_eq!(wrong, before);

    // Without add_meta nothing is recorded.
    let mut quiet = MSChromatogram::new();
    quiet.merge_peaks(&other, false).unwrap();
    assert!(quiet.merged_chromatogram_mzs().unwrap().is_empty());
}

/// The annotation-array policy: refuse, drop, or reproduce the source's
/// untouched-and-possibly-misaligned arrays.
#[test]
fn merge_annotation_array_policies() {
    let annotated_pair = || {
        let mut destination = MSChromatogram::from_peaks(vec![peak(1.0, 1.0)]);
        destination
            .integer_data_arrays
            .push(DataArray::new("i", vec![7]));
        let other = MSChromatogram::from_peaks(vec![peak(2.0, 1.0)]);
        (destination, other)
    };

    // Default: refused, nothing changed.
    let (mut destination, other) = annotated_pair();
    let before = destination.clone();
    assert!(matches!(
        destination.merge_peaks(&other, false),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(destination, before);

    // An array with no entries is a placeholder and does not block the merge; it
    // survives untouched.
    let mut placeholder = MSChromatogram::from_peaks(vec![peak(1.0, 1.0)]);
    placeholder
        .float_data_arrays
        .push(DataArray::new("empty", vec![]));
    placeholder.merge_peaks(&other, false).unwrap();
    assert_eq!(placeholder.len(), 2);
    assert_eq!(placeholder.float_data_arrays.len(), 1);
    assert!(placeholder.float_data_arrays[0].data.is_empty());

    // Drop: arrays are removed and the result validates.
    let (mut destination, other) = annotated_pair();
    destination
        .merge_peaks_with_options(
            &other,
            ChromatogramMergeOptions {
                data_arrays: MergedDataArrays::Drop,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(destination.len(), 2);
    assert!(destination.integer_data_arrays.is_empty());
    destination.validate().unwrap();

    // Source: arrays are left alone and the chromatogram no longer validates,
    // which is exactly the state the source header's @note warns about.
    let (mut destination, other) = annotated_pair();
    destination
        .merge_peaks_with_options(
            &other,
            ChromatogramMergeOptions::source().with_add_meta(true),
        )
        .unwrap();
    assert_eq!(destination.len(), 2);
    assert_eq!(destination.integer_data_arrays[0].data, [7]);
    assert!(destination.validate().is_err());
    assert!(destination.sort_by_position().is_err());
    assert_eq!(destination.merged_chromatogram_mzs().unwrap(), [0.0]);
}

/// Every rejection leaves both chromatograms untouched: unsorted input, an
/// invalid record, a non-finite summed intensity and each ceiling.
#[test]
fn merge_rejections_are_transactional() {
    let sorted = MSChromatogram::from_peaks(vec![peak(1.0, 1.0), peak(2.0, 1.0)]);

    let mut unsorted = MSChromatogram::from_peaks(vec![peak(2.0, 1.0), peak(1.0, 1.0)]);
    let before = unsorted.clone();
    assert!(matches!(
        unsorted.merge_peaks(&sorted, false),
        Err(Error::UnsortedData)
    ));
    assert_eq!(unsorted, before);

    let mut destination = sorted.clone();
    let before = destination.clone();
    assert!(matches!(
        destination.merge_peaks(&unsorted, false),
        Err(Error::UnsortedData)
    ));
    assert_eq!(destination, before);

    // A non-finite coordinate fails validate() before anything is built.
    let mut destination = sorted.clone();
    let before = destination.clone();
    assert!(
        destination
            .merge_peaks(
                &MSChromatogram::from_peaks(vec![peak(f64::NAN, 1.0)]),
                false
            )
            .is_err()
    );
    assert_eq!(destination, before);

    // The source lets a summed intensity overflow to infinity; this reports it.
    let mut destination = MSChromatogram::from_peaks(vec![peak(1.0, f32::MAX)]);
    let before = destination.clone();
    assert!(
        destination
            .merge_peaks(
                &MSChromatogram::from_peaks(vec![peak(1.0, f32::MAX)]),
                false
            )
            .is_err()
    );
    assert_eq!(destination, before);

    for limits in [
        ChromatogramMergeLimits {
            max_peaks: 3,
            ..Default::default()
        },
        ChromatogramMergeLimits {
            max_bytes: 0,
            ..Default::default()
        },
        ChromatogramMergeLimits {
            max_merged_mzs: 0,
            ..Default::default()
        },
    ] {
        let mut destination = sorted.clone();
        let before = destination.clone();
        assert!(
            destination
                .merge_peaks_with_options(
                    &sorted,
                    ChromatogramMergeOptions {
                        add_meta: true,
                        limits,
                        ..Default::default()
                    },
                )
                .is_err()
        );
        assert_eq!(destination, before);
    }

    // A non-finite product m/z cannot enter the metadata list.
    let mut other = sorted.clone();
    other.product.mz = f64::INFINITY;
    let mut destination = sorted.clone();
    let before = destination.clone();
    assert!(destination.merge_peaks(&other, true).is_err());
    assert_eq!(destination, before);
}

// --------------------------------------------------------------------------
// Mobilogram_test.cpp
// --------------------------------------------------------------------------

// Mobilogram_test.cpp:30-41 -- the shared p1..p3 fixture peaks.
const M1: MobilityPeak1D = MobilityPeak1D::new(2.0, 1.0);
const M2: MobilityPeak1D = MobilityPeak1D::new(10.0, 2.0);

// Mobilogram_test.cpp:945-967 -- the 21-peak search fixture.
fn mobilogram_fixture() -> Mobilogram {
    Mobilogram::from_peaks(
        search_fixture()
            .peaks
            .iter()
            .map(|p| mobility_peak(p.rt, p.intensity))
            .collect(),
    )
}

// Mobilogram_test.cpp:718-727 -- six peaks at mobility 1..6.
fn mobilogram_ladder() -> Mobilogram {
    Mobilogram::from_peaks(vec![
        mobility_peak(1.0, 29.0),
        mobility_peak(2.0, 60.0),
        mobility_peak(3.0, 34.0),
        mobility_peak(4.0, 29.0),
        mobility_peak(5.0, 37.0),
        mobility_peak(6.0, 31.0),
    ])
}

// Mobilogram_test.cpp:450-466 equivalent at 547-559 -- tied keys with identities.
fn mobilogram_tied() -> Mobilogram {
    let mut mobilogram = Mobilogram::from_peaks(vec![
        mobility_peak(1.0, 10.0),
        mobility_peak(2.0, 20.0),
        mobility_peak(1.0, 20.0),
        mobility_peak(2.0, 10.0),
    ]);
    mobilogram
        .string_data_arrays
        .push(text_array("s", &["A", "B", "C", "D"]));
    mobilogram
        .integer_data_arrays
        .push(DataArray::new("i", vec![1, 2, 3, 4]));
    mobilogram
}

fn mobilogram_identities(mobilogram: &Mobilogram) -> (Vec<&str>, Vec<i32>) {
    (
        mobilogram.string_data_arrays[0]
            .data
            .iter()
            .map(String::as_str)
            .collect(),
        mobilogram.integer_data_arrays[0].data.clone(),
    )
}

/// Sections `Mobilogram()` (48), `~Mobilogram()` (55), `[EXTRA] Mobilogram()`
/// (61), `getRT() const` (75), `setRT(double)` (82), `getDriftTimeUnit() const`
/// (91), `getDriftTimeUnitAsString() const` (98) and `setDriftTimeUnit` (105).
#[test]
fn mobilogram_source_construction_and_scalar_accessors() {
    let mobilogram = Mobilogram::new();
    assert_eq!(mobilogram, Mobilogram::default());
    drop(mobilogram);

    let mut mobilogram = Mobilogram::new();
    mobilogram.peaks.push(mobility_peak(47.11, 0.0));
    assert_eq!(mobilogram.len(), 1);
    assert_eq!(mobilogram.peaks[0].mobility, 47.11);

    let mut mobilogram = Mobilogram::new();
    assert_eq!(mobilogram.rt, -1.0);
    mobilogram.rt = 0.451;
    assert_eq!(mobilogram.rt, 0.451);
    assert_eq!(Mobilogram::new().drift_time_unit, DriftTimeUnit::None);
    assert_eq!(Mobilogram::new().drift_time_unit_as_str(), "<NONE>");
    mobilogram.drift_time_unit = DriftTimeUnit::Millisecond;
    assert_eq!(mobilogram.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(mobilogram.drift_time_unit_as_str(), "ms");
}

/// Section `virtual void updateRanges()` (118).
#[test]
fn mobilogram_source_update_ranges() {
    let mut mobilogram = Mobilogram::from_peaks(vec![M1, M2, M1]);
    let first = mobilogram.range_manager().unwrap();
    let second = mobilogram.range_manager().unwrap();
    assert_eq!(first, second);
    assert_eq!(second.max_intensity().unwrap(), 2.0);
    assert_eq!(second.min_intensity().unwrap(), 1.0);
    assert_eq!(second.max_mobility().unwrap(), 10.0);
    assert_eq!(second.min_mobility().unwrap(), 2.0);

    mobilogram.clear();
    mobilogram.peaks.push(M1);
    let ranges = mobilogram.range_manager().unwrap();
    assert_eq!(ranges.max_intensity().unwrap(), 1.0);
    assert_eq!(ranges.min_intensity().unwrap(), 1.0);
    assert_eq!(ranges.max_mobility().unwrap(), 2.0);
    assert_eq!(ranges.min_mobility().unwrap(), 2.0);
    assert_eq!(
        mobilogram.ranges().unwrap().mobility,
        Some(NumericRange { min: 2.0, max: 2.0 })
    );
}

/// Sections `Mobilogram(const Mobilogram&)` (148),
/// `Mobilogram(const Mobilogram&&)` (167), `operator=(const Mobilogram&)` (200)
/// and `operator=(const Mobilogram&&)` (227).
#[test]
fn mobilogram_source_copy_move_and_assignment() {
    let mut original = Mobilogram::from_peaks(vec![mobility_peak(47.11, 0.0)]);
    original.rt = 7.0;
    original.drift_time_unit = DriftTimeUnit::Millisecond;
    let copy = original.clone();
    assert_eq!(copy.rt, 7.0);
    assert_eq!(copy.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(copy.len(), 1);
    assert_eq!(copy.peaks[0].mobility, 47.11);

    // Move construction. The C++ section also asserts the constructor is
    // noexcept and the moved-from value is empty; neither has a runtime
    // counterpart in Rust.
    let mut original =
        Mobilogram::from_peaks(vec![mobility_peak(47.11, 0.0), mobility_peak(48.11, 0.0)]);
    original.rt = 9.0;
    original.drift_time_unit = DriftTimeUnit::InverseReducedMobility;
    let reference = original.clone();
    let moved = original;
    assert_eq!(moved, reference);
    assert_eq!(moved.rt, 9.0);
    assert_eq!(moved.drift_time_unit, DriftTimeUnit::InverseReducedMobility);
    assert_eq!(moved.len(), 2);
    assert_eq!(moved.peaks[0].mobility, 47.11);
    assert_eq!(moved.peaks[1].mobility, 48.11);

    // Copy assignment and assignment of a default value.
    let mut source = Mobilogram::from_peaks(vec![mobility_peak(47.11, 0.0)]);
    source.rt = 7.0;
    source.drift_time_unit = DriftTimeUnit::Millisecond;
    let mut target = Mobilogram::new();
    assert_eq!(target.len(), 0);
    target = source.clone();
    assert_eq!(target.rt, 7.0);
    assert_eq!(target.drift_time_unit, DriftTimeUnit::Millisecond);
    assert_eq!(target.len(), 1);
    assert_eq!(target.peaks[0].mobility, 47.11);
    target = Mobilogram::new();
    assert_eq!(target.rt, -1.0);
    assert_eq!(target.drift_time_unit, DriftTimeUnit::None);
    assert_eq!(target.len(), 0);

    // Move assignment.
    let mut target = Mobilogram::new();
    assert_eq!(target.len(), 0);
    target = moved.clone();
    assert_eq!(target, moved);
    assert_eq!(target.rt, 9.0);
    assert_eq!(
        target.drift_time_unit,
        DriftTimeUnit::InverseReducedMobility
    );
    assert_eq!(target.len(), 2);
    assert_eq!(target.peaks[0].mobility, 47.11);
    assert_eq!(target.peaks[1].mobility, 48.11);
    target = Mobilogram::default();
    assert_eq!(target.rt, -1.0);
    assert_eq!(target.drift_time_unit, DriftTimeUnit::None);
    assert_eq!(target.len(), 0);
}

/// Sections `bool operator==(const Mobilogram&) const` (274) and
/// `bool operator!=(const Mobilogram&) const` (301).
#[test]
fn mobilogram_source_equality() {
    let empty = Mobilogram::new();
    assert!(Mobilogram::new().source_equal(&empty));
    assert_eq!(Mobilogram::new(), empty);

    let mut edit = empty.clone();
    edit.peaks.resize(1, MobilityPeak1D::default());
    assert!(!edit.source_equal(&empty));

    let mut edit = empty.clone();
    edit.drift_time_unit = DriftTimeUnit::Millisecond;
    assert!(!empty.source_equal(&edit));

    let mut edit = empty.clone();
    edit.rt = 5.0;
    assert!(!empty.source_equal(&edit));

    let mut edit = empty.clone();
    edit.peaks.push(M1);
    edit.peaks.push(M2);
    edit.range_manager().unwrap();
    edit.clear();
    assert!(empty.source_equal(&edit));
    assert_eq!(empty, edit);
}

/// Sections `void sortByIntensity(bool reverse = false)` (332) and
/// `void sortByPosition()` (402).
#[test]
fn mobilogram_source_sorts() {
    let mobilities = [
        420.130, 412.824, 423.269, 415.287, 413.800, 419.113, 416.293, 418.232, 414.301, 412.321,
    ];
    let intensities = [
        201.0f32, 60.0, 56.0, 37.0, 34.0, 31.0, 31.0, 31.0, 29.0, 29.0,
    ];
    let mut mobilogram = Mobilogram::from_peaks(
        (0..mobilities.len())
            .map(|i| mobility_peak(mobilities[i], intensities[i]))
            .collect(),
    );
    mobilogram.sort_by_intensity(false).unwrap();
    let mut ascending = intensities;
    ascending.sort_by(f32::total_cmp);
    for (peak, expected) in mobilogram.peaks.iter().zip(ascending) {
        assert_eq!(peak.intensity, expected);
    }

    // sortByPosition, over the descending fixture at Mobilogram_test.cpp:404-405.
    let descending = [
        423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321,
    ];
    let paired = [
        56.0f32, 201.0, 31.0, 31.0, 31.0, 37.0, 29.0, 34.0, 60.0, 29.0,
    ];
    let mut mobilogram = Mobilogram::from_peaks(
        (0..descending.len())
            .map(|i| mobility_peak(descending[i], paired[i]))
            .collect(),
    );
    mobilogram.sort_by_position().unwrap();
    for (peak, expected) in mobilogram.peaks.iter().zip(paired.iter().rev()) {
        assert_eq!(peak.intensity, *expected);
    }
}

/// Section `[EXTRA] sorting reorders all data arrays alongside the peaks` (442),
/// the 34-assertion regression section.
#[test]
fn mobilogram_source_sorting_regression() {
    let make = || {
        let mut mobilogram = Mobilogram::from_peaks(vec![
            mobility_peak(3.0, 10.0),
            mobility_peak(1.0, 30.0),
            mobility_peak(2.0, 20.0),
        ]);
        mobilogram
            .float_data_arrays
            .push(DataArray::new("f1", vec![3.5, 1.5, 2.5]));
        mobilogram
            .string_data_arrays
            .push(text_array("s1", &["mb3", "mb1", "mb2"]));
        mobilogram
            .integer_data_arrays
            .push(DataArray::new("i1", vec![3, 1, 2]));
        mobilogram
    };

    let mut mobilogram = make();
    mobilogram.sort_by_position().unwrap();
    assert_eq!(mobilogram.peaks[0].mobility, 1.0);
    assert_eq!(mobilogram.peaks[1].mobility, 2.0);
    assert_eq!(mobilogram.peaks[2].mobility, 3.0);
    assert_eq!(mobilogram.float_data_arrays[0].data, [1.5, 2.5, 3.5]);
    assert_eq!(mobilogram.string_data_arrays[0].data, ["mb1", "mb2", "mb3"]);
    assert_eq!(mobilogram.integer_data_arrays[0].data, [1, 2, 3]);
    assert_eq!(mobilogram.float_data_arrays[0].name, "f1");
    assert_eq!(mobilogram.string_data_arrays[0].name, "s1");
    assert_eq!(mobilogram.integer_data_arrays[0].name, "i1");

    let mut mobilogram = make();
    mobilogram.sort_by_intensity(false).unwrap();
    assert_eq!(mobilogram.peaks[0].intensity, 10.0);
    assert_eq!(mobilogram.peaks[2].intensity, 30.0);
    assert_eq!(mobilogram.string_data_arrays[0].data[0], "mb3");
    assert_eq!(mobilogram.string_data_arrays[0].data[2], "mb1");
    assert_eq!(mobilogram.integer_data_arrays[0].data[0], 3);
    assert_eq!(mobilogram.integer_data_arrays[0].data[2], 1);

    let mut mobilogram = make();
    mobilogram.sort_by_intensity(true).unwrap();
    assert_eq!(mobilogram.peaks[0].intensity, 30.0);
    assert_eq!(mobilogram.peaks[2].intensity, 10.0);
    assert_eq!(mobilogram.string_data_arrays[0].data[0], "mb1");
    assert_eq!(mobilogram.string_data_arrays[0].data[2], "mb3");
    assert_eq!(mobilogram.integer_data_arrays[0].data[0], 1);
    assert_eq!(mobilogram.integer_data_arrays[0].data[2], 3);

    let mut mobilogram = make();
    mobilogram.integer_data_arrays[0].data.push(99);
    let before = mobilogram.clone();
    assert!(mobilogram.sort_by_position().is_err());
    assert_eq!(mobilogram, before);

    let mut mobilogram = mobilogram_tied();
    mobilogram.sort_by_position().unwrap();
    assert_eq!(mobilogram.peaks[0].mobility, 1.0);
    assert_eq!(mobilogram.peaks[1].mobility, 1.0);
    assert_eq!(mobilogram.peaks[2].mobility, 2.0);
    assert_eq!(mobilogram.peaks[3].mobility, 2.0);
    assert_eq!(
        mobilogram_identities(&mobilogram),
        (vec!["A", "C", "B", "D"], vec![1, 3, 2, 4])
    );

    let mut mobilogram = mobilogram_tied();
    mobilogram.sort_by_intensity(false).unwrap();
    assert_eq!(
        mobilogram_identities(&mobilogram),
        (vec!["A", "D", "B", "C"], vec![1, 4, 2, 3])
    );

    let mut mobilogram = mobilogram_tied();
    mobilogram.sort_by_intensity(true).unwrap();
    assert_eq!(
        mobilogram_identities(&mobilogram),
        (vec!["B", "C", "A", "D"], vec![2, 3, 1, 4])
    );
}

/// Sections `bool isSorted() const` (567) and
/// `template<class Predicate> bool isSorted(const Predicate&) const` (591).
#[test]
fn mobilogram_source_is_sorted_and_predicate_is_sorted() {
    let mut mobilogram = Mobilogram::from_peaks(vec![
        mobility_peak(1000.0, 1.0),
        mobility_peak(1001.0, 1.0),
        mobility_peak(1002.0, 1.0),
    ]);
    assert!(mobilogram.is_sorted().unwrap());
    mobilogram.peaks.reverse();
    assert!(!mobilogram.is_sorted().unwrap());

    let mobilities = [
        423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321,
    ];
    let intensities = [
        56.0f32, 201.0, 31.0, 31.0, 31.0, 37.0, 29.0, 34.0, 60.0, 29.0,
    ];
    let mut mobilogram = Mobilogram::from_peaks(
        (0..mobilities.len())
            .map(|i| mobility_peak(mobilities[i], intensities[i]))
            .collect(),
    );
    mobilogram.sort_by_position().unwrap();
    assert!(
        mobilogram
            .is_sorted_by(|m, a, b| m.peaks[a].mobility < m.peaks[b].mobility)
            .unwrap()
    );
    assert!(mobilogram.is_sorted().unwrap());

    mobilogram.sort_by_intensity(false).unwrap();
    assert!(
        mobilogram
            .is_sorted_by(|m, a, b| m.peaks[a].intensity < m.peaks[b].intensity)
            .unwrap()
    );
    assert!(
        !mobilogram
            .is_sorted_by(|m, a, b| m.peaks[a].mobility < m.peaks[b].mobility)
            .unwrap()
    );
    assert!(!mobilogram.is_sorted().unwrap());
}

/// Section `template<class Predicate> void sort(const Predicate& lambda)` (614).
#[test]
fn mobilogram_source_predicate_sort() {
    let mut mobilogram = Mobilogram::from_peaks(vec![
        mobility_peak(1.0, 10.0),
        mobility_peak(2.0, 20.0),
        mobility_peak(3.0, 30.0),
    ]);
    mobilogram
        .integer_data_arrays
        .push(DataArray::new("rank", vec![2, 0, 1]));
    mobilogram
        .string_data_arrays
        .push(text_array("", &["a", "b", "c"]));
    mobilogram
        .sort_by(|m, i, j| m.integer_data_arrays[0].data[i] < m.integer_data_arrays[0].data[j])
        .unwrap();
    assert_eq!(mobilogram.integer_data_arrays[0].data, [0, 1, 2]);
    assert_eq!(mobilogram.peaks[0].mobility, 2.0);
    assert_eq!(mobilogram.peaks[1].mobility, 3.0);
    assert_eq!(mobilogram.peaks[2].mobility, 1.0);
    assert_eq!(mobilogram.string_data_arrays[0].data, ["b", "c", "a"]);
}

/// Section `Mobilogram& select(const std::vector<Size>& indices)` (648).
#[test]
fn mobilogram_source_select() {
    let mut mobilogram = Mobilogram::from_peaks(vec![
        mobility_peak(1.0, 10.0),
        mobility_peak(2.0, 20.0),
        mobility_peak(3.0, 30.0),
    ]);
    mobilogram
        .float_data_arrays
        .push(DataArray::new("f1", vec![1.5, 2.5, 3.5]));
    mobilogram
        .string_data_arrays
        .push(text_array("", &["a", "b", "c"]));
    mobilogram
        .integer_data_arrays
        .push(DataArray::new("", vec![1, 2, 3]));
    mobilogram.select(&[2, 0]).unwrap();
    assert_eq!(mobilogram.len(), 2);
    assert_eq!(mobilogram.peaks[0].mobility, 3.0);
    assert_eq!(mobilogram.peaks[1].mobility, 1.0);
    assert_eq!(mobilogram.float_data_arrays[0].data, [3.5, 1.5]);
    assert_eq!(mobilogram.string_data_arrays[0].data, ["c", "a"]);
    assert_eq!(mobilogram.integer_data_arrays[0].data, [3, 1]);
    assert_eq!(mobilogram.float_data_arrays[0].name, "f1");

    let mut reversed =
        Mobilogram::from_peaks(vec![mobility_peak(1.0, 10.0), mobility_peak(2.0, 20.0)]);
    reversed
        .string_data_arrays
        .push(text_array("", &["a", "b"]));
    reversed.select(&[1, 0]).unwrap();
    assert_eq!(reversed.len(), 2);
    assert_eq!(reversed.string_data_arrays[0].data, ["b", "a"]);

    let mut bad = Mobilogram::from_peaks(vec![mobility_peak(1.0, 10.0), mobility_peak(2.0, 20.0)]);
    bad.integer_data_arrays
        .push(DataArray::new("", vec![1, 2, 3]));
    let before = bad.clone();
    assert!(bad.select(&[0, 1]).is_err());
    assert_eq!(bad, before);
    assert_eq!(bad.len(), 2);
    assert_eq!(bad.peaks[0].mobility, 1.0);
    assert_eq!(bad.peaks[1].mobility, 2.0);
    assert_eq!(bad.integer_data_arrays[0].data.len(), 3);

    let mut out_of_range =
        Mobilogram::from_peaks(vec![mobility_peak(1.0, 10.0), mobility_peak(2.0, 20.0)]);
    let before = out_of_range.clone();
    assert!(out_of_range.select(&[0, 5]).is_err());
    assert_eq!(out_of_range, before);
    assert_eq!(out_of_range.len(), 2);
}

/// Sections `Iterator MBEnd(CoordinateType)` (731) through
/// `ConstIterator MBBegin(CoordinateType) const` (826): eight C++ overloads over
/// the mobility 1..6 ladder.
#[test]
fn mobilogram_source_mb_bounds() {
    let mobilogram = mobilogram_ladder();
    let at = |index: usize| mobilogram.peaks[index].mobility;

    // MBEnd(mb), mutable and const.
    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_end(4.5).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_end(5.0).unwrap()), 6.0);
        assert_eq!(at(mobilogram.mobility_end(5.5).unwrap()), 6.0);
    }
    // MBBegin(mb), mutable and const.
    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_begin(4.5).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_begin(5.0).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_begin(5.5).unwrap()), 6.0);
    }
    // MBBegin(begin, mb, end), mutable and const.
    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_begin_in(4.5, 0..6).unwrap()), 5.0);
        assert_eq!(mobilogram.mobility_begin_in(4.5, 0..0).unwrap(), 0);
    }
    // MBEnd(begin, mb, end), mutable and const.
    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_end_in(4.5, 0..6).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_end_in(5.0, 0..6).unwrap()), 6.0);
        assert_eq!(mobilogram.mobility_end_in(4.5, 0..0).unwrap(), 0);
    }
}

/// Sections `Iterator PosBegin(CoordinateType)` (841) through
/// `ConstIterator PosEnd(ConstIterator, CoordinateType, ConstIterator) const`
/// (931): the documented aliases of MBBegin()/MBEnd().
#[test]
fn mobilogram_source_pos_bounds_are_mb_bounds() {
    let mobilogram = mobilogram_ladder();
    let at = |index: usize| mobilogram.peaks[index].mobility;
    let last = mobilogram.len() - 1;

    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_begin(4.5).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_begin(5.0).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_begin(5.5).unwrap()), 6.0);
    }
    assert_eq!(at(mobilogram.mobility_begin_in(4.5, 0..6).unwrap()), 5.0);
    assert_eq!(at(mobilogram.mobility_begin_in(5.5, 0..6).unwrap()), 6.0);
    assert_eq!(mobilogram.mobility_begin_in(4.5, 0..0).unwrap(), 0);
    assert_eq!(mobilogram.mobility_begin_in(8.0, 0..6).unwrap(), 6);
    assert_eq!(
        at(mobilogram.mobility_begin_in(8.0, 0..6).unwrap() - 1),
        at(last)
    );
    assert_eq!(at(mobilogram.mobility_begin_in(3.5, 0..6).unwrap()), 4.0);

    for _ in 0..2 {
        assert_eq!(at(mobilogram.mobility_end(4.5).unwrap()), 5.0);
        assert_eq!(at(mobilogram.mobility_end(5.0).unwrap()), 6.0);
        assert_eq!(at(mobilogram.mobility_end(5.5).unwrap()), 6.0);
    }
    assert_eq!(at(mobilogram.mobility_end_in(3.5, 0..6).unwrap()), 4.0);
    assert_eq!(at(mobilogram.mobility_end_in(4.0, 0..6).unwrap()), 5.0);
    assert_eq!(mobilogram.mobility_end_in(4.5, 0..0).unwrap(), 0);
    assert_eq!(at(mobilogram.mobility_end_in(4.5, 0..6).unwrap()), 5.0);
    assert_eq!(at(mobilogram.mobility_end_in(5.0, 0..6).unwrap()), 6.0);
}

/// Sections `Size findNearest(CoordinateType)` (971),
/// `findNearest(CoordinateType, CoordinateType)` (993),
/// `findNearest(CoordinateType, CoordinateType, CoordinateType)` (1017) and
/// `findHighestInWindow(...)` (1045).
#[test]
fn mobilogram_source_searches() {
    let mobilogram = mobilogram_fixture();
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
        assert_eq!(mobilogram.find_nearest(query).unwrap(), Some(index));
    }
    assert_eq!(Mobilogram::new().find_nearest(427.3).unwrap(), None);

    for (query, tolerance, index) in [
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
            mobilogram
                .find_nearest_with_tolerance(query, tolerance)
                .unwrap(),
            index
        );
        assert_eq!(
            mobilogram
                .find_nearest_in_window(query, tolerance, tolerance)
                .unwrap(),
            index
        );
        assert_eq!(
            mobilogram
                .find_highest_in_window(query, tolerance, tolerance)
                .unwrap(),
            index
        );
    }
    assert_eq!(
        Mobilogram::new()
            .find_nearest_in_window(427.3, 1.0, 1.0)
            .unwrap(),
        None
    );

    for (query, left, right, index) in [
        (427.3, 0.1, 0.001, Some(11)),
        (427.3, 0.001, 1.01, None),
        (427.3, 0.001, 1.1, Some(12)),
    ] {
        assert_eq!(
            mobilogram
                .find_nearest_in_window(query, left, right)
                .unwrap(),
            index
        );
        assert_eq!(
            mobilogram
                .find_highest_in_window(query, left, right)
                .unwrap(),
            index
        );
    }
    assert_eq!(
        mobilogram.find_highest_in_window(427.3, 9.0, 4.0).unwrap(),
        Some(8)
    );
    assert_eq!(
        mobilogram
            .find_highest_in_window(430.25, 1.9, 1.01)
            .unwrap(),
        Some(13)
    );
    assert_eq!(
        Mobilogram::new()
            .find_highest_in_window(427.3, 1.0, 1.0)
            .unwrap(),
        None
    );
}

/// Sections `ConstIterator getBasePeak() const` (1077), `Iterator getBasePeak()`
/// (1088) and `PeakType::IntensityType calculateTIC() const` (1099).
#[test]
fn mobilogram_source_base_peak_and_tic() {
    let mobilogram = mobilogram_fixture();
    assert_eq!(mobilogram.base_peak().unwrap().unwrap().intensity, 201.0);
    assert_eq!(mobilogram.base_peak_index().unwrap(), Some(8));
    assert_eq!(Mobilogram::new().base_peak().unwrap(), None);

    let mut mutable = mobilogram_fixture();
    let base = mutable.base_peak_mut().unwrap().unwrap();
    base.intensity += 0.0;
    assert_eq!(base.intensity, 201.0);
    assert_eq!(mutable.base_peak_index().unwrap(), Some(8));

    assert_eq!(mobilogram.calculate_tic().unwrap(), 1032.0);
    assert_eq!(Mobilogram::new().calculate_tic().unwrap(), 0.0);
}

/// Section `void clear()` (1108).
#[test]
fn mobilogram_source_clear() {
    let mut edit = Mobilogram::new();
    edit.peaks.resize(1, MobilityPeak1D::default());
    edit.rt = 5.0;
    edit.drift_time_unit = DriftTimeUnit::Millisecond;
    edit.float_data_arrays.resize(3, DataArray::default());
    edit.integer_data_arrays.resize(2, DataArray::default());
    edit.string_data_arrays.resize(1, DataArray::default());

    edit.clear();
    assert_eq!(edit.len(), 0);
    assert_ne!(edit, Mobilogram::new());
    assert!(edit.is_empty());
    assert!(edit.float_data_arrays.is_empty());
    assert!(edit.integer_data_arrays.is_empty());
    assert!(edit.string_data_arrays.is_empty());
    // The retained scalars are why it is not equal to a default mobilogram.
    assert_eq!(edit.rt, 5.0);
    assert_eq!(edit.drift_time_unit, DriftTimeUnit::Millisecond);
}

/// Section `[Mobilogram::RTLess] bool operator()(const Mobilogram&, const
/// Mobilogram&) const` (1130). The comparator maps onto the public `rt` field.
#[test]
fn mobilogram_source_rt_less() {
    let at = |rt: f64| {
        let mut mobilogram = Mobilogram::new();
        mobilogram.rt = rt;
        mobilogram
    };
    let mut list = [at(3.0), at(2.0), at(1.0)];
    list.sort_by(|a, b| a.rt.total_cmp(&b.rt));
    assert_eq!(list[0].rt, 1.0);
    assert_eq!(list[1].rt, 2.0);
    assert_eq!(list[2].rt, 3.0);

    let (first, second) = (at(0.451), at(0.5));
    let rt_less = |a: &Mobilogram, b: &Mobilogram| a.rt < b.rt;
    assert!(rt_less(&first, &second));
    assert!(!rt_less(&second, &first));
    assert!(!rt_less(&second, &second));
}

/// Section `[EXTRA] std::ostream& operator<<(std::ostream&, const Mobilogram&)`
/// (1165).
#[test]
fn mobilogram_source_stream_layout() {
    let mut mobilogram = Mobilogram::from_peaks(mobilogram_fixture().peaks[..11].to_vec());
    mobilogram.rt = 7.0;
    assert_eq!(
        mobilogram.to_string(),
        "-- MOBILOGRAM BEGIN --\n\
         POS: 412.321 INT: 29\n\
         POS: 412.824 INT: 60\n\
         POS: 413.8 INT: 34\n\
         POS: 414.301 INT: 29\n\
         POS: 415.287 INT: 37\n\
         POS: 416.293 INT: 31\n\
         POS: 418.232 INT: 31\n\
         POS: 419.113 INT: 31\n\
         POS: 420.13 INT: 201\n\
         POS: 423.269 INT: 56\n\
         POS: 426.292 INT: 34\n\
         -- MOBILOGRAM END --\n"
    );
}

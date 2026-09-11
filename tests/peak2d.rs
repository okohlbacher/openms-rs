// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::HasUniqueId;
use openms::kernel::{MobilityPeak2D, Peak2D, RichPeak2D};
use openms::metadata::{MetaValue, Unit};
use std::hash::{DefaultHasher, Hash, Hasher};
fn hash(value: &impl Hash) -> u64 {
    let mut h = DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}

#[test]
fn source_peak_construction_copy_and_mutable_position_literals() {
    let zero = Peak2D::default();
    assert_eq!(zero.position, [0., 0.]);
    assert_eq!(zero.intensity, 0.);
    let mut p = Peak2D::from_position([21.21, 22.22], 123.456f32);
    let copy = p;
    assert_eq!(copy, Peak2D::new(21.21, 22.22, 123.456f32));
    p.set_rt(12345.);
    p.set_mz(12345.);
    p.intensity = 17.8;
    p.position[0] = 876.;
    assert_eq!(p.position, [876., 12345.]);
    assert_eq!(p.intensity, 17.8);
    assert_ne!(p, copy);
    let mut im = MobilityPeak2D::from_position([21.21, 22.22], 123.456f32);
    let im_copy = im;
    im.set_mobility(876.);
    im.set_mz(12345.);
    im.intensity = 17.8;
    assert_eq!((im.mobility(), im.mz(), im.intensity), (876., 12345., 17.8));
    im.position = [10., 12.];
    assert_eq!(im.position, [10., 12.]);
    assert_ne!(im, im_copy);
}

#[test]
fn all_source_dimension_constants_names_units_and_checked_indices() {
    assert_eq!((Peak2D::RT, Peak2D::MZ, Peak2D::DIMENSION), (0, 1, 2));
    assert_eq!(
        (
            MobilityPeak2D::IM,
            MobilityPeak2D::MZ,
            MobilityPeak2D::DIMENSION
        ),
        (0, 1, 2)
    );
    for (d, short, full, unit, full_unit) in [
        (0, "RT", "retention time", "sec", "Seconds"),
        (1, "MZ", "mass-to-charge", "Th", "Thomson"),
    ] {
        assert_eq!(Peak2D::short_dimension_name(d).unwrap(), short);
        assert_eq!(Peak2D::full_dimension_name(d).unwrap(), full);
        assert_eq!(Peak2D::short_dimension_unit(d).unwrap(), unit);
        assert_eq!(Peak2D::full_dimension_unit(d).unwrap(), full_unit);
    }
    assert_eq!(
        [
            Peak2D::short_dimension_name_rt(),
            Peak2D::short_dimension_name_mz()
        ],
        ["RT", "MZ"]
    );
    assert_eq!(
        [
            Peak2D::full_dimension_name_rt(),
            Peak2D::full_dimension_name_mz()
        ],
        ["retention time", "mass-to-charge"]
    );
    assert_eq!(
        [
            Peak2D::short_dimension_unit_rt(),
            Peak2D::short_dimension_unit_mz()
        ],
        ["sec", "Th"]
    );
    assert_eq!(
        [
            Peak2D::full_dimension_unit_rt(),
            Peak2D::full_dimension_unit_mz()
        ],
        ["Seconds", "Thomson"]
    );
    for (d, short, full, unit, full_unit) in [
        (0, "IM", "ion mobility", "?", "?"),
        (1, "MZ", "mass-to-charge", "Th", "Thomson"),
    ] {
        assert_eq!(MobilityPeak2D::short_dimension_name(d).unwrap(), short);
        assert_eq!(MobilityPeak2D::full_dimension_name(d).unwrap(), full);
        assert_eq!(MobilityPeak2D::short_dimension_unit(d).unwrap(), unit);
        assert_eq!(MobilityPeak2D::full_dimension_unit(d).unwrap(), full_unit);
    }
    assert_eq!(
        [
            MobilityPeak2D::short_dimension_name_im(),
            MobilityPeak2D::short_dimension_name_mz()
        ],
        ["IM", "MZ"]
    );
    assert_eq!(
        [
            MobilityPeak2D::full_dimension_name_im(),
            MobilityPeak2D::full_dimension_name_mz()
        ],
        ["ion mobility", "mass-to-charge"]
    );
    assert_eq!(
        [
            MobilityPeak2D::short_dimension_unit_im(),
            MobilityPeak2D::short_dimension_unit_mz()
        ],
        ["?", "Th"]
    );
    assert_eq!(
        [
            MobilityPeak2D::full_dimension_unit_im(),
            MobilityPeak2D::full_dimension_unit_mz()
        ],
        ["?", "Thomson"]
    );
    for d in [2, usize::MAX] {
        assert!(Peak2D::short_dimension_name(d).is_err());
        assert!(Peak2D::full_dimension_name(d).is_err());
        assert!(Peak2D::short_dimension_unit(d).is_err());
        assert!(Peak2D::full_dimension_unit(d).is_err());
        assert!(MobilityPeak2D::short_dimension_name(d).is_err());
        assert!(MobilityPeak2D::full_dimension_name(d).is_err());
        assert!(MobilityPeak2D::short_dimension_unit(d).is_err());
        assert!(MobilityPeak2D::full_dimension_unit(d).is_err());
    }
}

#[test]
fn source_comparator_overloads_map_to_scalar_and_array_comparison() {
    let mut p = [
        Peak2D::new(3., 2.5, 2.5),
        Peak2D::new(2., 3.5, 3.5),
        Peak2D::new(1., 1.5, 1.5),
    ];
    p.sort_by(|a, b| a.intensity.partial_cmp(&b.intensity).unwrap());
    assert_eq!(p.map(|p| p.intensity), [1.5, 2.5, 3.5]);
    p.sort_by(|a, b| a.rt().partial_cmp(&b.rt()).unwrap());
    assert_eq!(p.map(|p| p.position), [[1., 1.5], [2., 3.5], [3., 2.5]]);
    p.sort_by(|a, b| a.position.partial_cmp(&b.position).unwrap());
    assert_eq!(p.map(|p| p.position), [[1., 1.5], [2., 3.5], [3., 2.5]]);
    p.sort_by(|a, b| a.mz().partial_cmp(&b.mz()).unwrap());
    assert_eq!(p.map(|p| p.position), [[1., 1.5], [3., 2.5], [2., 3.5]]);
    let a = Peak2D::new(10., 10., 10.);
    let b = Peak2D::new(12., 12., 12.);
    for (left, right) in [
        (a.rt(), b.rt()),
        (a.mz(), b.mz()),
        (f64::from(a.intensity), f64::from(b.intensity)),
    ] {
        assert!(left < right);
        assert!(right > left);
        assert_eq!(left.partial_cmp(&left), Some(std::cmp::Ordering::Equal));
    }
    assert!(a.position < b.position);
    assert!(a.position < [12., 12.]);
    assert!([10., 10.] < b.position);
    // C++20 std::array ordering and Rust array ordering both preserve unordered NaNs.
    assert!([f64::NAN, 0.].partial_cmp(&[1., 1.]).is_none());
    let mut im = p.map(|p| MobilityPeak2D::from_position(p.position, p.intensity));
    im.sort_by(|a, b| a.position.partial_cmp(&b.position).unwrap());
    assert_eq!(im.map(|p| p.position), [[1., 1.5], [2., 3.5], [3., 2.5]]);
    assert!(MobilityPeak2D::new(1., 2., 3.).mobility() < 2.);
}

#[test]
fn source_display_layout_native_precision_and_equal_zero_hashes() {
    assert_eq!(Peak2D::new(1., 2., 3.).to_string(), "RT: 1 MZ: 2 INT: 3");
    assert_eq!(
        MobilityPeak2D::new(1., 2., 3.).to_string(),
        "IM: 1 MZ: 2 INT: 3"
    );
    assert_eq!(
        format!("{:.2}", Peak2D::new(1., 2., 3.)),
        "RT: 1.00 MZ: 2.00 INT: 3.00"
    );
    for (a, b) in [
        (Peak2D::new(0., 0., 0.), Peak2D::new(-0., -0., -0.)),
        (Peak2D::new(1., 2., 3.), Peak2D::new(1., 2., 3.)),
    ] {
        assert_eq!(a, b);
        assert_eq!(hash(&a), hash(&b));
    }
    assert_eq!(
        hash(&MobilityPeak2D::new(0., 0., 0.)),
        hash(&MobilityPeak2D::new(-0., -0., -0.))
    );
    let p = Peak2D::new(1., 2., 3.);
    for q in [
        Peak2D::new(2., 2., 3.),
        Peak2D::new(1., 3., 3.),
        Peak2D::new(1., 2., 4.),
    ] {
        assert_ne!(p, q);
        assert_ne!(hash(&p), hash(&q));
    }
    let nan = Peak2D::new(f64::NAN, 0., 0.);
    assert_ne!(nan, nan);
    assert!(
        MobilityPeak2D::new(f64::INFINITY, 0., 0.)
            .mobility()
            .is_infinite()
    );
}

#[test]
fn rich_source_copy_metadata_and_plain_assignment_ownership() {
    let mut p = RichPeak2D::new(21.21, 22.22, 123.456f32);
    p.metadata.insert("cluster_id".into(), 4711i64.into());
    p.metadata.insert("text".into(), "bla".into());
    p.unique_id = 17;
    let copy = p.clone();
    p.metadata.insert("text".into(), "bluff".into());
    assert_eq!(copy.metadata["text"].as_str().unwrap(), "bla");
    assert_ne!(p, copy);
    p = copy;
    let ptr = p.metadata["text"].as_str().unwrap().as_ptr();
    let old = p.replace_from_peak(Peak2D::new(1., 2., 3.));
    assert_eq!(old.metadata["text"].as_str().unwrap().as_ptr(), ptr);
    assert_eq!(old.unique_id, 17);
    assert_eq!(old.intensity, 123.456f32);
    assert_eq!(p.peak, Peak2D::new(1., 2., 3.));
    assert!(p.metadata.is_empty());
    assert_eq!(p.unique_id, 0);
    assert_eq!(p.to_string(), "RT: 1 MZ: 2 INT: 3");
    let from_ref = RichPeak2D::from(&p.peak);
    assert_eq!(from_ref, p);
    let mut other = old;
    other.unique_id = 9;
    p.unique_id = 8;
    p.swap_unique_id(&mut other);
    assert_eq!((p.unique_id, other.unique_id), (9, 8));
    assert!(p.metadata.is_empty());
    assert!(p.has_valid_unique_id());
    assert!(!p.has_invalid_unique_id());
    assert_eq!(p.clear_unique_id(), 1);
    assert_eq!(p.clear_unique_id(), 0);
    assert!(p.has_invalid_unique_id());
}

#[test]
fn rich_identity_and_hash_include_metadata_units_and_id() {
    let mut a = RichPeak2D::new(0., 0., 0.);
    let mut b = RichPeak2D::new(-0., -0., -0.);
    a.metadata
        .insert("v".into(), MetaValue::try_from(0.).unwrap());
    b.metadata
        .insert("v".into(), MetaValue::try_from(-0.).unwrap());
    assert_eq!(a, b);
    assert_eq!(hash(&a), hash(&b));
    b.unique_id = 1;
    assert_ne!(a, b);
    assert_ne!(hash(&a), hash(&b));
    b.unique_id = 0;
    b.metadata.insert(
        "v".into(),
        MetaValue::try_from(0.)
            .unwrap()
            .with_unit(Unit::new("UO:0000010", "second", "UO").unwrap())
            .unwrap(),
    );
    assert_ne!(a, b);
    assert_ne!(hash(&a), hash(&b));
}

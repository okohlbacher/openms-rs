// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::metadata::{ContactPerson, Gradient, HPLC};

#[test]
fn source_contact_name_literals_and_owned_metadata() {
    let mut contact = ContactPerson::default();
    for (name, first, last) in [
        ("Diddl Maus", "Diddl", "Maus"),
        ("Normal, Otto", "Otto", "Normal"),
        ("Meiser, Hans F.", "Hans F.", "Meiser"),
    ] {
        contact.set_name(name).unwrap();
        assert_eq!(contact.first_name, first);
        assert_eq!(contact.last_name, last);
    }
    contact.email = "address@example.org".into();
    contact.institution = "Institute".into();
    contact.contact_info = "Information".into();
    contact.url = "https://example.org".into();
    contact.address = "Street".into();
    contact
        .metadata
        .insert("label".into(), "owned value".into());
    let mut copy = contact.clone();
    assert_eq!(copy, contact);
    copy.metadata.clear();
    assert_ne!(copy, contact);
    assert_eq!(contact.metadata.len(), 1);
    copy = contact.clone();
    copy.first_name.clear();
    assert_ne!(copy, contact);
}

#[test]
fn source_contact_delimiter_precedence_empty_fields_and_preserved_first_name() {
    let mut contact = ContactPerson {
        first_name: "retained".into(),
        ..Default::default()
    };
    for (name, first, last) in [
        ("Only", "retained", "Only"),
        ("", "retained", ""),
        (" A B", "", "A"),
        ("A  B", "A", ""),
        ("A B C", "A", "B"),
        ("last, first, ignored", "first", "last"),
        (" \tL\r,\nF \t", "F", "L"),
        (",", "", ""),
        ("F\tL", "", "F\tL"),
        (
            "\u{a0}L\u{a0},\u{a0}F\u{a0}",
            "\u{a0}F\u{a0}",
            "\u{a0}L\u{a0}",
        ),
    ] {
        contact.set_name(name).unwrap();
        assert_eq!(
            (&*contact.first_name, &*contact.last_name),
            (first, last),
            "{name:?}"
        );
    }
    let before = contact.clone();
    assert!(contact.set_name(&"x".repeat(4 * 1024 * 1024 + 1)).is_err());
    assert_eq!(contact, before);
}

fn source_gradient() -> Gradient {
    let mut gradient = Gradient::new();
    for time in [5, 7] {
        gradient.add_timepoint(time).unwrap();
    }
    for eluent in ["A", "B", "C"] {
        gradient.add_eluent(eluent).unwrap();
    }
    for (eluent, values) in [("A", [90, 30]), ("B", [10, 50]), ("C", [0, 20])] {
        for (time, value) in [5, 7].into_iter().zip(values) {
            gradient.set_percentage(eluent, time, value).unwrap();
        }
    }
    gradient
}

#[test]
fn source_gradient_table_and_column_validation_literals() {
    let mut gradient = Gradient::new();
    assert!(gradient.is_valid().unwrap());
    gradient.add_timepoint(5).unwrap();
    assert!(!gradient.is_valid().unwrap());
    gradient = source_gradient();
    assert_eq!(gradient.eluents(), ["A", "B", "C"]);
    assert_eq!(gradient.timepoints(), [5, 7]);
    assert_eq!(
        gradient.percentages(),
        [vec![90, 30], vec![10, 50], vec![0, 20]]
    );
    assert_eq!(gradient.percentage("B", 7).unwrap(), 50);
    assert!(gradient.is_valid().unwrap());
    gradient.set_percentage("A", 5, 91).unwrap();
    assert!(!gradient.is_valid().unwrap());
    gradient.set_percentage("B", 5, 9).unwrap();
    assert!(gradient.is_valid().unwrap());
    gradient.clear_percentages().unwrap();
    assert_eq!(gradient.percentages(), [vec![0, 0], vec![0, 0], vec![0, 0]]);
    assert!(!gradient.is_valid().unwrap());
}

#[test]
fn gradient_signed_order_exact_names_and_rejected_updates_are_atomic() {
    let mut gradient = Gradient::new();
    for time in [i32::MIN, -1, 0, i32::MAX] {
        gradient.add_timepoint(time).unwrap();
    }
    for eluent in ["", "A", "a", "α"] {
        gradient.add_eluent(eluent).unwrap();
    }
    let before = gradient.clone();
    assert!(gradient.add_timepoint(i32::MAX).is_err());
    assert!(gradient.add_timepoint(-2).is_err());
    assert!(gradient.add_eluent("α").is_err());
    assert!(gradient.set_percentage("missing", 0, 20).is_err());
    assert!(gradient.set_percentage("A", 4, 20).is_err());
    assert!(gradient.set_percentage("A", 0, 101).is_err());
    assert!(gradient.percentage("missing", 0).is_err());
    assert!(gradient.percentage("A", 4).is_err());
    assert_eq!(gradient, before);
    gradient.set_percentage("", i32::MIN, 100).unwrap();
    assert_ne!(gradient, before);
    assert_eq!(gradient.percentage("", i32::MIN).unwrap(), 100);
}

#[test]
fn source_clear_axes_preserve_stale_percentage_values_and_explicit_repair() {
    let mut gradient = source_gradient();
    gradient.clear_timepoints();
    assert!(gradient.timepoints().is_empty());
    assert_eq!(gradient.percentages()[0], [90, 30]);
    assert!(gradient.is_valid().unwrap());
    gradient.add_timepoint(-10).unwrap();
    assert_eq!(gradient.percentages()[0], [90, 30, 0]);
    assert_eq!(gradient.percentage("A", -10).unwrap(), 90);
    gradient.clear_eluents();
    assert!(gradient.eluents().is_empty());
    assert_eq!(gradient.percentages().len(), 3);
    gradient.add_eluent("replacement").unwrap();
    assert_eq!(gradient.percentages().len(), 4);
    assert_eq!(gradient.percentage("replacement", -10).unwrap(), 90);
    gradient.clear_percentages().unwrap();
    assert_eq!(gradient.percentages(), [vec![0]]);
    assert_eq!(gradient.percentage("replacement", -10).unwrap(), 0);
}

#[test]
fn stale_source_shape_that_would_index_out_of_bounds_is_checked() {
    let mut gradient = Gradient::new();
    gradient.add_eluent("old empty row").unwrap();
    gradient.clear_eluents();
    gradient.add_timepoint(1).unwrap();
    gradient.add_eluent("new nonempty row").unwrap();
    assert_eq!(gradient.percentages(), [vec![], vec![0]]);
    let before = gradient.clone();
    assert!(gradient.percentage("new nonempty row", 1).is_err());
    assert!(gradient.set_percentage("new nonempty row", 1, 100).is_err());
    assert!(gradient.is_valid().is_err());
    assert_eq!(gradient, before);
    gradient.clear_percentages().unwrap();
    gradient.set_percentage("new nonempty row", 1, 100).unwrap();
    assert!(gradient.is_valid().unwrap());
}

#[test]
fn every_small_two_eluent_composition_matches_independent_total_rule() {
    let mut gradient = Gradient::new();
    gradient.add_timepoint(-3).unwrap();
    gradient.add_eluent("water").unwrap();
    gradient.add_eluent("organic").unwrap();
    for left in 0..=100 {
        for right in 0..=100 {
            gradient.set_percentage("water", -3, left).unwrap();
            gradient.set_percentage("organic", -3, right).unwrap();
            assert_eq!(gradient.is_valid().unwrap(), left + right == 100);
        }
    }
}

#[test]
fn gradient_resource_errors_precede_value_changes() {
    let mut gradient = source_gradient();
    let before = gradient.clone();
    assert!(
        gradient
            .add_eluent(&"x".repeat(Gradient::MAX_BYTES / 4))
            .is_err()
    );
    assert_eq!(gradient, before);
}

#[test]
fn complete_hplc_defaults_fields_copy_and_unsigned_extremes() {
    let mut hplc = HPLC::default();
    assert_eq!((hplc.temperature, hplc.pressure, hplc.flux), (21, 0, 0));
    assert!(hplc.instrument.is_empty() && hplc.column.is_empty() && hplc.comment.is_empty());
    assert!(hplc.gradient.is_valid().unwrap());
    hplc.instrument = "instrument".into();
    hplc.column = "column".into();
    hplc.temperature = i32::MIN;
    hplc.pressure = u32::MAX;
    hplc.flux = u32::MAX;
    hplc.comment = "comment".into();
    hplc.gradient = source_gradient();
    let mut copy = hplc.clone();
    assert_eq!(copy, hplc);
    copy.gradient.clear_timepoints();
    assert_ne!(copy, hplc);
    assert_eq!(hplc.gradient.timepoints(), [5, 7]);
}

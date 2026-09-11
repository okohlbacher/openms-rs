use openms::metadata::{CVTerm, CVTermList, MetaInfo, MetaValue, MetaValueData, Product, Unit};
use std::collections::{HashSet, hash_map::DefaultHasher};
use std::hash::{Hash, Hasher};

fn digest(value: &impl Hash) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

// Observe the bytes supplied to a Hasher, without pinning an implementation's
// numeric digest. These tests ensure individual fields are actually included.
#[derive(Default)]
struct HashInput(Vec<u8>);
impl Hasher for HashInput {
    fn finish(&self) -> u64 {
        0
    }
    fn write(&mut self, bytes: &[u8]) {
        self.0.extend_from_slice(bytes);
    }
}
fn input(value: &impl Hash) -> Vec<u8> {
    let mut hasher = HashInput::default();
    value.hash(&mut hasher);
    hasher.0
}

fn unit() -> Unit {
    Unit::new("UO:0000010", "second", "UO").unwrap()
}
fn value(data: MetaValueData) -> MetaValue {
    MetaValue::new(data).unwrap()
}
fn term(name: &str, data: MetaValueData) -> CVTerm {
    CVTerm {
        accession: "MS:1000016".into(),
        name: name.into(),
        cv_ref: "MS".into(),
        value: value(data),
    }
}

#[test]
fn source_product_hash_shape_and_each_scalar_field() {
    // Product_test.cpp:150: identical 100.5/10/20 products hash equally; a
    // product at 200.5 differs. No source/platform digest is asserted.
    let product = Product {
        mz: 100.5,
        isolation_window_lower_offset: 10.0,
        isolation_window_upper_offset: 20.0,
        ..Default::default()
    };
    assert_eq!(digest(&product), digest(&product.clone()));
    for changed in [
        Product {
            mz: 200.5,
            ..product.clone()
        },
        Product {
            isolation_window_lower_offset: 11.0,
            ..product.clone()
        },
        Product {
            isolation_window_upper_offset: 21.0,
            ..product.clone()
        },
    ] {
        assert_ne!(product, changed);
        assert_ne!(input(&product), input(&changed));
    }
}

#[test]
fn scalar_list_and_nested_product_zeros_obey_equality() {
    for (a, b) in [
        (MetaValueData::Float(0.0), MetaValueData::Float(-0.0)),
        (
            MetaValueData::FloatList(vec![0.0, -0.0]),
            MetaValueData::FloatList(vec![-0.0, 0.0]),
        ),
    ] {
        assert_eq!(a, b);
        assert_eq!(input(&a), input(&b));
        let a = value(a).with_unit(unit()).unwrap();
        let b = value(b).with_unit(unit()).unwrap();
        assert_eq!(a, b);
        assert_eq!(digest(&a), digest(&b));
    }
    let mut base = Product::default();
    base.cv_terms
        .metadata
        .insert("zero".into(), value(MetaValueData::Float(0.0)));
    base.cv_terms
        .add(term("scan time", MetaValueData::FloatList(vec![0.0, 0.0])))
        .unwrap();
    for mz in [0.0, -0.0] {
        for lower in [0.0, -0.0] {
            for upper in [0.0, -0.0] {
                let mut other = Product {
                    mz,
                    isolation_window_lower_offset: lower,
                    isolation_window_upper_offset: upper,
                    cv_terms: base.cv_terms.clone(),
                };
                other
                    .cv_terms
                    .metadata
                    .insert("zero".into(), value(MetaValueData::Float(-0.0)));
                other
                    .cv_terms
                    .replace(term(
                        "scan time",
                        MetaValueData::FloatList(vec![-0.0, -0.0]),
                    ))
                    .unwrap();
                assert_eq!(base, other);
                assert_eq!(input(&base), input(&other));
            }
        }
    }
}

#[test]
fn all_value_types_and_list_boundaries_are_hashed() {
    let values = [
        MetaValueData::Empty,
        MetaValueData::String(String::new()),
        MetaValueData::Integer(0),
        MetaValueData::Float(0.0),
        MetaValueData::StringList(vec![]),
        MetaValueData::IntegerList(vec![]),
        MetaValueData::FloatList(vec![]),
    ];
    let streams: HashSet<_> = values.iter().map(input).collect();
    assert_eq!(streams.len(), values.len());
    for data in values {
        let a = value(data);
        assert_eq!(input(&a), input(&a.clone()));
    }
    for (a, b) in [
        (
            MetaValueData::StringList(vec!["ab".into(), "c".into()]),
            MetaValueData::StringList(vec!["a".into(), "bc".into()]),
        ),
        (
            MetaValueData::IntegerList(vec![1, 2]),
            MetaValueData::IntegerList(vec![2, 1]),
        ),
        (
            MetaValueData::FloatList(vec![1.0, 2.0]),
            MetaValueData::FloatList(vec![2.0, 1.0]),
        ),
        (
            MetaValueData::FloatList(vec![1.0]),
            MetaValueData::FloatList(vec![1.0, 0.0]),
        ),
    ] {
        assert_ne!(a, b);
        assert_ne!(input(&a), input(&b));
    }
}

#[test]
fn hashes_use_typed_content_instead_of_lossy_display_or_epsilon_rounding() {
    let integer: MetaValue = 1_i64.into();
    let float = value(MetaValueData::Float(1.0));
    let string: MetaValue = "1".into();
    assert_eq!(integer.to_string(), float.to_string());
    assert_eq!(float.to_string(), string.to_string());
    assert_ne!(input(&integer), input(&float));
    assert_ne!(input(&float), input(&string));
    let nearby = value(MetaValueData::Float(1.0 + 1e-7));
    assert_ne!(float, nearby); // Existing native equality is exact.
    assert_ne!(input(&float), input(&nearby));
    let one = value(MetaValueData::StringList(vec!["a, b".into()]));
    let two = value(MetaValueData::StringList(vec!["a".into(), "b".into()]));
    assert_eq!(one.to_string(), two.to_string());
    assert_ne!(input(&one), input(&two));
}

#[test]
fn unit_presence_and_all_three_identity_fields_contribute() {
    let base = unit();
    let mut set = HashSet::new();
    set.insert(base.clone());
    set.insert(base.clone());
    assert_eq!(set.len(), 1);
    for other in [
        Unit::new("UO:0000011", "second", "UO").unwrap(),
        Unit::new("UO:0000010", "seconds", "UO").unwrap(),
        Unit::new("UO:0000010", "second", "other").unwrap(),
    ] {
        assert_ne!(base, other);
        assert_ne!(input(&base), input(&other));
        let a = MetaValue::default().with_unit(base.clone()).unwrap();
        let b = MetaValue::default().with_unit(other).unwrap();
        assert_ne!(input(&a), input(&b));
    }
    let no_unit = MetaValue::default();
    let with_unit = no_unit.clone().with_unit(base).unwrap();
    assert!(with_unit.is_empty());
    assert_ne!(no_unit, with_unit);
    assert_ne!(input(&no_unit), input(&with_unit));
}

#[test]
fn every_cv_term_field_changes_hash_input() {
    let base = term("scan time", MetaValueData::Integer(4));
    let mut changes = vec![base.clone(); 5];
    changes[0].accession = "MS:1000017".into();
    changes[1].name = "other name".into();
    changes[2].cv_ref = "other vocabulary".into();
    changes[3].value = 5_i64.into();
    changes[4].value = base.value.clone().with_unit(unit()).unwrap();
    for changed in changes {
        assert_ne!(base, changed);
        assert_ne!(input(&base), input(&changed));
        let mut a = Product::default();
        let mut b = Product::default();
        a.cv_terms.add(base.clone()).unwrap();
        b.cv_terms.add(changed).unwrap();
        assert_ne!(input(&a), input(&b));
    }
}

#[test]
fn cv_bucket_keys_order_duplicates_and_ordinary_metadata_are_retained() {
    let a = term("a", MetaValueData::Integer(1));
    let b = term("b", MetaValueData::Integer(2));
    let mut base = CVTermList::new();
    base.add_terms(&[a.clone(), b.clone()]).unwrap();
    let mut reversed = CVTermList::new();
    reversed.add_terms(&[b, a.clone()]).unwrap();
    assert_ne!(base, reversed);
    assert_ne!(input(&base), input(&reversed));
    let mut duplicate = base.clone();
    duplicate.add(a).unwrap();
    assert_ne!(input(&base), input(&duplicate));
    let mut empty_key = base.clone();
    empty_key.replace_accession("MS:9999999", vec![]).unwrap();
    assert_ne!(input(&base), input(&empty_key));
    let mut metadata = base.clone();
    metadata.metadata.insert("label".into(), "retained".into());
    assert_ne!(input(&base), input(&metadata));
    let mut only_metadata = CVTermList::new();
    only_metadata
        .metadata
        .insert("label".into(), "retained".into());
    assert!(only_metadata.is_empty()); // Predicate ignores ordinary metadata.
    assert_ne!(input(&CVTermList::new()), input(&only_metadata));
}

#[test]
fn equal_maps_hash_equally_regardless_of_insertion_history() {
    let mut a = CVTermList::new();
    let mut b = CVTermList::new();
    let terms = [
        term("time", MetaValueData::Float(2.0)),
        CVTerm::new("MS:1000511", "MS level", "MS"),
    ];
    for term in &terms {
        a.add(term.clone()).unwrap();
    }
    for term in terms.iter().rev() {
        b.add(term.clone()).unwrap();
    }
    let entries = [("first", 1_i64), ("second", 2_i64)];
    for (key, value) in entries {
        a.metadata.insert(key.into(), value.into());
    }
    for (key, value) in entries.into_iter().rev() {
        b.metadata.insert(key.into(), value.into());
    }
    assert_eq!(a, b);
    assert_eq!(input(&a), input(&b));
    assert_eq!(digest(&a), digest(&b));
    let a = Product {
        mz: -5.0,
        cv_terms: a,
        ..Default::default()
    };
    let b = Product {
        mz: -5.0,
        cv_terms: b,
        ..Default::default()
    };
    assert_eq!(digest(&a), digest(&b));
    let mut metadata: MetaInfo = a.cv_terms.metadata.clone();
    metadata.insert("third".into(), MetaValue::default());
    assert_ne!(input(&a.cv_terms.metadata), input(&metadata));
}

#[test]
fn raw_nan_capable_values_stay_partial_and_hash_without_validation() {
    let a = MetaValueData::Float(f64::from_bits(0x7ff8_0000_0000_0001));
    let b = MetaValueData::Float(f64::from_bits(0x7ff8_0000_0000_0002));
    assert_ne!(a, a.clone());
    assert_eq!(input(&a), input(&a.clone()));
    assert_ne!(input(&a), input(&b));
    assert!(MetaValue::new(a).is_err());
    let product = Product {
        mz: f64::NAN,
        isolation_window_lower_offset: -1.0,
        isolation_window_upper_offset: f64::INFINITY,
        ..Default::default()
    };
    assert_ne!(product, product.clone());
    assert!(product.validate().is_err());
    assert_eq!(input(&product), input(&product.clone()));
    let values = MetaValueData::FloatList(vec![f64::NEG_INFINITY, f64::NAN]);
    assert_eq!(input(&values), input(&values.clone()));
}

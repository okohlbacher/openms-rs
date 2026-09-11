use openms::chemistry::{
    AASequence, ProteinProteinCrossLink as Link, ProteinProteinCrossLinkType as Kind,
    TermSpecificity,
};
use std::{
    collections::{HashSet, hash_map::DefaultHasher},
    hash::{Hash, Hasher},
    sync::Arc,
};
fn hash(x: &Link) -> u64 {
    let mut h = DefaultHasher::new();
    x.hash(&mut h);
    h.finish()
}
#[test]
fn source_type_and_default_fields() {
    let mut x = Link::default();
    assert_eq!(x.get_type(), Kind::Loop);
    assert_eq!(x.cross_link_position, (0, 0));
    assert_eq!(x.precursor_correction, 0);
    assert_eq!(x.term_spec_alpha, TermSpecificity::Anywhere);
    x.alpha = Some(Arc::new(AASequence::parse("PEPTIDE").unwrap()));
    x.beta = Some(Arc::new(AASequence::parse("EDEPITPEPE").unwrap()));
    x.cross_link_position = (3, 5);
    assert_eq!(x.get_type(), Kind::Cross);
    x.beta = None;
    assert_eq!(x.get_type(), Kind::Loop);
    x.cross_link_position.1 = -1;
    assert_eq!(x.get_type(), Kind::Mono);
    x.beta = Some(Arc::new(AASequence::default()));
    assert_eq!(x.get_type(), Kind::Mono);
    assert_eq!(Kind::COUNT, 3);
}
#[test]
fn identity_keys_and_every_source_field() {
    let mut x = Link::new(150.).unwrap();
    x.alpha = Some(Arc::new(AASequence::parse("PEPTIDE").unwrap()));
    x.beta = x.alpha.clone();
    x.cross_linker_name = "DSS".into();
    let copy = x.clone();
    assert_eq!(x, copy);
    assert_eq!(hash(&x), hash(&copy));
    let mut set = HashSet::new();
    assert!(set.insert(x.clone()));
    assert!(!set.insert(copy));
    let mut different = x.clone();
    different.alpha = Some(Arc::new((**x.alpha.as_ref().unwrap()).clone()));
    assert_ne!(x, different);
    assert!(set.insert(different));
    for field in 0..7 {
        let mut y = x.clone();
        match field {
            0 => y.beta = None,
            1 => y.cross_link_position.0 = 9,
            2 => y.set_cross_linker_mass(151.).unwrap(),
            3 => y.cross_linker_name.push('x'),
            4 => y.term_spec_alpha = TermSpecificity::NTerm,
            5 => y.term_spec_beta = TermSpecificity::CTerm,
            _ => y.precursor_correction = 1,
        };
        assert_ne!(x, y);
        assert!(set.insert(y));
    }
    let handle = x.alpha.as_ref().unwrap().clone();
    drop(x);
    assert_eq!(handle.as_str(), "PEPTIDE");
}
#[test]
fn finite_mass_and_signed_zero_hash_contract() {
    let mut x = Link::new(-0.).unwrap();
    let y = Link::new(0.).unwrap();
    assert_eq!(x, y);
    assert_eq!(hash(&x), hash(&y));
    x.set_cross_linker_mass(-100.).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(x.set_cross_linker_mass(bad).is_err());
        assert_eq!(x.cross_linker_mass(), -100.);
        assert!(Link::new(bad).is_err());
    }
}

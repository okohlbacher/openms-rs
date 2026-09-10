// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::Result;
use openms::metadata::{MetaInfo, MetaValue, MetaValueData, Unit, UnitOntology};
use openms::param::{DefaultParamHandler as Handler, Param, ParamValue as V};

fn set(p: &mut Param, key: &str, value: V, description: &str) {
    p.set_value(key, value, description, &[]).unwrap();
}
fn source_handler() -> Handler {
    let mut defaults = Param::new();
    set(&mut defaults, "int", V::from(0), "intdesc");
    set(&mut defaults, "string", V::from("default"), "stingdesc");
    let mut handler = Handler::new("dummy").unwrap();
    handler.set_defaults(defaults).unwrap();
    handler.set_subsections(&["ignore".into()]).unwrap();
    assert!(handler.defaults_to_parameters().unwrap().is_empty());
    handler
}

#[test]
fn source_construction_defaults_names_and_partial_updates() {
    let mut empty = Handler::new("dummy2").unwrap();
    assert_eq!(empty.name(), "dummy2");
    assert!(empty.defaults().is_empty());
    assert!(empty.parameters().is_empty());
    assert!(empty.subsections().is_empty());
    assert!(empty.check_defaults());
    assert!(empty.warn_empty_defaults());
    empty.set_name("SetName").unwrap();
    assert_eq!(empty.name(), "SetName");
    let mut handler = source_handler();
    assert_eq!(handler.defaults().size(), 2);
    assert_eq!(handler.parameters().value("int").unwrap(), &V::from(0));
    assert_eq!(
        handler.parameters().value("string").unwrap(),
        &V::from("default")
    );
    let mut supplied = Param::new();
    set(&mut supplied, "int", V::from(1), "");
    set(&mut supplied, "string", V::from("test"), "");
    set(&mut supplied, "ignore:bli", V::from(4711), "");
    assert!(handler.set_parameters(&supplied).unwrap().is_empty());
    assert_eq!(handler.parameters().value("int").unwrap(), &V::from(1));
    assert_eq!(
        handler.parameters().value("string").unwrap(),
        &V::from("test")
    );
    assert_eq!(
        handler.parameters().value("ignore:bli").unwrap(),
        &V::from(4711)
    );
    // Every setParameters starts from the supplied tree, not the previous state.
    let mut partial = Param::new();
    set(&mut partial, "int", V::from(2), "");
    handler.set_parameters(&partial).unwrap();
    assert_eq!(
        handler.parameters().value("string").unwrap(),
        &V::from("default")
    );
    assert!(!handler.parameters().exists("ignore:bli").unwrap());
}

#[test]
fn registered_subsections_skip_validation_but_remain_in_current_parameters() {
    let mut handler = source_handler();
    let mut defaults = handler.defaults().checked_clone().unwrap();
    set(
        &mut defaults,
        "ignore:known",
        V::from(5),
        "ignored delegated option",
    );
    handler.set_defaults(defaults).unwrap();
    let mut input = Param::new();
    set(
        &mut input,
        "ignore:known",
        V::from("wrong type is delegated"),
        "",
    );
    set(&mut input, "ignore:unknown", V::from(0), "");
    set(&mut input, "ignore_other:unknown", V::from(0), "");
    set(&mut input, "unknown", V::from(0), "");
    let warnings = handler.set_parameters(&input).unwrap();
    assert_eq!(warnings.len(), 2);
    assert!(warnings.iter().any(|w| w.contains("ignore_other:unknown")));
    assert!(warnings.iter().any(|w| w.contains("'unknown'")));
    assert_eq!(
        handler.parameters().value("ignore:known").unwrap(),
        &V::from("wrong type is delegated")
    );
}

#[derive(Debug, PartialEq)]
struct Settings {
    count: i32,
    label: String,
}
fn settings(p: &Param) -> Result<Settings> {
    let label = p.value("string")?.as_str()?;
    if label == "callback-error" {
        return Err(openms::Error::InvalidValue("derived invariant".into()));
    }
    Ok(Settings {
        count: p.value("int")?.to_i32()?,
        label: label.into(),
    })
}
#[derive(Debug, PartialEq)]
struct Algorithm {
    handler: Handler,
    settings: Settings,
}
impl Algorithm {
    fn new() -> Self {
        let mut handler = source_handler();
        let (settings, _) = handler.defaults_to_parameters_with(settings).unwrap();
        Self { handler, settings }
    }
    fn update(&mut self, p: &Param) -> Result<Vec<String>> {
        let (settings, warnings) = self.handler.set_parameters_with(p, settings)?;
        self.settings = settings;
        Ok(warnings)
    }
}

#[test]
fn derived_state_lifecycle_is_atomic_on_validation_and_callback_failure() {
    let mut algorithm = Algorithm::new();
    assert_eq!(algorithm.settings.label, "default");
    let mut defaults = algorithm.handler.defaults().checked_clone().unwrap();
    defaults.set_min_int("int", 0).unwrap();
    defaults.set_max_int("int", 10).unwrap();
    algorithm.handler.set_defaults(defaults).unwrap();
    let mut good = Param::new();
    set(&mut good, "int", V::from(7), "");
    set(&mut good, "string", V::from("test"), "");
    algorithm.update(&good).unwrap();
    let before = Algorithm {
        handler: algorithm.handler.checked_clone().unwrap(),
        settings: Settings {
            count: 7,
            label: "test".into(),
        },
    };
    let mut bad = good.checked_clone().unwrap();
    set(&mut bad, "int", V::from(-1), "");
    let called = std::cell::Cell::new(false);
    assert!(
        algorithm
            .handler
            .set_parameters_with(&bad, |_| {
                called.set(true);
                Ok(())
            })
            .is_err()
    );
    assert!(!called.get());
    assert_eq!(algorithm, before);
    set(&mut bad, "int", V::from("not numeric"), "");
    assert!(algorithm.update(&bad).is_err());
    assert_eq!(algorithm, before);
    set(&mut bad, "int", V::from(7), "");
    set(&mut bad, "string", V::from("callback-error"), "");
    assert!(algorithm.update(&bad).is_err());
    assert_eq!(algorithm, before);
}

#[test]
fn defaults_to_parameters_preserves_values_and_warns_only_first_missing_description() {
    let mut handler = source_handler();
    let mut input = Param::new();
    set(&mut input, "string", V::from("keep"), "");
    handler.set_parameters(&input).unwrap();
    let mut defaults = handler.defaults().checked_clone().unwrap();
    set(&mut defaults, "first_missing", V::from(0), "");
    set(&mut defaults, "second_missing", V::from(1), "");
    // Source initialization does not check restrictions on default values.
    defaults.set_min_int("first_missing", 5).unwrap();
    handler.set_defaults(defaults).unwrap();
    let (settings, warnings) = handler.defaults_to_parameters_with(settings).unwrap();
    assert_eq!(settings.label, "keep");
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("'first_missing,'"));
    assert!(!warnings[0].contains("second_missing"));
    let before = handler.checked_clone().unwrap();
    assert!(
        handler
            .defaults_to_parameters_with::<()>(|_| Err(openms::Error::InvalidValue(
                "failed update".into()
            )))
            .is_err()
    );
    assert_eq!(handler, before);
}

#[test]
fn checking_and_empty_default_warnings_are_independent_flags() {
    let mut handler = Handler::new("empty").unwrap();
    let mut input = Param::new();
    set(&mut input, "custom", V::from(1), "");
    assert_eq!(handler.set_parameters(&input).unwrap().len(), 2);
    handler.set_warn_empty_defaults(false);
    assert_eq!(handler.set_parameters(&input).unwrap().len(), 1);
    handler.set_check_defaults(false);
    assert!(handler.set_parameters(&input).unwrap().is_empty());
    let mut defaults = Param::new();
    set(&mut defaults, "custom", V::from("string"), "default");
    set(&mut defaults, "filled", V::from(2), "default");
    handler.set_defaults(defaults).unwrap();
    handler.set_parameters(&input).unwrap();
    assert_eq!(handler.parameters().value("custom").unwrap(), &V::from(1));
    assert_eq!(handler.parameters().value("filled").unwrap(), &V::from(2));
}

#[test]
fn source_equality_and_owned_copy_preserve_policy_fields() {
    let mut a = source_handler();
    let b = a.checked_clone().unwrap();
    assert!(a.source_equal(&b).unwrap());
    assert_eq!(a, b);
    a.set_name("other").unwrap();
    assert!(!a.source_equal(&b).unwrap());
    a = b.checked_clone().unwrap();
    a.set_check_defaults(false);
    assert!(!a.source_equal(&b).unwrap());
    a = b.checked_clone().unwrap();
    a.set_warn_empty_defaults(false);
    assert!(!a.source_equal(&b).unwrap());
    a = b.checked_clone().unwrap();
    a.set_subsections(&["other".into()]).unwrap();
    assert!(!a.source_equal(&b).unwrap());
    a = b.checked_clone().unwrap();
    let mut defaults = a.defaults().checked_clone().unwrap();
    set(&mut defaults, "int", V::from(0), "different description");
    a.set_defaults(defaults).unwrap();
    assert_ne!(a, b);
    assert!(a.source_equal(&b).unwrap());
}

#[test]
fn metadata_uses_source_leaf_keys_prefix_and_last_collision_wins() {
    let mut p = Param::new();
    set(&mut p, "int", V::from(1), "");
    set(&mut p, "string", V::from("test"), "");
    set(&mut p, "ignore:bli", V::from(4711), "");
    set(&mut p, "a:same", V::from(1), "");
    set(&mut p, "b:same", V::from(2), "");
    let mut metadata = MetaInfo::new();
    metadata.insert(
        "keep".into(),
        MetaValue::from("original")
            .with_unit(Unit::from_ontology(UnitOntology::Unit, 10))
            .unwrap(),
    );
    Handler::write_parameters_to_meta_values(&p, &mut metadata, "").unwrap();
    Handler::write_parameters_to_meta_values(&p, &mut metadata, "prefix").unwrap();
    for prefix in ["", "prefix:"] {
        assert_eq!(metadata[&format!("{prefix}int")].as_i64().unwrap(), 1);
        assert_eq!(
            metadata[&format!("{prefix}string")].as_str().unwrap(),
            "test"
        );
        assert_eq!(metadata[&format!("{prefix}bli")].as_i64().unwrap(), 4711);
        assert_eq!(metadata[&format!("{prefix}same")].as_i64().unwrap(), 2);
    }
    assert!(metadata["keep"].unit().is_some());
    let before = metadata.clone();
    Handler::write_parameters_to_meta_values(&p, &mut metadata, "prefix:").unwrap();
    assert_eq!(metadata, before);
    assert!(!metadata.contains_key("ignore:bli"));
}

#[test]
fn all_metadata_value_types_and_nonfinite_conversion_failure_are_atomic() {
    let mut p = Param::new();
    for (key, value) in [
        ("empty", V::Empty),
        ("s", V::from("x\0y")),
        ("i", V::from(i64::MAX)),
        ("f", V::from(1.25)),
        ("sl", V::from(vec![String::from("a")])),
        ("il", V::from(vec![i32::MIN, i32::MAX])),
        ("fl", V::from(vec![1.25, 2.5])),
    ] {
        set(&mut p, key, value, "");
    }
    let mut metadata = MetaInfo::new();
    Handler::write_parameters_to_meta_values(&p, &mut metadata, "").unwrap();
    assert_eq!(metadata["empty"].data(), &MetaValueData::Empty);
    assert_eq!(metadata["s"].as_str().unwrap(), "x\0y");
    assert_eq!(metadata["i"].as_i64().unwrap(), i64::MAX);
    assert_eq!(
        metadata["il"].as_integer_list().unwrap(),
        &[i64::from(i32::MIN), i64::from(i32::MAX)]
    );
    assert_eq!(metadata["fl"].as_float_list().unwrap(), &[1.25, 2.5]);
    let before = metadata.clone();
    set(&mut p, "later", V::from(vec![f64::NAN]), "");
    assert!(Handler::write_parameters_to_meta_values(&p, &mut metadata, "").is_err());
    assert_eq!(metadata, before);
}

#[test]
fn oversized_subsection_update_is_rejected_atomically() {
    let mut handler = source_handler();
    let before = handler.checked_clone().unwrap();
    assert!(
        handler
            .set_subsections(&vec![String::new(); 100_001])
            .is_err()
    );
    assert_eq!(handler, before);
    // Duplicates and a trailing colon are stored verbatim, matching protected
    // source configuration rather than silently normalizing it.
    handler
        .set_subsections(&["ignore".into(), "ignore".into(), "ignore:".into()])
        .unwrap();
    assert_eq!(handler.subsections(), &["ignore", "ignore", "ignore:"]);
}

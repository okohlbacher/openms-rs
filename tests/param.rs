// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::param::{
    CommandLineOptions, MAX_PARAM_DEPTH, MAX_PARAM_ENTRIES, Param, ParamEntry, ParamIterator,
    ParamNode, ParamUpdateOptions, ParamValue as V,
};
fn strings(xs: &[&str]) -> Vec<String> {
    xs.iter().map(|s| (*s).into()).collect()
}
fn set(p: &mut Param, k: &str, v: impl Into<V>) {
    p.set_value(k, v.into(), "", &[]).unwrap();
}
fn text<'a>(p: &'a Param, k: &str) -> &'a str {
    p.value(k).unwrap().as_str().unwrap()
}
fn entry(name: &str, v: impl Into<V>) -> ParamEntry {
    ParamEntry::new(name, v.into(), "", &[]).unwrap()
}

#[test]
fn source_entry_restrictions_and_file_tag_asymmetry() {
    let mut p = Param::new();
    set(&mut p, "int", 5);
    p.set_min_int("int", 5).unwrap();
    p.set_max_int("int", 8).unwrap();
    assert!(p.entry("int").unwrap().is_valid().unwrap());
    set(&mut p, "int", 10);
    assert!(!p.entry("int").unwrap().is_valid().unwrap());
    set(&mut p, "float", 5.1);
    p.set_min_float("float", 5.1).unwrap();
    p.set_max_float("float", 8.1).unwrap();
    assert!(p.entry("float").unwrap().is_valid().unwrap());
    set(&mut p, "float", 10.1);
    assert!(!p.entry("float").unwrap().is_valid().unwrap());
    set(&mut p, "str", "bli");
    p.set_valid_strings("str", &strings(&["bla", "bluff"]))
        .unwrap();
    assert!(!p.entry("str").unwrap().is_valid().unwrap());
    p.add_tag("str", "output prefix").unwrap();
    assert!(p.entry("str").unwrap().is_valid().unwrap());
    set(&mut p, "list", strings(&["bli"]));
    p.set_valid_strings("list", &strings(&["bla"])).unwrap();
    p.add_tag("list", "output prefix").unwrap();
    assert!(!p.entry("list").unwrap().is_valid().unwrap());
    p.add_tag("list", "input file").unwrap();
    assert!(p.entry("list").unwrap().is_valid().unwrap());
    assert!(p.set_min_int("str", 1).is_err());
    assert!(p.valid_strings("int").is_err());
}
#[test]
fn sentinel_bounds_nan_and_checked_integer_narrowing() {
    let mut e = entry("i", i32::MIN);
    assert_eq!(e.min_int, -i32::MAX);
    assert!(e.is_valid().unwrap());
    e.min_int = 0;
    assert!(!e.is_valid().unwrap());
    e.value = V::Integer(i64::MAX);
    assert!(e.is_valid().is_err());
    let mut e = entry("f", f64::INFINITY);
    assert!(e.is_valid().unwrap());
    e.max_float = 4.0;
    assert!(!e.is_valid().unwrap());
    e.value = V::Float(f64::NAN);
    assert!(e.is_valid().unwrap());
    e.min_float = f64::NAN;
    e.max_float = f64::NAN;
    assert!(e.is_valid().unwrap());
    let mut e = entry("ints", vec![1i32, 2, 3]);
    e.max_int = 2;
    assert!(!e.is_valid().unwrap());
    let mut e = entry("floats", vec![1.0, 2.0, 3.0]);
    e.min_float = 2.0;
    assert!(!e.is_valid().unwrap());
}
#[test]
fn source_equality_is_explicit_and_retains_duplicate_quirk() {
    let a = entry("name", 1);
    let mut b = a.clone();
    b.description = "changed".into();
    b.max_int = 1;
    b.tags.insert("advanced".into());
    assert_ne!(a, b);
    assert!(a.source_equal(&b).unwrap());
    let mut left = ParamNode::new("n", "left").unwrap();
    left.entries = vec![entry("a", 1), entry("b", 2)];
    let mut right = ParamNode::new("n", "right").unwrap();
    right.entries = vec![entry("b", 2), entry("a", 1)];
    assert_ne!(left, right);
    assert!(left.source_equal(&right).unwrap());
    left.entries = vec![entry("a", 1), entry("a", 1)];
    assert!(left.source_equal(&right).unwrap());
    assert!(!right.source_equal(&left).unwrap());
}
#[test]
fn source_node_queries_paths_and_insertion_order() {
    let mut root = ParamNode::new("A", "").unwrap();
    root.entries.push(entry("B", 1));
    let mut c = ParamNode::new("C", "").unwrap();
    c.entries = vec![entry("D", 2), entry("E", 3)];
    root.nodes.push(c);
    let mut b = ParamNode::new("B", "").unwrap();
    b.entries.push(entry("G", 4));
    root.nodes.push(b);
    assert_eq!(root.size().unwrap(), 4);
    assert_eq!(root.find_parent_of("C:D").unwrap().unwrap().name, "C");
    assert_eq!(
        root.find_entry_recursive("B:G").unwrap().unwrap().value,
        V::Integer(4)
    );
    assert!(root.find_parent_of("H:C:").unwrap().is_none());
    let mut n = ParamNode::default();
    n.entries.push(entry("H", 5));
    root.insert_node(&n, "F:Z:").unwrap();
    assert_eq!(
        root.find_entry_recursive("F:Z::H").unwrap().unwrap().value,
        V::Integer(5)
    );
    n.name = "W".into();
    root.insert_node(&n, "Q").unwrap();
    assert!(root.find_node("QW").unwrap().is_some());
    root.insert_entry(&entry("H", 7), "FD:ZD:D").unwrap();
    assert!(root.find_entry_recursive("FD:ZD:DH").unwrap().is_some());
    assert_eq!(ParamNode::suffix(":AB"), "AB");
}
#[test]
fn replacement_keeps_restrictions_description_and_replaces_tags() {
    let mut p = Param::new();
    p.set_value("a:x", 1.into(), "old", &strings(&["old"]))
        .unwrap();
    p.set_min_int("a:x", 3).unwrap();
    p.set_value("a:x", 2.into(), "", &strings(&["new"]))
        .unwrap();
    let e = p.entry("a:x").unwrap();
    assert_eq!(e.description, "old");
    assert_eq!(e.min_int, 3);
    assert!(e.tags.contains("new"));
    assert!(!e.tags.contains("old"));
    let mut e = entry("x", 4);
    e.min_int = 100;
    e.description = "newdesc".into();
    p.insert_entry(e, "a:").unwrap();
    assert_eq!(p.entry("a:x").unwrap().min_int, 3);
    assert_eq!(p.description("a:x").unwrap(), "newdesc");
    p.add_section("a", "oldsection").unwrap();
    p.add_section("a", "").unwrap();
    assert_eq!(p.section_description("a").unwrap(), "oldsection");
    let before = p.clone();
    assert!(p.add_section("a:x", "collision").is_err());
    assert_eq!(p, before);
    // Source permits an intermediate node and a leaf to have the same name.
    set(&mut p, "a:x:y", 8);
    assert_eq!(p.value("a:x").unwrap(), &V::Integer(4));
    assert_eq!(p.value("a:x:y").unwrap(), &V::Integer(8));
}
#[test]
fn source_section_prefix_and_trailing_colon_queries() {
    let mut p = Param::new();
    set(&mut p, "test:test:x", 1);
    p.set_section_description("test", "desc").unwrap();
    for key in [
        "test",
        "test:",
        "test:test",
        "test:test:",
        "tes",
        "test:test:x",
    ] {
        assert!(p.has_section(key).unwrap());
    }
    assert!(!p.has_section("missing").unwrap());
    assert!(p.has_section("").is_err());
    assert_eq!(p.section_description("test:").unwrap(), "");
    assert_eq!(p.section_description("test").unwrap(), "desc");
    assert!(p.set_section_description("test:", "no").is_err());
    p.add_section("empty", "emptydesc").unwrap();
    assert!(p.has_section("empty").unwrap());
    assert!(p.copy("empty:", false).unwrap().is_empty());
    assert_eq!(p.copy("empty:", false).unwrap().root().nodes.len(), 0);
}
#[test]
fn source_iterator_traces_empty_sections_and_fused_end() {
    let mut p = Param::new();
    set(&mut p, "A", "1");
    set(&mut p, "r:s:B", "2");
    set(&mut p, "r:s:C", "3");
    p.set_section_description("r:s", "s_desc").unwrap();
    p.add_section("r:", "empty").unwrap();
    set(&mut p, "t:D", "4");
    let mut it = p.iter().unwrap();
    let a = it.next().unwrap();
    assert_eq!(a.key, "A");
    assert!(a.trace.is_empty());
    let b = it.next().unwrap();
    assert_eq!(b.key, "r:s:B");
    assert_eq!(
        b.trace
            .iter()
            .map(|t| (t.name.as_str(), t.opened))
            .collect::<Vec<_>>(),
        vec![("r", true), ("s", true)]
    );
    assert_eq!(b.trace[1].description, "s_desc");
    assert!(it.next().unwrap().trace.is_empty());
    let d = it.next().unwrap();
    assert_eq!(d.key, "t:D");
    assert_eq!(
        d.trace
            .iter()
            .map(|t| (t.name.as_str(), t.opened))
            .collect::<Vec<_>>(),
        vec![("s", false), ("r", false), ("t", true)]
    );
    assert!(it.next().is_none());
    assert!(it.next().is_none());
    assert_eq!(it.end_trace().len(), 1);
    assert_eq!(it.end_trace()[0].name, "t");
    assert!(!it.end_trace()[0].opened);
    let empty = ParamNode::default();
    assert_eq!(ParamIterator::new(&empty).unwrap().count(), 0);
    assert!(p.to_text().unwrap().contains("\"r:s|B\" -> \"2\""));
}
#[test]
fn find_leaf_excludes_root_and_rejects_foreign_cursor() {
    let mut p = Param::new();
    for key in [
        "leaf",
        "a:b:leaf",
        "b:a:leaf",
        "a:c:leaf",
        "a:c:another-leaf",
    ] {
        set(&mut p, key, key);
    }
    let first = p.find_first("leaf").unwrap().unwrap();
    assert_eq!(first.key, "a:b:leaf");
    let next = p.find_next("leaf", first.entry).unwrap().unwrap();
    assert_eq!(next.key, "a:c:leaf");
    assert_eq!(
        p.find_next("leaf", next.entry).unwrap().unwrap().key,
        "b:a:leaf"
    );
    let other = p.clone();
    assert!(
        p.find_next("leaf", other.entry("a:b:leaf").unwrap())
            .is_err()
    );
}
fn source_tree() -> Param {
    let mut p = Param::new();
    for (key, v) in [("test:float", 17.4), ("test2:float", 17.5)] {
        set(&mut p, key, v);
    }
    set(&mut p, "test:string", "test,test,test");
    set(&mut p, "test:int", 17);
    set(&mut p, "test2:string", "test2");
    set(&mut p, "test2:int", 18);
    p.set_section_description("test", "sectiondesc").unwrap();
    p
}
#[test]
fn source_copy_remove_prefix_and_subsets() {
    let p = source_tree();
    assert_eq!(p.copy("test:", false).unwrap().size(), 3);
    let stripped = p.copy("test:", true).unwrap();
    assert_eq!(stripped.value("int").unwrap(), &V::Integer(17));
    assert_eq!(p.copy("test", false).unwrap().size(), 6);
    let stripped = p.copy("test", true).unwrap();
    assert_eq!(stripped.value(":int").unwrap(), &V::Integer(17));
    assert_eq!(stripped.value("2:int").unwrap(), &V::Integer(18));
    let mut subset = Param::new();
    set(&mut subset, "test:irrelevant", 0);
    set(&mut subset, "absent", 0);
    let (copied, warnings) = p.copy_subset_with_messages(&subset).unwrap();
    assert_eq!(copied.size(), 3);
    assert_eq!(warnings.len(), 1);
    let mut p = p.clone();
    set(&mut p, "test:string2", "x");
    p.remove("test").unwrap();
    p.remove("test:strin").unwrap();
    assert_eq!(p.size(), 7);
    p.remove("test:string").unwrap();
    p.remove("test:string2").unwrap();
    p.remove("test:float").unwrap();
    p.remove("test:int").unwrap();
    assert_eq!(p.size(), 3);
    assert!(p.has_section("test").unwrap()); // test2 also starts with test: source prefix predicate!
}
#[test]
fn removal_prunes_only_affected_paths_and_empty_root_terminates() {
    let mut p = source_tree();
    p.add_section("kept", "empty").unwrap();
    p.remove_all("test:float").unwrap();
    assert_eq!(p.size(), 5);
    p.remove_all("test:").unwrap();
    assert_eq!(p.size(), 3);
    p.remove_all("test").unwrap();
    assert!(p.is_empty());
    assert_eq!(p.root().nodes.len(), 1);
    p.remove_all("").unwrap();
    assert!(p.root().nodes.is_empty());
    set(&mut p, "", 1);
    p.remove("").unwrap();
    assert!(p.is_empty());
}
#[test]
fn source_defaults_restrictions_and_merge_descriptions() {
    let mut defaults = Param::new();
    set(&mut defaults, "float", 1.0);
    set(&mut defaults, "float2", 2.0);
    set(&mut defaults, "PATH:onlyfordescription", 45.2);
    defaults
        .set_section_description("PATH", "PATHdesc")
        .unwrap();
    defaults.set_min_float("float2", 0.0).unwrap();
    defaults.add_tag("float2", "advanced").unwrap();
    let mut p = Param::new();
    set(&mut p, "float", -2.0);
    set(&mut p, "PATH:float", -1.0);
    let messages = p.set_defaults(&defaults, "", true).unwrap();
    assert_eq!(messages.len(), 2);
    assert_eq!(p.value("float").unwrap(), &V::Float(-2.0));
    assert_eq!(p.entry("float2").unwrap().min_float, 0.0);
    assert_eq!(p.section_description("PATH").unwrap(), "PATHdesc");
    p.set_defaults(&defaults, "PATH", false).unwrap();
    assert_eq!(p.section_description("PATH:PATH").unwrap(), "PATHdesc");
    assert_eq!(p.value("PATH:float").unwrap(), &V::Float(-1.0));
    let mut merge = Param::new();
    set(&mut merge, "PATH:float", 100.0);
    merge
        .set_section_description("PATH", "replacement")
        .unwrap();
    p.merge(&merge).unwrap();
    assert_eq!(p.value("PATH:float").unwrap(), &V::Float(-1.0));
    assert_eq!(p.section_description("PATH").unwrap(), "PATHdesc");
    let mut isolated = Param::new();
    set(&mut isolated, "PATH:float", 1.0);
    isolated.set_section_description("PATH", "old").unwrap();
    isolated.merge(&merge).unwrap();
    assert_eq!(isolated.section_description("PATH").unwrap(), "replacement");
}
#[test]
fn check_defaults_preserves_prefix_lookup_quirk() {
    let mut defaults = Param::new();
    set(&mut defaults, "int", 5);
    defaults.set_min_int("int", 0).unwrap();
    let mut p = Param::new();
    set(&mut p, "int", -1);
    assert!(p.check_defaults("tool", &defaults, "").is_err());
    p.clear();
    set(&mut p, "pref:int", -1);
    assert!(
        p.check_defaults("tool", &defaults, "pref")
            .unwrap()
            .is_empty()
    );
    set(&mut defaults, "pref:int", 5);
    defaults.set_min_int("pref:int", 0).unwrap();
    assert!(p.check_defaults("tool", &defaults, "pref").is_err());
    set(&mut p, "other", "unknown");
    assert!(
        !p.check_defaults("tool", &Param::new(), "")
            .unwrap()
            .is_empty()
    );
}
#[test]
fn update_mapping_protected_values_partial_success_and_invalids() {
    let mut current = Param::new();
    for k in [
        "good",
        "limited",
        "root:type",
        "tool:1:type",
        "tool:version",
        "version",
        "new:leaf",
    ] {
        set(&mut current, k, 1);
    }
    current.set_max_int("limited", 2).unwrap();
    let mut old = Param::new();
    for k in [
        "good",
        "limited",
        "root:type",
        "tool:1:type",
        "tool:version",
        "version",
        "old:leaf",
        "unknown",
    ] {
        set(&mut old, k, 3);
    }
    let report = current
        .update_with_options(
            &old,
            ParamUpdateOptions {
                fail_on_invalid_values: true,
                fail_on_unknown_parameters: true,
                add_unknown: true,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(!report.success);
    for k in ["good", "root:type", "version", "new:leaf"] {
        assert_eq!(current.value(k).unwrap(), &V::Integer(3));
    }
    for k in ["limited", "tool:1:type", "tool:version"] {
        assert_eq!(current.value(k).unwrap(), &V::Integer(1));
    }
    assert!(!current.exists("unknown").unwrap());
    let mut ambiguous = Param::new();
    set(&mut ambiguous, "a:leaf", 1);
    set(&mut ambiguous, "b:leaf", 2);
    set(&mut ambiguous, "leaf", 9);
    let mut old = Param::new();
    set(&mut old, "old:leaf", 4);
    ambiguous.update(&old, true).unwrap();
    assert_eq!(ambiguous.value("old:leaf").unwrap(), &V::Integer(4));
    assert_eq!(ambiguous.value("a:leaf").unwrap(), &V::Integer(1));
}
#[test]
fn update_equal_invalid_value_is_not_revalidated_and_errors_are_atomic() {
    let mut current = Param::new();
    set(&mut current, "same", 9);
    current.set_max_int("same", 0).unwrap();
    let old = current.clone();
    assert!(
        current
            .update_with_options(
                &old,
                ParamUpdateOptions {
                    fail_on_invalid_values: true,
                    ..Default::default()
                }
            )
            .unwrap()
            .success
    );
    let mut old = Param::new();
    set(&mut old, "same", i64::MAX);
    let before = current.clone();
    assert!(current.update(&old, false).is_err());
    assert_eq!(current, before);
}
#[test]
fn literal_source_cli_both_overloads_and_option_boundaries() {
    let mut p = Param::new();
    p.parse_command_line(
        &strings(&[
            "executable",
            "-a",
            "-1.0",
            "-b",
            "bv",
            "-c",
            "cv",
            "rv1",
            "rv2",
            "-1.0",
        ]),
        "test4",
    )
    .unwrap();
    assert_eq!(text(&p, "test4:-a"), "-1.0");
    assert_eq!(
        p.value("test4:misc").unwrap(),
        &V::StringList(strings(&["rv1", "rv2", "-1.0"]))
    );
    let mut opts = CommandLineOptions::default();
    opts.one_argument.insert("-a".into(), "a".into());
    opts.no_argument.insert("-b".into(), "b".into());
    let mut p = Param::new();
    p.parse_command_line_mapped(
        &strings(&["exe", "-a", "av", "-b", "bv", "-c", "cv", "rv1", "rv2"]),
        &opts,
    )
    .unwrap();
    assert_eq!(text(&p, "a"), "av");
    assert_eq!(text(&p, "b"), "true");
    assert_eq!(
        p.value("misc").unwrap(),
        &V::StringList(strings(&["bv", "cv", "rv1", "rv2"]))
    );
    assert_eq!(
        p.value("unknown").unwrap(),
        &V::StringList(strings(&["-c"]))
    );
    for k in ["d", "e", "f", "g"] {
        opts.multiple_arguments.insert(format!("-{k}"), k.into());
    }
    p.clear();
    p.parse_command_line_mapped(
        &strings(&["mult", "-d", "1.333", "2.23", "3", "-e", "4", "-f", "-g"]),
        &opts,
    )
    .unwrap();
    assert_eq!(
        p.value("d").unwrap(),
        &V::StringList(strings(&["1.333", "2.23", "3"]))
    );
    assert_eq!(p.value("f").unwrap(), &V::StringList(vec![]));
    assert_eq!(p.value("g").unwrap(), &V::StringList(vec![]));
    p.clear();
    p.parse_command_line(&strings(&["exe", "-", "-3", "-.3", "--", "-last"]), "")
        .unwrap();
    assert_eq!(
        p.value("misc").unwrap(),
        &V::StringList(strings(&["-", "-3"]))
    );
    assert_eq!(text(&p, "-.3"), "");
    assert_eq!(text(&p, "--"), "");
    assert_eq!(text(&p, "-last"), "");
}
#[test]
fn invalid_tags_late_cli_errors_and_depth_are_atomic() {
    let mut p = Param::new();
    p.set_value("x", 1.into(), "", &strings(&["comma,allowed-by-setvalue"]))
        .unwrap();
    let before = p.clone();
    assert!(p.add_tags("x", &strings(&["valid", "bad,tag"])).is_err());
    assert_eq!(p, before);
    assert!(p.set_valid_strings("x", &strings(&["a"])).is_err());
    set(&mut p, "misc", 1);
    let before = p.clone();
    assert!(
        p.parse_command_line(&strings(&["exe", "-new", "value", "plain"]), "")
            .is_err()
    );
    assert_eq!(p, before);
    let too_deep = "a:".repeat(MAX_PARAM_DEPTH + 1);
    assert!(p.set_value(&too_deep, 1.into(), "", &[]).is_err());
    assert_eq!(p, before);
    let mut node = ParamNode::default();
    for _ in 0..MAX_PARAM_DEPTH + 1 {
        let mut parent = ParamNode::default();
        parent.nodes.push(node);
        node = parent;
    }
    assert!(Param::from_root(node).is_err());
    let mut node = ParamNode::default();
    node.entries = vec![entry("x", 0); MAX_PARAM_ENTRIES + 1];
    assert!(Param::from_root(node).is_err());
}

#[test]
fn source_update_requeries_unresolvable_literal_colon_draft_names() {
    let mut root = ParamNode::default();
    root.entries.push(entry("a:b", 1));
    let old = Param::from_root(root).unwrap();
    let mut current = Param::new();
    let before = current.clone();
    assert!(current.update(&old, true).is_err());
    assert_eq!(current, before);
}

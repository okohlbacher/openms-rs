// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::system::file::{
    self as f, CopyOptions, FileContext, MatchingFileListsStatus as Match, ResolvedDataPath,
    TempDir, TempFile,
};
use std::{
    fs,
    path::{Path, PathBuf},
    sync::{Arc, Barrier},
    time::{Duration, SystemTime, UNIX_EPOCH},
};

fn temp() -> TempDir {
    TempDir::new_in(std::env::temp_dir(), false).unwrap()
}
fn context(base: &Path) -> FileContext {
    let mut c = FileContext::new(base.join("bin"), base.join("home"), base.join("temp")).unwrap();
    fs::create_dir_all(&c.executable_directory).unwrap();
    fs::create_dir_all(&c.temporary_directory).unwrap();
    c.data_candidates.clear();
    c
}
fn data(base: &Path) {
    fs::create_dir_all(base.join("CHEMISTRY")).unwrap();
    fs::write(
        base.join("CHEMISTRY/unimod.xml"),
        "marker only; source checks existence",
    )
    .unwrap();
}

#[test]
fn source_filename_and_compound_extension_literals() {
    for (value, parent, base) in [
        (
            "/source/config/bla/bluff.h",
            "/source/config/bla",
            "bluff.h",
        ),
        (r"c:\config\bla\tuff.h", r"c:\config\bla", "tuff.h"),
        ("filename_only.h", ".", "filename_only.h"),
        ("/path/only/", "/path/only", ""),
        ("/", "", ""),
        ("", ".", ""),
    ] {
        assert_eq!(f::path(value), parent);
        assert_eq!(f::basename(value), base);
    }
    for (value, stem, extension) in [
        ("/path/to/sample.mzML", "sample", ".mzML"),
        ("/path/to/sample.mzML.gz", "sample", ".mzML.gz"),
        ("/path/to/file.txt", "file", ".txt"),
        ("/path/to/file.txt.tgz", "file.txt", ".tgz"),
        ("/path/to/file", "file", ""),
        ("experiment.featureXML", "experiment", ".featureXML"),
        ("", "", ""),
        ("/home.with.dot/filename", "filename", ""),
        (r"c:\data\sample.idXML", "sample", ".idXML"),
        (".mzML", "", ".mzML"),
        ("Δ/test.mzML.GZ", "test", ".mzML.GZ"),
    ] {
        assert_eq!(f::stem_name(value), stem, "{value}");
        assert_eq!(f::extension(value), extension, "{value}");
    }
}
#[test]
fn source_file_queries_exact_fifteen_byte_fixture_and_signed_epoch() {
    let d = temp();
    let text = d.path().join("text");
    let empty = d.path().join("empty");
    let missing = d.path().join("missing");
    fs::write(&text, include_bytes!("data/system_file_text.txt")).unwrap();
    fs::write(&empty, include_bytes!("data/system_file_empty.txt")).unwrap();
    assert_eq!(f::file_size(&text).unwrap(), 15);
    assert_eq!(f::file_size(&empty).unwrap(), 0);
    assert!(f::file_size(&missing).is_err());
    assert!(f::file_size(d.path()).is_err());
    assert!(f::exists(&text));
    assert!(!f::exists(&missing));
    assert!(!f::empty(&text));
    assert!(f::empty(&empty));
    assert!(f::empty(&missing));
    assert!(f::empty(d.path()));
    assert!(f::readable(&text));
    assert!(f::readable(d.path()));
    assert!(!f::readable(&missing));
    assert!(f::is_directory(d.path()));
    assert!(!f::is_directory(&text));
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs() as i64;
    assert!((f::get_modification_time(&text).unwrap() - now).abs() <= 5);
    fs::File::options()
        .write(true)
        .open(&text)
        .unwrap()
        .set_times(fs::FileTimes::new().set_modified(UNIX_EPOCH - Duration::from_millis(500)))
        .unwrap();
    assert_eq!(f::get_modification_time(&text).unwrap(), -1);
    assert!(f::get_modification_time(&missing).is_err());
    assert!(f::absolute_path(".").unwrap().is_absolute());
    assert!(f::get_executable_path().unwrap().is_dir());
}
#[test]
fn native_permissions_probe_never_creates_or_truncates_requested_path() {
    let d = temp();
    let present = d.path().join("present");
    let absent = d.path().join("absent");
    fs::write(&present, b"payload").unwrap();
    assert!(f::writable(&present));
    assert_eq!(fs::read(&present).unwrap(), b"payload");
    assert!(f::writable(&absent));
    assert!(!absent.exists());
    assert!(f::writable(d.path()));
    assert!(!f::writable(d.path().join("absent_parent/file")));
    assert!(!f::writable(""));
    assert_eq!(fs::read_dir(d.path()).unwrap().count(), 1);
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&present, fs::Permissions::from_mode(0o600)).unwrap();
        assert!(!f::executable(&present));
        fs::set_permissions(&present, fs::Permissions::from_mode(0o601)).unwrap();
        assert!(f::executable(&present));
    }
}
#[test]
fn concurrent_writable_probes_preserve_real_writer_output() {
    let d = temp();
    let path = Arc::new(d.path().join("output"));
    for _ in 0..32 {
        let barrier = Arc::new(Barrier::new(3));
        let mut threads = Vec::new();
        for _ in 0..2 {
            let b = barrier.clone();
            let p = path.clone();
            threads.push(std::thread::spawn(move || {
                b.wait();
                assert!(f::writable(&*p));
            }));
        }
        barrier.wait();
        fs::write(&*path, b"important").unwrap();
        for thread in threads {
            thread.join().unwrap();
        }
        assert_eq!(fs::read(&*path).unwrap(), b"important");
        fs::remove_file(&*path).unwrap();
    }
    assert_eq!(fs::read_dir(d.path()).unwrap().count(), 0);
}
#[test]
fn native_copy_and_moves_preserve_destination_on_failures() {
    let d = temp();
    let a = d.path().join("a");
    let b = d.path().join("b");
    fs::write(&a, b"source").unwrap();
    fs::write(&b, b"target").unwrap();
    assert!(f::copy(&a, &b).is_err());
    assert!(f::rename(&a, &b, false).is_err());
    assert_eq!(fs::read(&b).unwrap(), b"target");
    assert_eq!(fs::read(&a).unwrap(), b"source");
    assert!(f::rename(d.path().join("missing"), &b, true).is_err());
    assert_eq!(fs::read(&b).unwrap(), b"target");
    f::rename(&a, &a, false).unwrap();
    f::rename(&a, &b, true).unwrap();
    assert!(!a.exists());
    assert_eq!(fs::read(&b).unwrap(), b"source");
    f::copy(&b, &a).unwrap();
    assert_eq!(fs::read(&a).unwrap(), b"source");
    f::remove(&a).unwrap();
    f::remove(&a).unwrap();
    assert!(!a.exists());
}
#[test]
fn concurrent_no_clobber_moves_publish_exactly_one_complete_winner() {
    let d = temp();
    let destination = Arc::new(d.path().join("winner"));
    let barrier = Arc::new(Barrier::new(12));
    let mut threads = Vec::new();
    for i in 0..12 {
        let from = d.path().join(format!("input-{i}"));
        fs::write(&from, vec![i as u8; 4096]).unwrap();
        let to = destination.clone();
        let go = barrier.clone();
        threads.push(std::thread::spawn(move || {
            go.wait();
            let result = f::rename(&from, &*to, false);
            (i, from, result)
        }));
    }
    let mut winner = None;
    for t in threads {
        let (i, path, result) = t.join().unwrap();
        if result.is_ok() {
            assert!(winner.replace(i).is_none());
            assert!(!path.exists());
        } else {
            assert!(path.exists());
        }
    }
    assert_eq!(
        fs::read(&*destination).unwrap(),
        vec![winner.unwrap() as u8; 4096]
    );
}
#[test]
fn recursive_copy_modes_preflight_conflicts_and_source_contents() {
    let d = temp();
    let a = d.path().join("source");
    let b = d.path().join("target");
    f::make_dir(a.join("pdata/1")).unwrap();
    fs::write(a.join("pdata/1/proc"), b"original").unwrap();
    f::copy_dir_recursively(&a, &b, CopyOptions::Overwrite).unwrap();
    assert_eq!(fs::read(b.join("pdata/1/proc")).unwrap(), b"original");
    fs::write(b.join("pdata/1/proc"), b"modified").unwrap();
    fs::write(a.join("new"), b"new").unwrap();
    assert!(f::copy_dir_recursively(&a, &b, CopyOptions::Cancel).is_err());
    assert!(!b.join("new").exists());
    f::copy_dir_recursively(&a, &b, CopyOptions::Skip).unwrap();
    assert_eq!(fs::read(b.join("pdata/1/proc")).unwrap(), b"modified");
    assert!(b.join("new").exists());
    f::copy_dir_recursively(&a, &b, CopyOptions::Overwrite).unwrap();
    assert_eq!(fs::read(b.join("pdata/1/proc")).unwrap(), b"original");
    assert!(f::copy_dir_recursively(&a, &a, CopyOptions::Overwrite).is_err());
    assert!(f::copy_dir_recursively(&a, a.join("new-child"), CopyOptions::Overwrite).is_err());
    assert!(!a.join("new-child").exists());
    f::remove_dir(&b).unwrap();
    assert!(!b.exists());
    f::remove_dir_recursively(&b).unwrap();
}
#[test]
#[cfg(unix)]
fn source_directory_links_are_followed_for_copy_but_never_recursive_removal() {
    use std::os::unix::fs::symlink;
    let d = temp();
    let source = d.path().join("source");
    let destination = d.path().join("dest");
    f::make_dir(&source).unwrap();
    symlink(&source, source.join("loop")).unwrap();
    assert!(f::copy_dir_recursively(&source, &destination, CopyOptions::Overwrite).is_err());
    assert!(!destination.exists());
    f::remove(source.join("loop")).unwrap();
    let external = d.path().join("external");
    fs::create_dir(&external).unwrap();
    fs::write(external.join("keep"), b"yes").unwrap();
    symlink(&external, source.join("link")).unwrap();
    f::copy_dir_recursively(&source, &destination, CopyOptions::Overwrite).unwrap();
    assert_eq!(fs::read(destination.join("link/keep")).unwrap(), b"yes");
    f::remove_dir_recursively(&source).unwrap();
    assert!(external.join("keep").exists());
    let dangling = d.path().join("dangling");
    symlink("absent", &dangling).unwrap();
    f::remove(&dangling).unwrap();
    assert!(fs::symlink_metadata(dangling).is_err());
}
#[test]
fn directory_lists_wildcards_and_sorted_results() {
    let d = temp();
    for name in ["subB", "subA"] {
        f::make_dir(d.path().join(name)).unwrap();
    }
    for name in ["a.txt", "b.txt", "c.csv", "Δ.txt", ".hidden", "x[1].txt"] {
        fs::write(d.path().join(name), b"x").unwrap();
    }
    let dirs = f::list_directories(d.path()).unwrap();
    assert_eq!(dirs, vec![d.path().join("subA"), d.path().join("subB")]);
    for (pattern, expected) in [
        ("[ab].txt", vec!["a.txt", "b.txt"]),
        ("[!a-c].txt", vec!["Δ.txt"]),
        (r"x\[1\].txt", vec!["x[1].txt"]),
        ("?.csv", vec!["c.csv"]),
        (".*", vec![".hidden"]),
        ("missing*", vec![]),
    ] {
        assert_eq!(
            f::file_list(d.path(), pattern, false).unwrap(),
            expected.into_iter().map(PathBuf::from).collect::<Vec<_>>()
        );
    }
    assert_eq!(
        f::file_list(d.path(), "a.txt", true).unwrap(),
        vec![d.path().join("a.txt")]
    );
    assert!(f::file_list(d.path().join("absent"), "*", false).is_err());
}
#[test]
fn owned_temporary_guards_keep_close_and_alternative_behavior() {
    let d = temp();
    let child = TempDir::new_in(d.path().join("missing/parents"), false).unwrap();
    let path = child.path().to_owned();
    assert!(path.is_dir());
    drop(child);
    assert!(!path.exists());
    let kept = TempDir::new_in(d.path(), true).unwrap();
    let path = kept.path().to_owned();
    drop(kept);
    assert!(path.is_dir());
    let file = TempFile::new_in(d.path()).unwrap();
    let path = file.path().to_owned();
    fs::write(&path, b"temporary").unwrap();
    drop(file);
    assert!(!path.exists());
    let kept = TempFile::new_in(d.path()).unwrap().keep();
    assert!(kept.exists());
    let unowned = TempFile::alternative(d.path().join("not-created")).unwrap();
    let path = unowned.path().to_owned();
    drop(unowned);
    assert!(!path.exists());
    let unowned = TempFile::alternative(&kept).unwrap();
    unowned.close().unwrap();
    assert!(kept.exists());
    TempFile::new_in(d.path()).unwrap().close().unwrap();
    TempDir::new_in(d.path(), false).unwrap().close().unwrap();
    let a = f::get_unique_name(false).unwrap();
    let b = f::get_unique_name(false).unwrap();
    assert_ne!(a, b);
    assert_eq!(a.split('_').count(), 4);
    assert_eq!(a.split('_').next().unwrap().len(), 8);
}
#[test]
fn source_matching_filenames_including_equal_length_set_multiplicity_quirk() {
    let list = |s: &[&str]| s.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    for (a, b, expected) in [
        (
            vec!["file1.txt", "file2.txt"],
            vec!["file1.txt", "file2.txt"],
            Match::Match,
        ),
        (
            vec!["file1.txt", "file2.txt"],
            vec!["file2.txt", "file1.txt"],
            Match::OrderMismatch,
        ),
        (
            vec!["file1.txt", "file2.txt"],
            vec!["file1.txt", "file3.txt"],
            Match::SetMismatch,
        ),
        (
            vec!["file1.txt", "file2.txt"],
            vec!["file1.txt"],
            Match::SetMismatch,
        ),
        (
            vec!["/a/file1.txt", "/b/file2.mzML"],
            vec!["/x/file1.mzML", "/y/file2.txt"],
            Match::Match,
        ),
        (vec![], vec![], Match::Match),
        (
            vec!["a", "a", "b"],
            vec!["a", "b", "b"],
            Match::OrderMismatch,
        ),
    ] {
        assert_eq!(
            f::validate_matching_file_names(&list(&a), &list(&b), true, true).unwrap(),
            expected
        );
    }
    assert_eq!(
        f::validate_matching_file_names(
            &list(&["/a/file.txt"]),
            &list(&["/b/file.txt"]),
            false,
            false
        )
        .unwrap(),
        Match::SetMismatch
    );
}
#[test]
fn explicit_data_context_precedence_cache_retry_and_resource_search() {
    let d = temp();
    let mut c = context(d.path());
    let one = d.path().join("one");
    let two = d.path().join("two");
    data(&one);
    data(&two);
    c.data_candidates = vec![
        ResolvedDataPath {
            path: one.clone(),
            source: "first configured".into(),
        },
        ResolvedDataPath {
            path: two.clone(),
            source: "second configured".into(),
        },
    ];
    c.data_override = Some(d.path().join("missing"));
    assert!(c.get_openms_data_path().is_err());
    c.data_override = Some(two.clone());
    assert_eq!(c.get_openms_data_path().unwrap(), two);
    assert!(
        c.get_openms_data_path_source()
            .unwrap()
            .contains("override")
    );
    c.data_override = None;
    assert_eq!(c.get_openms_data_path().unwrap(), two);
    c.clear_data_cache();
    assert_eq!(c.get_openms_data_path().unwrap(), one);
    fs::write(one.join("resource"), b"one").unwrap();
    fs::write(two.join("resource"), b"two").unwrap();
    assert_eq!(
        c.find("resource", std::slice::from_ref(&two)).unwrap(),
        two.join("resource")
    );
    assert_eq!(c.find("resource", &[]).unwrap(), one.join("resource"));
    let found = c.find("CHEMISTRY/unimod.xml", &[]).unwrap();
    assert_eq!(c.find(&found, &[]).unwrap(), found);
    let message = c.find("missing", &[]).unwrap_err().to_string();
    assert!(message.contains("OPENMS_DATA_PATH"));
    assert!(message.contains("first configured"));
    assert!(c.find(" \t", &[]).is_err());
    c.clear_data_cache();
    c.data_override = Some(d.path().join("missing"));
    assert!(c.find("resource", &[two]).is_err());
    assert!(c.find(found, &[]).is_ok());
}
#[test]
fn documentation_search_is_explicit_and_needs_data_like_source() {
    let d = temp();
    let mut c = context(d.path());
    let data_path = d.path().join("data");
    data(&data_path);
    c.data_override = Some(data_path);
    let docs = d.path().join("docs");
    fs::create_dir(&docs).unwrap();
    fs::write(docs.join("guide.txt"), b"guide").unwrap();
    c.documentation_directories.push(docs.clone());
    assert_eq!(c.find_doc("guide.txt").unwrap(), docs.join("guide.txt"));
    assert!(c.find_doc("missing").is_err());
}
#[test]
fn source_path_split_and_executable_search_are_not_permission_checks() {
    let separator = if cfg!(windows) { ';' } else { ':' };
    let split =
        f::get_path_locations(&format!("/usr/bin{separator}{separator}relative\\tools")).unwrap();
    assert_eq!(split, vec!["/usr/bin/", "/", "relative/tools/"]);
    assert!(f::get_path_locations("").unwrap().is_empty());
    let d = temp();
    let mut c = context(d.path());
    let other = d.path().join("other");
    fs::create_dir(&other).unwrap();
    let name = if cfg!(windows) { "echo.exe" } else { "echo" };
    fs::write(other.join(name), b"not executable").unwrap();
    c.search_path.push(other.clone());
    assert_eq!(c.find_executable(name).unwrap(), Some(other.join(name)));
    assert!(c.find_sibling_topp_executable(name).is_err());
    assert!(c.find_executable("absent").unwrap().is_none());
    fs::write(c.executable_directory.join(name), b"sibling").unwrap();
    assert_eq!(
        c.find_sibling_topp_executable(name).unwrap(),
        c.executable_directory.join(name)
    );
    assert_eq!(
        c.find_executable(other.join(name)).unwrap(),
        Some(other.join(name))
    );
    assert_eq!(c.get_executable_path(), c.executable_directory);
}
#[test]
fn missing_system_configuration_returns_five_defaults_without_creating_files() {
    let d = temp();
    let mut c = context(d.path());
    let p = c.get_system_parameters().unwrap();
    assert_eq!(p.size(), 5);
    assert_eq!(
        p.value("version").unwrap().as_str().unwrap(),
        openms::CORE_SDK_VERSION
    );
    assert_eq!(p.value("threads").unwrap().to_i32().unwrap(), 1);
    assert!(
        p.value("id_db_dir")
            .unwrap()
            .as_string_list()
            .unwrap()
            .is_empty()
    );
    assert!(!c.config_directory.exists());
    assert_eq!(c.get_openms_home_path(), c.home_directory);
    assert_eq!(c.get_openms_config_dir(), c.config_directory);
    assert_eq!(c.get_temp_directory().unwrap(), c.temporary_directory);
    assert_eq!(c.get_user_directory().unwrap(), c.home_directory);
    c.temporary_override = Some(d.path().join("override-temp"));
    c.user_override = Some(d.path().join("override-home"));
    assert_eq!(
        c.get_temp_directory().unwrap(),
        c.temporary_override.clone().unwrap()
    );
    assert_eq!(
        c.get_user_directory().unwrap(),
        c.user_override.clone().unwrap()
    );
    let alternative = d.path().join("alternative");
    let guard = c.get_temporary_file(Some(&alternative)).unwrap();
    assert_eq!(guard.path(), alternative);
    drop(guard);
    assert!(!alternative.exists());
    c.user_override = Some(PathBuf::new());
    assert_eq!(c.get_user_directory().unwrap(), Path::new("/"));
    c.temporary_override = None;
    let guard = c.get_temporary_file(None).unwrap();
    assert!(guard.path().exists());
}
#[test]
#[cfg(feature = "paramxml")]
fn source_stale_config_returns_original_tree_with_only_version_updated() {
    use openms::param::{Param, ParamValue};
    let d = temp();
    let mut c = context(d.path());
    fs::create_dir_all(&c.config_directory).unwrap();
    let config = c.config_directory.join("OpenMS.ini");
    let mut p = Param::new();
    p.set_value("version", ParamValue::from("old"), "", &[])
        .unwrap();
    p.set_value("custom", ParamValue::from(4711), "", &[])
        .unwrap();
    openms::format::paramxml::store(&config, &p).unwrap();
    let before = fs::read(&config).unwrap();
    let (actual, warnings) = c.get_system_parameters_with_warnings().unwrap();
    assert_eq!(actual.size(), 2);
    assert_eq!(actual.value("custom").unwrap().to_i32().unwrap(), 4711);
    assert_eq!(
        actual.value("version").unwrap().as_str().unwrap(),
        openms::CORE_SDK_VERSION
    );
    assert!(!actual.exists("threads").unwrap());
    assert!(warnings.iter().any(|s| s.contains("deprecated")));
    assert_eq!(fs::read(&config).unwrap(), before);
    let db = d.path().join("db");
    fs::create_dir(&db).unwrap();
    fs::write(db.join("protein.fasta"), b">p\nAA\n").unwrap();
    let data_path = d.path().join("data");
    data(&data_path);
    c.data_override = Some(data_path);
    p.set_value(
        "id_db_dir",
        ParamValue::StringList(vec![db.to_str().unwrap().into()]),
        "",
        &[],
    )
    .unwrap();
    p.set_value("temp_dir", ParamValue::from(" custom temp "), "", &[])
        .unwrap();
    p.set_value("home_dir", ParamValue::from(" custom home "), "", &[])
        .unwrap();
    openms::format::paramxml::store(&config, &p).unwrap();
    assert_eq!(
        c.find_database("protein.fasta").unwrap(),
        db.join("protein.fasta")
    );
    assert_eq!(c.get_temp_directory().unwrap(), Path::new(" custom temp "));
    assert_eq!(c.get_user_directory().unwrap(), Path::new(" custom home "));
    fs::write(&config, "broken XML").unwrap();
    c.temporary_override = Some(d.path().join("override"));
    assert!(c.get_temp_directory().is_err());
    assert!(c.get_user_directory().is_err());
}
#[test]
#[cfg(not(feature = "paramxml"))]
fn existing_config_without_xml_feature_is_an_explicit_error() {
    let d = temp();
    let c = context(d.path());
    fs::create_dir_all(&c.config_directory).unwrap();
    fs::write(c.config_directory.join("OpenMS.ini"), "<PARAMETERS/>").unwrap();
    assert!(matches!(
        c.get_system_parameters(),
        Err(openms::Error::Unsupported(_))
    ));
}
#[test]
fn native_resource_bounds_are_checked_before_directory_mutation() {
    let d = temp();
    let too_deep = d.path().join(
        std::iter::repeat_n("x", f::MAX_DEPTH + 1)
            .collect::<Vec<_>>()
            .join("/"),
    );
    assert!(f::make_dir(&too_deep).is_err());
    assert!(!d.path().join("x").exists());
    let long = "x".repeat(f::MAX_PATH_BYTES + 1);
    assert!(f::file_list(d.path(), &long, false).is_err());
    assert!(f::absolute_path(&long).is_err());
    let a = vec!["x".to_owned(); f::MAX_ENTRIES + 1];
    assert!(f::validate_matching_file_names(&a, &a, false, false).is_err());
    let c = context(d.path());
    assert!(c.find(&long, &[]).is_err());
}

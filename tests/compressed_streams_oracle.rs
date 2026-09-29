// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Differential test against the executed Release build.
//!
//! `tests/data/compressed_streams/f7b_oracle/records.tsv` is the output of
//! `../oracle/F7b-streams` (driver compiled against the Linux x86_64 Release
//! install `openms4-release-bc9cc12-c19e494-174b576`, run twice on ibminode06
//! with identical records; zlib 1.3.2, libbzip2 1.0.8, libzip 1.11.2). Each
//! stream case is replayed step by step and compared on the byte count,
//! `isOpen`, `streamEnd`, the bytes (hex, or length and FNV-1a 64) and the
//! exception class. The cases where the port deliberately differs are listed
//! in `DIVERGENT` and checked against the documented Rust behaviour instead.

#![cfg(feature = "file-compression")]

use openms::Error;
use openms::format::bzip2_ifstream::Bzip2Ifstream;
use openms::format::gzip_ifstream::GzipIfstream;
use openms::format::string_manager::StringManager;
use openms::format::zip_archive_file::ZipArchiveFile;
use openms::format::zip_ifstream::ZipIfstream;
use std::collections::BTreeMap;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/compressed_streams")
}

fn file(name: &str) -> PathBuf {
    let oracle = dir().join("f7b_oracle").join(name);
    if oracle.exists() {
        oracle
    } else {
        dir().join(name)
    }
}

fn records() -> BTreeMap<(String, String), String> {
    let text = std::fs::read_to_string(dir().join("f7b_oracle/records.tsv")).unwrap();
    let mut map = BTreeMap::new();
    for line in text.lines() {
        let mut fields = line.splitn(3, '\t');
        let case = fields.next().unwrap().to_owned();
        let step = fields.next().unwrap().to_owned();
        let result = fields.next().unwrap_or("").to_owned();
        assert!(
            map.insert((case, step), result).is_none(),
            "duplicate record"
        );
    }
    map
}

fn bytes(data: &[u8]) -> String {
    if data.len() <= 64 {
        return data.iter().map(|b| format!("{b:02x}")).collect();
    }
    let mut hash: u64 = 1_469_598_103_934_665_603;
    for &b in data {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(1_099_511_628_211);
    }
    format!("{}:{hash:016x}", data.len())
}

/// The C++ exception class a Rust error corresponds to.
fn class(error: &Error) -> &'static str {
    match error {
        Error::InvalidValue(_) => "IllegalArgument|InvalidValue",
        Error::Parse { .. } => "Parse Error",
        Error::Io(error) if error.kind() == ErrorKind::NotFound => "FileNotFound",
        Error::Io(error)
            if matches!(
                error.kind(),
                ErrorKind::InvalidData | ErrorKind::PermissionDenied
            ) =>
        {
            "FileNotReadable"
        }
        Error::Io(error) if error.kind() == ErrorKind::UnexpectedEof => "FileEmpty",
        other => panic!("unmapped error {other:?}"),
    }
}

fn same_class(oracle: &str, error: &Error) -> bool {
    let name = oracle
        .strip_prefix("ERR:")
        .and_then(|rest| rest.split(':').next())
        .unwrap_or("");
    class(error).split('|').any(|c| c == name)
}

/// The operations of one stream type, as the driver's `stream_case` uses them.
trait Stream: Sized {
    fn open_path(path: &Path) -> Result<Self, Error>;
    fn read_into(&mut self, buf: &mut [u8]) -> Result<usize, Error>;
    fn reopen(&mut self, path: &Path) -> Result<(), Error>;
    fn close_it(&mut self);
    fn is_open_(&self) -> bool;
    fn stream_end_(&self) -> bool;
}

macro_rules! stream {
    ($t:ty) => {
        impl Stream for $t {
            fn open_path(path: &Path) -> Result<Self, Error> {
                <$t>::from_path(path)
            }
            fn read_into(&mut self, buf: &mut [u8]) -> Result<usize, Error> {
                <$t>::read(self, buf)
            }
            fn reopen(&mut self, path: &Path) -> Result<(), Error> {
                self.open(path)
            }
            fn close_it(&mut self) {
                self.close()
            }
            fn is_open_(&self) -> bool {
                self.is_open()
            }
            fn stream_end_(&self) -> bool {
                self.stream_end()
            }
        }
    };
}
stream!(GzipIfstream);
stream!(Bzip2Ifstream);
stream!(ZipIfstream);

fn flags<S: Stream>(s: &S) -> String {
    format!(
        "open={} end={}",
        u8::from(s.is_open_()),
        u8::from(s.stream_end_())
    )
}

/// Replay a case; returns the Rust outcome of every step in the driver's
/// format, with errors as `ERR:<class>`.
fn replay<S: Stream>(
    name: &str,
    reads: &[i64],
    prepare: impl Fn(&mut S),
) -> Vec<(String, Result<String, Error>)> {
    let path = file(name);
    let mut out = Vec::new();
    let mut stream = match S::open_path(&path) {
        Ok(mut s) => {
            prepare(&mut s);
            out.push(("open".to_owned(), Ok(flags(&s))));
            s
        }
        Err(error) => {
            out.push(("open".to_owned(), Err(error)));
            return out;
        }
    };
    for (i, &n) in reads.iter().enumerate() {
        let step = i + 1;
        match n {
            -1 => {
                stream.close_it();
                out.push((format!("{step}:close"), Ok(flags(&stream))));
            }
            -2 => {
                let result = stream.reopen(&path).map(|()| flags(&stream));
                out.push((format!("{step}:reopen"), result));
            }
            n => {
                let mut buf = vec![0u8; n as usize];
                let result = stream
                    .read_into(&mut buf)
                    .map(|k| format!("n={k} {} data={}", flags(&stream), bytes(&buf[..k])));
                out.push((format!("{step}:read{n}"), result));
            }
        }
    }
    out
}

const BIG: i64 = 200_000;

/// (case, file, reads) exactly as `driver.cpp` runs them.
const GZIP: &[(&str, &str, &[i64])] = &[
    ("g_fixture", "GzipIfStream_1.gz", &[29, 10, 10]),
    ("g_exact", "GzipIfStream_1.gz", &[30, 10]),
    ("g_zero", "GzipIfStream_1.gz", &[0, 29, 0, 2]),
    ("g_corrupt", "GzipIfStream_1_corrupt.gz", &[10, 10, 10, 10]),
    ("g_corrupt30", "GzipIfStream_1_corrupt.gz", &[30, 30]),
    ("g_plain", "plain.txt", &[100]),
    ("g_one", "one.txt", &[10]),
    ("g_empty", "empty.bin", &[10]),
    ("g_members", "members.gz", &[1000]),
    ("g_payload", "payload.gz", &[BIG]),
    ("g_payload_exact", "payload.gz", &[80000, 1]),
    ("g_trunc_trailer", "trunc_trailer.gz", &[BIG, 10]),
    ("g_trunc_body", "trunc_body.gz", &[BIG, 10]),
    ("g_bad_crc", "bad_crc.gz", &[BIG, 10]),
    ("g_bad_crc_chunks", "bad_crc.gz", &[40000, 40000, 10]),
    ("g_bad_len", "bad_len.gz", &[BIG, 10]),
    ("g_hcrc_good", "hcrc_good.gz", &[100]),
    ("g_hcrc_bad", "hcrc_bad.gz", &[100, 10]),
    ("g_extra", "extra.gz", &[100]),
    ("g_bad_second", "bad_second.gz", &[100, 10]),
    ("g_reserved_flag", "reserved_flag.gz", &[100, 10]),
    ("g_missing", "missing.gz", &[10]),
    ("g_close", "GzipIfStream_1.gz", &[-1, 10, -2, 29]),
];

const BZIP2: &[(&str, &str, &[i64])] = &[
    ("b_fixture", "Bzip2IfStream_1.bz2", &[29, 10, 10]),
    ("b_exact", "Bzip2IfStream_1.bz2", &[30, 10]),
    ("b_zero", "Bzip2IfStream_1.bz2", &[0, 29]),
    ("b_corrupt", "Bzip2IfStream_1_corrupt.bz2", &[10, 10]),
    ("b_payload", "payload.bz2", &[BIG]),
    ("b_payload_exact", "payload.bz2", &[120000, 1]),
    ("b_payload_chunks", "payload.bz2", &[60000, 60000, 1]),
    ("b_two", "two.bz2", &[BIG]),
    ("b_trunc", "trunc.bz2", &[BIG, 10]),
    ("b_trunc_chunks", "trunc.bz2", &[100000, 100000]),
    ("b_plain", "plain.bz2", &[10]),
    ("b_empty", "empty.bin", &[10]),
    ("b_missing", "missing.bz2", &[10]),
];

const ZIP: &[(&str, &str, &[i64])] = &[
    ("z_fixture", "Deisotoper_input1.mzML.zip", &[28008, 10, 10]),
    ("z_short", "Deisotoper_input1.mzML.zip", &[30000, 10]),
    ("z_zero", "Deisotoper_input1.mzML.zip", &[0, -1, 0, 10]),
    ("z_bad_crc", "featurexml_bad_crc.zip", &[BIG, 10]),
    (
        "z_bad_crc_chunks",
        "featurexml_bad_crc.zip",
        &[1000, BIG, 10],
    ),
    ("z_bzip2", "featurexml_bzip2.zip", &[BIG, 10]),
    ("z_stored", "featurexml_stored.zip", &[BIG, 10]),
    ("z_lzma", "featurexml_lzma.zip", &[10]),
];

/// Cases the port answers differently on purpose; each is checked by its own
/// test below.
const DIVERGENT: &[&str] = &[
    // zlib 1.3.2 decodes up to 1 MiB ahead and reports corrupt deflate data
    // on the call after the one that decoded it; flate2 is not asked to
    // decode ahead, so the byte count delivered before the error differs.
    "g_corrupt",
    "g_corrupt30",
    // Truncation is an error by default; the source behaviour is opt-in.
    "g_trunc_trailer",
    "g_trunc_body",
    // The Release build drops a CRC-32 or ISIZE mismatch when the member is
    // decoded in one fetch and the read runs to the end; the port reports it.
    "g_bad_crc",
    "g_bad_crc_chunks",
    "g_bad_len",
    // A bad header CRC reads as an empty stream there, reserved header flags
    // as a plain file; the port refuses both.
    "g_hcrc_bad",
    "g_reserved_flag",
];

fn compare<S: Stream>(
    cases: &[(&str, &str, &[i64])],
    oracle: &BTreeMap<(String, String), String>,
) -> usize {
    let mut compared = 0;
    for &(case, name, reads) in cases {
        if DIVERGENT.contains(&case) {
            continue;
        }
        for (step, result) in replay::<S>(name, reads, |_| {}) {
            let expected = oracle
                .get(&(case.to_owned(), step.clone()))
                .unwrap_or_else(|| panic!("no oracle record for {case} {step}"));
            match result {
                Ok(text) => assert_eq!(&text, expected, "{case} {step}"),
                Err(error) => assert!(
                    same_class(expected, &error),
                    "{case} {step}: oracle {expected}, port {error:?}"
                ),
            }
            compared += 1;
        }
    }
    compared
}

#[test]
fn gzip_streams_match_the_release_build() {
    let oracle = records();
    assert_eq!(compare::<GzipIfstream>(GZIP, &oracle), 38);
}

#[test]
fn bzip2_streams_match_the_release_build() {
    let oracle = records();
    assert_eq!(compare::<Bzip2Ifstream>(BZIP2, &oracle), 35);
}

#[test]
fn zip_streams_match_the_release_build() {
    let oracle = records();
    assert_eq!(compare::<ZipIfstream>(ZIP, &oracle), 26);
}

fn case_reads(case: &str) -> (&'static str, &'static [i64]) {
    GZIP.iter()
        .find(|(c, _, _)| *c == case)
        .map(|&(_, f, r)| (f, r))
        .unwrap()
}

#[test]
fn gzip_truncation_matches_the_release_build_when_selected() {
    let oracle = records();
    for case in ["g_trunc_trailer", "g_trunc_body"] {
        let (name, reads) = case_reads(case);
        for (step, result) in replay::<GzipIfstream>(name, reads, |s| s.set_accept_truncated(true))
        {
            let expected = &oracle[&(case.to_owned(), step.clone())];
            match result {
                Ok(text) => assert_eq!(&text, expected, "{case} {step}"),
                Err(error) => assert!(same_class(expected, &error), "{case} {step}"),
            }
        }
        // By default the same input fails, and nothing claims a clean end.
        let steps = replay::<GzipIfstream>(name, reads, |_| {});
        assert!(
            steps
                .iter()
                .any(|(_, r)| matches!(r, Err(Error::Parse { .. }))),
            "{case}"
        );
    }
}

#[test]
fn gzip_integrity_failures_are_reported_where_the_release_build_is_silent() {
    let oracle = records();
    for case in ["g_bad_crc", "g_bad_len"] {
        let (name, reads) = case_reads(case);
        let steps = replay::<GzipIfstream>(name, reads, |_| {});
        // The same bytes are delivered first ...
        let expected = &oracle[&(case.to_owned(), "1:read200000".to_owned())];
        let data = |s: &str| s.split(" data=").nth(1).map(str::to_owned);
        let first = steps[1].1.as_ref().unwrap();
        assert_eq!(data(first), data(expected), "{case}");
        // ... then the mismatch is an error rather than a clean end.
        assert!(matches!(steps[2].1, Err(Error::Parse { .. })), "{case}");
    }
    // Chunked, the Release build reports it too, one call earlier.
    let (name, reads) = case_reads("g_bad_crc_chunks");
    let steps = replay::<GzipIfstream>(name, reads, |_| {});
    assert_eq!(
        steps[1].1.as_ref().unwrap(),
        &oracle[&("g_bad_crc_chunks".to_owned(), "1:read40000".to_owned())]
    );
    assert!(
        steps
            .iter()
            .any(|(_, r)| matches!(r, Err(Error::Parse { .. })))
    );

    for case in ["g_hcrc_bad", "g_reserved_flag", "g_corrupt", "g_corrupt30"] {
        let (name, reads) = case_reads(case);
        let steps = replay::<GzipIfstream>(name, reads, |_| {});
        assert!(
            steps
                .iter()
                .any(|(_, r)| matches!(r, Err(Error::Parse { .. }))),
            "{case}"
        );
    }
}

fn expect(oracle: &BTreeMap<(String, String), String>, step: &str) -> String {
    oracle[&("za".to_owned(), step.to_owned())].clone()
}

fn joined(entries: Vec<String>) -> String {
    entries.iter().map(|e| format!("{e}|")).collect()
}

fn one_line(text: &str) -> String {
    text.replace(['\n', '\t', '\r'], " ")
}

#[test]
fn zip_archive_file_matches_the_release_build() {
    let oracle = records();
    let base = std::env::temp_dir().join(format!("openms-f7b-za-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    let za = base.join("za");
    std::fs::create_dir_all(za.join("ztree/b")).unwrap();
    std::fs::create_dir_all(za.join("ztree/e")).unwrap();
    std::fs::write(za.join("ztree/a.parquet"), [1u8; 10]).unwrap();
    std::fs::write(za.join("ztree/b/c.parquet"), [2u8; 3]).unwrap();
    std::fs::write(za.join("src.txt"), b"replacement").unwrap();
    std::fs::write(za.join("text.zip"), b"plain").unwrap();
    std::fs::write(za.join("empty.zip"), b"").unwrap();
    let out = za.join("out.zip");

    ZipArchiveFile::zip_directory(za.join("ztree"), &out).unwrap();
    assert_eq!(
        joined(ZipArchiveFile::list_entries(&out)),
        expect(&oracle, "zipDirectory")
    );
    ZipArchiveFile::write_sidecar_index(&out).unwrap();
    assert_eq!(
        joined(ZipArchiveFile::list_entries(&out)),
        expect(&oracle, "writeSidecarIndex1")
    );
    let read_entry = |entry: &str| {
        let mut temp = None;
        let path = ZipArchiveFile::extract_entry_to_temp_file(&out, entry, &mut temp).unwrap();
        String::from_utf8(std::fs::read(path).unwrap()).unwrap()
    };
    assert_eq!(
        one_line(&read_entry("out.zip.idx.json")),
        expect(&oracle, "sidecar1")
    );
    ZipArchiveFile::write_sidecar_index(&out).unwrap();
    assert_eq!(
        joined(ZipArchiveFile::list_entries(&out)),
        expect(&oracle, "writeSidecarIndex2")
    );
    assert_eq!(
        one_line(&read_entry("out.zip.idx.json")),
        expect(&oracle, "sidecar2")
    );
    ZipArchiveFile::write_sidecar_index(za.join("ztree")).unwrap();
    assert_eq!(
        one_line(&std::fs::read_to_string(za.join("ztree.idx.json")).unwrap()),
        expect(&oracle, "sidecarDir")
    );

    let mut temp = None;
    let unpacked = ZipArchiveFile::unzip_directory(&out, &mut temp).unwrap();
    let mut lines = Vec::new();
    let mut pending = vec![unpacked.clone()];
    while let Some(d) = pending.pop() {
        for entry in std::fs::read_dir(&d).unwrap() {
            let entry = entry.unwrap();
            let rel = entry
                .path()
                .strip_prefix(&unpacked)
                .unwrap()
                .to_string_lossy()
                .replace('\\', "/");
            if entry.file_type().unwrap().is_dir() {
                lines.push(format!("{rel}/"));
                pending.push(entry.path());
            } else {
                lines.push(format!("{rel}={}", entry.metadata().unwrap().len()));
            }
        }
    }
    lines.sort();
    assert_eq!(
        format!(
            "{} {}",
            unpacked.file_name().unwrap().to_string_lossy(),
            joined(lines)
        ),
        expect(&oracle, "unzipDirectory")
    );

    let class_of = |step: &str, result: Result<(), Error>| {
        let expected = expect(&oracle, step);
        let error = result.expect_err(step);
        assert!(
            same_class(&expected, &error),
            "{step}: {expected} vs {error:?}"
        );
    };
    let mut none = None;
    assert_eq!(
        ZipArchiveFile::unzip_directory(za.join("ztree"), &mut none).unwrap(),
        za.join("ztree")
    );
    assert!(none.is_none());
    assert!(expect(&oracle, "unzipDirectoryDir").ends_with(" temp=0"));
    class_of(
        "unzipMissing",
        ZipArchiveFile::unzip_directory(za.join("missing.zip"), &mut None).map(|_| ()),
    );
    class_of(
        "unzipText",
        ZipArchiveFile::unzip_directory(za.join("text.zip"), &mut None).map(|_| ()),
    );
    ZipArchiveFile::add_or_replace_from_file(&out, "a.parquet", za.join("src.txt")).unwrap();
    assert_eq!(
        format!(
            "{} {}",
            joined(ZipArchiveFile::list_entries(&out)),
            read_entry("a.parquet")
        ),
        expect(&oracle, "addReplace")
    );
    let new = za.join("new.zip");
    ZipArchiveFile::add_or_replace_from_file(&new, "x/y.bin", za.join("src.txt")).unwrap();
    assert_eq!(
        joined(ZipArchiveFile::list_entries(&new)),
        expect(&oracle, "addNewArchive")
    );
    // Divergent: the header documents FileNotFound for a missing source, the
    // implementation fails later with InvalidValue; the port follows the
    // header, and the archive is untouched either way.
    assert!(expect(&oracle, "addMissingSource").starts_with("ERR:InvalidValue:"));
    let error =
        ZipArchiveFile::add_or_replace_from_file(&new, "z", za.join("nope.txt")).unwrap_err();
    assert!(matches!(&error, Error::Io(e) if e.kind() == ErrorKind::NotFound));
    assert_eq!(joined(ZipArchiveFile::list_entries(&new)), "x/y.bin|");
    class_of(
        "addToText",
        ZipArchiveFile::add_or_replace_from_file(za.join("text.zip"), "a", za.join("src.txt")),
    );
    class_of(
        "addToEmpty",
        ZipArchiveFile::add_or_replace_from_file(za.join("empty.zip"), "a", za.join("src.txt")),
    );
    assert_eq!(
        format!(
            "[{}]",
            joined(ZipArchiveFile::list_entries(za.join("text.zip")))
        ),
        expect(&oracle, "listText")
    );
    assert_eq!(
        format!(
            "[{}]",
            joined(ZipArchiveFile::list_entries(za.join("missing.zip")))
        ),
        expect(&oracle, "listMissing")
    );
    class_of(
        "extractMissingEntry",
        ZipArchiveFile::extract_entry_to_temp_file(&out, "nope", &mut None).map(|_| ()),
    );
    class_of(
        "extractMissingArchive",
        ZipArchiveFile::extract_entry_to_temp_file(za.join("none.zip"), "a", &mut None).map(|_| ()),
    );
    let mut none = None;
    assert_eq!(
        ZipArchiveFile::extract_entry_to_temp_file(za.join("ztree"), "b/c.parquet", &mut none)
            .unwrap(),
        za.join("ztree/b/c.parquet")
    );
    assert!(none.is_none());
    assert!(expect(&oracle, "extractDirLayout").ends_with("za/ztree/b/c.parquet temp=0"));
    class_of(
        "extractDirMissing",
        ZipArchiveFile::extract_entry_to_temp_file(za.join("ztree"), "nope", &mut None).map(|_| ()),
    );
    // Divergent: the source returns the escaping path for an unpacked layout.
    assert_eq!(expect(&oracle, "extractDirEscape"), "za/src.txt");
    assert!(
        ZipArchiveFile::extract_entry_to_temp_file(za.join("ztree"), "../src.txt", &mut None)
            .is_err()
    );
    class_of(
        "unzipTraversal",
        ZipArchiveFile::unzip_directory(file("traversal.zip"), &mut None).map(|_| ()),
    );
    let _ = std::fs::remove_dir_all(&base);
}

#[test]
fn string_manager_matches_the_release_build() {
    let oracle = records();
    let cases: &[(&str, &str)] = &[
        ("42", "42"),
        ("ws-17", "  -17\t\n"),
        ("+0", "+0"),
        ("vt5", "\u{b}5"),
        ("5vt", "5\u{b}"),
        ("ff5", "\u{c}5"),
        ("2147483647", "2147483647"),
        ("2147483648", "2147483648"),
        ("-2147483649", "-2147483649"),
        ("4294967296", "4294967296"),
        ("99999999999999999999", "99999999999999999999"),
        ("empty", ""),
        ("spaces", "   "),
        ("4 2", "4 2"),
        ("0x10", "0x10"),
        ("arabic1", "\u{661}"),
        ("plus space", "+ 5"),
    ];
    for &(label, text) in cases {
        let units: Vec<u16> = text.encode_utf16().collect();
        let expected = &oracle[&("parseInt".to_owned(), label.to_owned())];
        match StringManager::parse_int_source(&units) {
            Ok(value) => assert_eq!(&value.to_string(), expected, "{label}"),
            Err(_) => assert!(expected.starts_with("ERR:"), "{label}: {expected}"),
        }
        // The native parse agrees except where the source's (int) cast wraps
        // a value outside i32, which it refuses.
        let wrapped = ["2147483648", "-2147483649", "4294967296"].contains(&label);
        if wrapped {
            assert!(matches!(
                StringManager::parse_int(&units),
                Err(Error::InvalidRange(_))
            ));
        } else if let Ok(value) = expected.parse::<i32>() {
            assert_eq!(StringManager::parse_int(&units).unwrap(), value, "{label}");
        } else {
            assert!(StringManager::parse_int(&units).is_err(), "{label}");
        }
    }
    let utf8 = StringManager::convert_utf16(
        &"M\u{fc}ller \u{3b1} \u{1f9ea}"
            .encode_utf16()
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(
        bytes(utf8.as_bytes()),
        oracle[&("convert".to_owned(), "utf8".to_owned())]
    );
    let narrow: String = StringManager::convert("M\u{fc}ller")
        .iter()
        .map(|u| format!("{u:04x} "))
        .collect();
    assert_eq!(narrow, oracle[&("convert".to_owned(), "narrow".to_owned())]);
    // Divergent: Xerces pairs the high surrogate with the following 'A' and
    // emits U+2441; the port refuses the unpaired surrogate.
    assert_eq!(
        oracle[&("convert".to_owned(), "unpaired".to_owned())],
        "e29181"
    );
    assert!(StringManager::convert_utf16(&[0xD800, 0x41, 0]).is_err());
}

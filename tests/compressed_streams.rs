// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! GzipIfstream, Bzip2Ifstream, ZipIfstream and ZipArchiveFile.
//!
//! The sections marked "class test" port `GzipIfstream_test.cpp`,
//! `Bzip2Ifstream_test.cpp` and `ZipIfstream_test.cpp` at core `bc9cc12`,
//! with their literals as expected values; their data files are the upstream
//! bytes (tests/data/compressed_streams/provenance.json). The rest checks
//! invariants of the formats: a round trip reproduces its input, and the
//! gzip and bzip2 framing rules of RFC 1952 and libbzip2.

#![cfg(feature = "file-compression")]

use openms::Error;
use openms::format::bzip2_ifstream::Bzip2Ifstream;
use openms::format::gzip_ifstream::GzipIfstream;
use openms::format::zip_archive_file::{UNPACK_DIRECTORY, ZipArchiveFile};
use openms::format::zip_ifstream::ZipIfstream;
use std::io::{ErrorKind, Read, Write};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/compressed_streams")
        .join(name)
}

const MISSING: &str = "ThisFileDoesNotExist";
const TEXT: &str = "Was decompression successful?";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Self {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "openms-compressed-streams-{tag}-{}-{id}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }

    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        std::fs::write(&path, bytes).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn io_kind(error: Error) -> ErrorKind {
    match error {
        Error::Io(error) => error.kind(),
        other => panic!("expected an I/O error, got {other:?}"),
    }
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

fn bzip2(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

// ---------------------------------------------------------------- GzipIfstream

#[test]
fn gzip_class_test_default_and_detailed_constructor() {
    let stream = GzipIfstream::new();
    assert!(!stream.is_open());
    assert!(stream.stream_end());

    assert_eq!(
        io_kind(GzipIfstream::from_path(data(MISSING)).unwrap_err()),
        ErrorKind::NotFound
    );
    let mut gzip = GzipIfstream::from_path(data("GzipIfStream_1.gz")).unwrap();
    assert!(!gzip.stream_end());
    assert!(gzip.is_open());
    let mut buffer = [0u8; 29];
    assert_eq!(gzip.read(&mut buffer).unwrap(), 29);
    assert_eq!(&buffer, TEXT.as_bytes());
}

#[test]
fn gzip_class_test_open() {
    let mut gzip = GzipIfstream::new();
    assert_eq!(
        io_kind(gzip.open(data(MISSING)).unwrap_err()),
        ErrorKind::NotFound
    );
    gzip.open(data("GzipIfStream_1.gz")).unwrap();
    assert!(!gzip.stream_end());
    assert!(gzip.is_open());
    let mut buffer = [0u8; 29];
    assert_eq!(gzip.read(&mut buffer).unwrap(), 29);
    assert_eq!(&buffer, TEXT.as_bytes());
}

#[test]
fn gzip_class_test_read() {
    let _corrupt = GzipIfstream::from_path(data("GzipIfStream_1_corrupt.gz")).unwrap();
    let mut buffer = [0u8; 30];
    let mut gzip2 = GzipIfstream::from_path(data("GzipIfStream_1.gz")).unwrap();
    assert!(gzip2.is_open());
    gzip2.read(&mut buffer[..29]).unwrap();
    assert_eq!(gzip2.read(&mut buffer[..10]).unwrap(), 1);
    assert!(!gzip2.is_open());
    assert!(gzip2.stream_end());

    gzip2.open(data("GzipIfStream_1_corrupt.gz")).unwrap();
    gzip2.close();
    assert!(!gzip2.is_open());
    assert!(gzip2.stream_end());
    assert!(matches!(
        gzip2.read(&mut buffer[..10]),
        Err(Error::InvalidValue(_))
    ));
    gzip2.close();
    assert!(!gzip2.is_open());
    assert!(gzip2.stream_end());
    assert!(matches!(
        gzip2.read(&mut buffer[..10]),
        Err(Error::InvalidValue(_))
    ));

    gzip2.open(data("GzipIfStream_1.gz")).unwrap();
    for start in [0, 5, 10, 15, 20] {
        assert_eq!(gzip2.read(&mut buffer[start..start + 5]).unwrap(), 5);
    }
    assert_eq!(gzip2.read(&mut buffer[25..29]).unwrap(), 4);
    let mut end_of_file = [0u8; 2];
    assert_eq!(gzip2.read(&mut end_of_file).unwrap(), 1);
    assert!(gzip2.stream_end());
    assert_eq!(&buffer[..29], TEXT.as_bytes());
}

#[test]
fn gzip_corrupt_member_fails_and_closes() {
    // The fixture's deflate data is damaged ("invalid distance too far back");
    // reading it to the end must fail rather than yield the text.
    let mut gzip = GzipIfstream::from_path(data("GzipIfStream_1_corrupt.gz")).unwrap();
    let mut out = Vec::new();
    let mut buffer = [0u8; 10];
    let error = loop {
        match gzip.read(&mut buffer) {
            Ok(n) => {
                out.extend_from_slice(&buffer[..n]);
                assert!(gzip.is_open() || n < buffer.len());
                if !gzip.is_open() {
                    panic!("corrupt gzip read to its end: {out:?}");
                }
            }
            Err(error) => break error,
        }
    };
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    assert!(!gzip.is_open());
    assert!(gzip.stream_end());
    assert_ne!(out, TEXT.as_bytes());
}

#[test]
fn gzip_exact_request_leaves_the_stream_open_until_the_next_read() {
    // gzeof becomes true only after a read tried to go past the end.
    let mut gzip = GzipIfstream::from_path(data("GzipIfStream_1.gz")).unwrap();
    let mut buffer = [0u8; 30];
    assert_eq!(gzip.read(&mut buffer).unwrap(), 30);
    assert_eq!(&buffer[..29], TEXT.as_bytes());
    assert_eq!(buffer[29], b'\n');
    assert!(gzip.is_open());
    assert!(!gzip.stream_end());
    assert_eq!(gzip.read(&mut buffer).unwrap(), 0);
    assert!(!gzip.is_open());
    assert!(gzip.stream_end());
    // An empty request is a no-op even on a closed stream's Read impl.
    assert_eq!(Read::read(&mut gzip, &mut buffer).unwrap(), 0);
}

#[test]
fn gzip_plain_input_members_and_trailing_garbage() {
    let scratch = Scratch::new("gzip");
    // Transparent: non-gzip bytes pass through, an empty file is empty.
    for plain in [&b""[..], b"x", b"plain text, not gzip\n"] {
        let path = scratch.write("plain", plain);
        let mut out = Vec::new();
        GzipIfstream::from_path(&path)
            .unwrap()
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, plain);
    }
    // Concatenated members decode in sequence; trailing garbage is ignored.
    let mut bytes = gzip(b"first member\n");
    bytes.extend(gzip(b"second member\n"));
    bytes.extend(b"\0\0\0 trailing garbage");
    let path = scratch.write("members.gz", &bytes);
    let mut out = Vec::new();
    GzipIfstream::from_path(&path)
        .unwrap()
        .read_to_end(&mut out)
        .unwrap();
    assert_eq!(out, b"first member\nsecond member\n");
}

#[test]
fn gzip_truncation_fails_unless_the_source_behaviour_is_selected() {
    let scratch = Scratch::new("gzip-truncated");
    let payload: Vec<u8> = (0..20_000u32).flat_map(|i| i.to_le_bytes()).collect();
    let full = gzip(&payload);
    let path = scratch.write("truncated.gz", &full[..full.len() - 5]);
    let mut strict = GzipIfstream::from_path(&path).unwrap();
    let mut out = Vec::new();
    assert!(strict.read_to_end(&mut out).is_err());

    let mut source = GzipIfstream::from_path(&path).unwrap();
    source.set_accept_truncated(true);
    assert!(source.accept_truncated());
    let mut out = Vec::new();
    source.read_to_end(&mut out).unwrap();
    assert_eq!(out, payload);
    assert!(!source.is_open());

    // A damaged CRC-32 is corruption, not truncation, either way.
    let mut damaged = full.clone();
    let n = damaged.len();
    damaged[n - 8] ^= 1;
    let path = scratch.write("crc.gz", &damaged);
    let mut stream = GzipIfstream::from_path(&path).unwrap();
    stream.set_accept_truncated(true);
    let mut out = Vec::new();
    assert!(stream.read_to_end(&mut out).is_err());
}

// --------------------------------------------------------------- Bzip2Ifstream

#[test]
fn bzip2_class_test_default_and_detailed_constructor() {
    let stream = Bzip2Ifstream::new();
    assert!(!stream.is_open());
    assert!(stream.stream_end());

    assert_eq!(
        io_kind(Bzip2Ifstream::from_path(data(MISSING)).unwrap_err()),
        ErrorKind::NotFound
    );
    let mut bzip = Bzip2Ifstream::from_path(data("Bzip2IfStream_1.bz2")).unwrap();
    assert!(!bzip.stream_end());
    assert!(bzip.is_open());
    let mut buffer = [0u8; 29];
    assert_eq!(bzip.read(&mut buffer).unwrap(), 29);
    assert_eq!(&buffer, TEXT.as_bytes());
}

#[test]
fn bzip2_class_test_open() {
    let mut bzip = Bzip2Ifstream::new();
    assert_eq!(
        io_kind(bzip.open(data(MISSING)).unwrap_err()),
        ErrorKind::NotFound
    );
    bzip.open(data("Bzip2IfStream_1.bz2")).unwrap();
    assert!(!bzip.stream_end());
    assert!(bzip.is_open());
    let mut buffer = [0u8; 29];
    assert_eq!(bzip.read(&mut buffer).unwrap(), 29);
    assert_eq!(&buffer, TEXT.as_bytes());
}

#[test]
fn bzip2_class_test_read() {
    let mut buffer = [0u8; 30];
    let mut bzip = Bzip2Ifstream::from_path(data("Bzip2IfStream_1_corrupt.bz2")).unwrap();
    assert!(matches!(
        bzip.read(&mut buffer[..10]),
        Err(Error::Parse { .. })
    ));

    let mut bzip2 = Bzip2Ifstream::from_path(data("Bzip2IfStream_1.bz2")).unwrap();
    bzip2.read(&mut buffer[..29]).unwrap();
    assert_eq!(bzip2.read(&mut buffer[..10]).unwrap(), 1);
    assert!(!bzip2.is_open());
    assert!(bzip2.stream_end());

    bzip2.open(data("Bzip2IfStream_1_corrupt.bz2")).unwrap();
    assert!(matches!(
        bzip2.read(&mut buffer[..10]),
        Err(Error::Parse { .. })
    ));
    bzip2.close();
    assert!(!bzip2.is_open());
    assert!(bzip2.stream_end());
    assert!(matches!(
        bzip2.read(&mut buffer[..10]),
        Err(Error::InvalidValue(_))
    ));
    bzip2.close();
    assert!(!bzip2.is_open());
    assert!(bzip2.stream_end());
    assert!(matches!(
        bzip2.read(&mut buffer[..10]),
        Err(Error::InvalidValue(_))
    ));
    bzip2.open(data("Bzip2IfStream_1.bz2")).unwrap();
    assert_eq!(bzip2.read(&mut buffer[..29]).unwrap(), 29);
    assert_eq!(&buffer[..29], TEXT.as_bytes());
}

#[test]
fn bzip2_end_of_stream_closes_even_an_exact_request() {
    // BZ2_bzRead returns BZ_STREAM_END as soon as the end-of-stream marker is
    // decoded, and the source closes on it.
    let mut bzip = Bzip2Ifstream::from_path(data("Bzip2IfStream_1.bz2")).unwrap();
    let mut buffer = [0u8; 30];
    assert_eq!(bzip.read(&mut buffer).unwrap(), 30);
    assert_eq!(&buffer[..29], TEXT.as_bytes());
    assert!(!bzip.is_open());
    assert!(bzip.stream_end());
    // An empty request on an open stream is a no-op.
    let mut bzip = Bzip2Ifstream::from_path(data("Bzip2IfStream_1.bz2")).unwrap();
    assert_eq!(bzip.read(&mut []).unwrap(), 0);
    assert!(bzip.is_open());
}

#[test]
fn bzip2_reads_only_the_first_stream_and_refuses_truncation() {
    let scratch = Scratch::new("bzip2");
    let payload: Vec<u8> = (0..30_000u32).flat_map(|i| (i * 7).to_le_bytes()).collect();
    let mut bytes = bzip2(&payload);
    let single = bytes.len();
    bytes.extend(bzip2(b"second stream"));
    let path = scratch.write("two.bz2", &bytes);
    let mut out = Vec::new();
    Bzip2Ifstream::from_path(&path)
        .unwrap()
        .read_to_end(&mut out)
        .unwrap();
    assert_eq!(out, payload);

    let path = scratch.write("cut.bz2", &bytes[..single - 3]);
    let mut out = Vec::new();
    assert!(
        Bzip2Ifstream::from_path(&path)
            .unwrap()
            .read_to_end(&mut out)
            .is_err()
    );
    let path = scratch.write("empty.bz2", b"");
    let mut stream = Bzip2Ifstream::from_path(&path).unwrap();
    assert!(matches!(
        stream.read(&mut [0u8; 4]),
        Err(Error::Parse { .. })
    ));
    assert!(!stream.is_open());
}

// ----------------------------------------------------------------- ZipIfstream

#[test]
fn zip_class_test_open_read_stream_end_is_open_close() {
    let zif = ZipIfstream::new();
    assert!(!zif.is_open());
    assert!(zif.stream_end());

    let mut zif = ZipIfstream::new();
    zif.open(data("Deisotoper_input1.mzML.zip")).unwrap();
    assert!(zif.is_open());
    assert!(!zif.stream_end());
    zif.close();

    let mut zif = ZipIfstream::from_path(data("Deisotoper_input1.mzML.zip")).unwrap();
    let mut buf = [0u8; 255];
    let n = zif.read(&mut buf).unwrap();
    assert!(n > 0);
    assert!(String::from_utf8_lossy(&buf[..n]).contains("<?xml"));

    let mut zif = ZipIfstream::from_path(data("Deisotoper_input1.mzML.zip")).unwrap();
    assert!(!zif.stream_end());
    let mut buf = [0u8; 4096];
    while !zif.stream_end() {
        zif.read(&mut buf).unwrap();
    }
    assert!(zif.stream_end());

    let mut zif = ZipIfstream::new();
    assert!(!zif.is_open());
    zif.open(data("Deisotoper_input1.mzML.zip")).unwrap();
    assert!(zif.is_open());
    zif.close();
    assert!(!zif.is_open());

    let mut zif = ZipIfstream::from_path(data("Deisotoper_input1.mzML.zip")).unwrap();
    assert!(zif.is_open());
    zif.close();
    assert!(!zif.is_open());
    assert!(zif.stream_end());
}

#[test]
fn zip_class_test_extra_sections() {
    // Directory entry plus one file entry.
    let mut zif = ZipIfstream::new();
    zif.open(data("DirEntry.mzML.zip")).unwrap();
    assert!(zif.is_open());
    let mut buf = [0u8; 255];
    assert!(zif.read(&mut buf).unwrap() > 0);
    zif.close();

    // Empty archive: FileEmpty.
    let mut zif = ZipIfstream::new();
    assert_eq!(
        io_kind(zif.open(data("Empty.zip")).unwrap_err()),
        ErrorKind::UnexpectedEof
    );
    // Two file entries: ParseError naming both.
    match zif.open(data("MultiEntry.mzML.zip")) {
        Err(Error::Parse { message, .. }) => {
            assert!(message.contains("Deisotoper_input1.mzML"), "{message}");
            assert!(message.contains("Deisotoper_input2.mzML"), "{message}");
        }
        other => panic!("expected a parse error, got {other:?}"),
    }
    // Missing file: FileNotFound.
    assert_eq!(
        io_kind(zif.open("/nonexistent/path/to/file.zip").unwrap_err()),
        ErrorKind::NotFound
    );
    // Garbage: FileNotReadable.
    let scratch = Scratch::new("zip-corrupt");
    let corrupt = scratch.write(
        "corrupt.zip",
        b"This is definitely not a valid ZIP archive -- just plain text bytes.",
    );
    assert_eq!(
        io_kind(zif.open(&corrupt).unwrap_err()),
        ErrorKind::InvalidData
    );
    assert!(!zif.is_open());
}

#[test]
fn zip_entry_matches_the_plain_file_and_closes_on_zero() {
    let plain = std::fs::read(data("Deisotoper_input1.mzML")).unwrap();
    let mut zif = ZipIfstream::from_path(data("Deisotoper_input1.mzML.zip")).unwrap();
    let mut out = Vec::new();
    let mut buf = vec![0u8; 10_000];
    loop {
        let n = zif.read(&mut buf).unwrap();
        out.extend_from_slice(&buf[..n]);
        if n == 0 {
            break;
        }
        // A short read leaves the stream open; only a read of 0 closes it.
        assert!(zif.is_open());
    }
    assert_eq!(out, plain);
    assert!(!zif.is_open());
    assert!(zif.stream_end());
    // n == 0 is a no-op even when closed; otherwise a closed read fails.
    assert_eq!(zif.read(&mut []).unwrap(), 0);
    assert!(matches!(zif.read(&mut buf), Err(Error::InvalidValue(_))));
}

#[test]
fn zip_entry_with_damaged_data_fails() {
    let scratch = Scratch::new("zip-damaged");
    let mut bytes = std::fs::read(data("Deisotoper_input1.mzML.zip")).unwrap();
    // The local header of the only entry is 30 bytes plus its 22-byte name;
    // flip a byte well inside the deflate data.
    bytes[30 + 22 + 4000] ^= 0x55;
    let path = scratch.write("damaged.zip", &bytes);
    let mut zif = ZipIfstream::from_path(&path).unwrap();
    let mut out = Vec::new();
    assert!(zif.read_to_end(&mut out).is_err());
}

#[cfg(feature = "mzml")]
#[test]
fn file_handler_loads_mzml_from_zip_like_the_plain_file() {
    // ZipIfstream_test "[EXTRA] End-to-end: load mzML from ZIP via FileHandler".
    use openms::format::file_handler::FileHandler;
    use openms::format::file_types::FileType;
    assert_eq!(
        openms::format::file_types::type_by_file_name("sample.mzML.zip"),
        FileType::MzMl
    );
    assert_eq!(
        openms::format::file_types::type_by_file_name("data.featureXML.zip"),
        FileType::FeatureXml
    );
    assert_eq!(
        openms::format::file_types::type_by_file_name("ids.idXML.zip"),
        FileType::IdXml
    );
    let zipped = FileHandler::load_experiment(data("Deisotoper_input1.mzML.zip"), &[]).unwrap();
    let plain = FileHandler::load_experiment(data("Deisotoper_input1.mzML"), &[]).unwrap();
    assert_eq!(zipped.spectra.len(), plain.spectra.len());
    assert!(!plain.spectra.is_empty());
    for (z, p) in zipped.spectra.iter().zip(&plain.spectra) {
        assert_eq!(z.peaks.len(), p.peaks.len());
        assert_eq!(z.ms_level, p.ms_level);
    }
}

// ------------------------------------------------------------- ZipArchiveFile

#[test]
fn zip_directory_round_trips_through_unzip_and_list() {
    let scratch = Scratch::new("archive");
    let source = scratch.0.join("layout.oswpq");
    scratch.write(
        "layout.oswpq/precursors_features.parquet",
        b"precursor bytes",
    );
    scratch.write("layout.oswpq/sub/transition_features.parquet", &[7u8; 5000]);
    std::fs::create_dir_all(source.join("empty")).unwrap();
    let archive = scratch.0.join("results.oswpq.zip");
    std::fs::write(&archive, b"stale").unwrap();
    ZipArchiveFile::zip_directory(&source, &archive).unwrap();

    assert_eq!(
        ZipArchiveFile::list_entries(&archive),
        vec![
            "precursors_features.parquet".to_owned(),
            "sub/transition_features.parquet".to_owned(),
        ]
    );
    let raw = std::fs::File::open(&archive).unwrap();
    let mut listing = zip::ZipArchive::new(raw).unwrap();
    for i in 0..listing.len() {
        assert_eq!(
            listing.by_index_raw(i).unwrap().compression(),
            zip::CompressionMethod::Stored
        );
    }

    let mut temp = None;
    let unpacked = ZipArchiveFile::unzip_directory(&archive, &mut temp).unwrap();
    assert!(temp.is_some());
    assert_eq!(unpacked.file_name().unwrap(), UNPACK_DIRECTORY);
    assert_eq!(
        std::fs::read(unpacked.join("precursors_features.parquet")).unwrap(),
        b"precursor bytes"
    );
    assert_eq!(
        std::fs::read(unpacked.join("sub/transition_features.parquet")).unwrap(),
        vec![7u8; 5000]
    );
    // A directory is returned unchanged and leaves temp_dir alone.
    let mut untouched = None;
    assert_eq!(
        ZipArchiveFile::unzip_directory(&source, &mut untouched).unwrap(),
        source
    );
    assert!(untouched.is_none());
    // Unreadable input is FileNotFound; a non-archive is InvalidValue.
    assert_eq!(
        io_kind(
            ZipArchiveFile::unzip_directory(scratch.0.join("missing.zip"), &mut None).unwrap_err()
        ),
        ErrorKind::NotFound
    );
    let text = scratch.write("text.zip", b"not an archive");
    assert!(matches!(
        ZipArchiveFile::unzip_directory(&text, &mut None),
        Err(Error::InvalidValue(_))
    ));
    assert!(ZipArchiveFile::list_entries(&text).is_empty());
    assert!(ZipArchiveFile::list_entries(scratch.0.join("missing.zip")).is_empty());
}

#[test]
fn add_or_replace_keeps_positions_and_other_entries() {
    let scratch = Scratch::new("replace");
    let archive = scratch.0.join("store.zip");
    let a = scratch.write("a.bin", b"alpha");
    let b = scratch.write("b.bin", b"bravo");
    let c = scratch.write("c.bin", b"charlie, replaced");
    ZipArchiveFile::add_or_replace_from_file(&archive, "a.parquet", &a).unwrap();
    ZipArchiveFile::add_or_replace_from_file(&archive, "dir/b.parquet", &b).unwrap();
    ZipArchiveFile::add_or_replace_from_file(&archive, "a.parquet", &c).unwrap();
    assert_eq!(
        ZipArchiveFile::list_entries(&archive),
        vec!["a.parquet".to_owned(), "dir/b.parquet".to_owned()]
    );
    let mut temp = None;
    let a_out =
        ZipArchiveFile::extract_entry_to_temp_file(&archive, "a.parquet", &mut temp).unwrap();
    assert_eq!(std::fs::read(&a_out).unwrap(), b"charlie, replaced");
    let b_out =
        ZipArchiveFile::extract_entry_to_temp_file(&archive, "dir/b.parquet", &mut temp).unwrap();
    assert_eq!(std::fs::read(&b_out).unwrap(), b"bravo");
    assert!(b_out.starts_with(temp.as_ref().unwrap().path()));

    // Missing source file, entry and archive.
    let before = std::fs::read(&archive).unwrap();
    assert_eq!(
        io_kind(
            ZipArchiveFile::add_or_replace_from_file(&archive, "x", scratch.0.join("nope"))
                .unwrap_err()
        ),
        ErrorKind::NotFound
    );
    assert_eq!(std::fs::read(&archive).unwrap(), before);
    assert_eq!(
        io_kind(
            ZipArchiveFile::extract_entry_to_temp_file(&archive, "nope", &mut temp).unwrap_err()
        ),
        ErrorKind::NotFound
    );
    assert_eq!(
        io_kind(
            ZipArchiveFile::extract_entry_to_temp_file(scratch.0.join("none.zip"), "a", &mut None)
                .unwrap_err()
        ),
        ErrorKind::NotFound
    );
    // An existing non-archive is refused and left alone.
    let text = scratch.write("text.zip", b"plain");
    assert!(matches!(
        ZipArchiveFile::add_or_replace_from_file(&text, "a", &a),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(std::fs::read(&text).unwrap(), b"plain");
    // An empty file is not an archive either (Release oracle za/addToEmpty).
    let empty = scratch.write("empty.zip", b"");
    assert!(matches!(
        ZipArchiveFile::add_or_replace_from_file(&empty, "a", &a),
        Err(Error::InvalidValue(_))
    ));
    assert_eq!(std::fs::read(&empty).unwrap(), b"");
}

#[test]
fn sidecar_index_in_archives_and_directories() {
    let scratch = Scratch::new("sidecar");
    let source = scratch.0.join("run.oswpq");
    scratch.write("run.oswpq/a.parquet", &[1u8; 10]);
    scratch.write("run.oswpq/b/c.parquet", &[2u8; 3]);
    ZipArchiveFile::write_sidecar_index(&source).unwrap();
    let mut external = source.as_os_str().to_owned();
    external.push(".idx.json");
    assert_eq!(
        std::fs::read_to_string(PathBuf::from(external)).unwrap(),
        "{\n  \"a.parquet\": 10,\n  \"b/c.parquet\": 3\n}\n"
    );

    let archive = scratch.0.join("run.zip");
    std::fs::remove_file(scratch.0.join("run.oswpq.idx.json")).unwrap();
    ZipArchiveFile::zip_directory(&source, &archive).unwrap();
    ZipArchiveFile::write_sidecar_index(&archive).unwrap();
    let mut temp = None;
    let index = ZipArchiveFile::extract_entry_to_temp_file(&archive, "run.zip.idx.json", &mut temp)
        .unwrap();
    let first = "{\n  \"a.parquet\": 10,\n  \"b/c.parquet\": 3\n}\n";
    assert_eq!(std::fs::read_to_string(&index).unwrap(), first);
    // A second index lists the first one and replaces it in place.
    ZipArchiveFile::write_sidecar_index(&archive).unwrap();
    assert_eq!(
        ZipArchiveFile::list_entries(&archive),
        vec![
            "a.parquet".to_owned(),
            "b/c.parquet".to_owned(),
            "run.zip.idx.json".to_owned(),
        ]
    );
    let index = ZipArchiveFile::extract_entry_to_temp_file(&archive, "run.zip.idx.json", &mut temp)
        .unwrap();
    assert_eq!(
        std::fs::read_to_string(&index).unwrap(),
        format!(
            "{{\n  \"a.parquet\": 10,\n  \"b/c.parquet\": 3,\n  \"run.zip.idx.json\": {}\n}}\n",
            first.len()
        )
    );
    // The unpacked-layout branch of extract_entry_to_temp_file.
    assert_eq!(
        ZipArchiveFile::extract_entry_to_temp_file(&source, "b/c.parquet", &mut None).unwrap(),
        source.join("b/c.parquet")
    );
    assert!(ZipArchiveFile::extract_entry_to_temp_file(&source, "../run.zip", &mut None).is_err());
}

/// Every loader that reads through `path_io::open` reads a single-file ZIP
/// archive exactly as it reads the plain file (SHARED-02). Invariant: the
/// archive holds the plain file's bytes, deflated.
fn zipped(scratch: &Scratch, plain: &Path) -> PathBuf {
    let name = plain.file_name().unwrap().to_string_lossy().into_owned();
    let path = scratch.0.join(format!("{name}.zip"));
    let mut writer = zip::ZipWriter::new(std::fs::File::create(&path).unwrap());
    writer
        .start_file(
            name,
            zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated),
        )
        .unwrap();
    writer.write_all(&std::fs::read(plain).unwrap()).unwrap();
    writer.finish().unwrap();
    path
}

fn repo(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// `Debug` of a loaded value with the file path it was loaded from masked.
fn masked(value: &impl std::fmt::Debug, path: &Path) -> String {
    format!("{value:?}").replace(&path.display().to_string(), "<path>")
}

#[cfg(feature = "consensusxml")]
#[test]
fn consensusxml_loads_from_zip() {
    let scratch = Scratch::new("zip-consensus");
    let plain = repo("tests/data/consensusxml/ConsensusXMLFile_1.consensusXML");
    let archive = zipped(&scratch, &plain);
    let a = openms::format::consensusxml::load(&plain).unwrap();
    let b = openms::format::consensusxml::load(&archive).unwrap();
    assert!(!a.is_empty());
    assert_eq!(masked(&a, &plain), masked(&b, &archive));
}

#[cfg(feature = "idxml")]
#[test]
fn idxml_loads_from_zip() {
    let scratch = Scratch::new("zip-idxml");
    let plain = repo("tests/data/idxml_upstream_whole.idXML");
    let archive = zipped(&scratch, &plain);
    let a = openms::format::idxml::load(&plain).unwrap();
    assert!(!a.peptide_identifications.is_empty());
    assert_eq!(a, openms::format::idxml::load(&archive).unwrap());
}

#[cfg(feature = "paramxml")]
#[test]
fn paramxml_loads_from_zip() {
    let scratch = Scratch::new("zip-paramxml");
    let plain = repo("tests/data/baseline_filter_edges_benchmark.ini");
    let archive = zipped(&scratch, &plain);
    let a = openms::format::paramxml::load(&plain).unwrap();
    assert_eq!(a, openms::format::paramxml::load(&archive).unwrap());
}

#[cfg(feature = "cv-mapping")]
#[test]
fn cv_mapping_loads_from_zip() {
    let scratch = Scratch::new("zip-cv-mapping");
    let plain = repo("tests/data/cv_mapping/ms-mapping.xml");
    let archive = zipped(&scratch, &plain);
    let reader = openms::format::cv_mapping::CVMappingFile {
        options: Default::default(),
    };
    let a = reader.load(&plain).unwrap();
    let b = reader.load(&archive).unwrap();
    assert_eq!(masked(&a, &plain), masked(&b, &archive));
}

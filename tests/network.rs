// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/Network.h`: the single section of `Network_test.cpp`, plus the naming
//! rules of the source's file-static `saveFileName_`, which that section only
//! touches through one file name.
//!
//! The upstream section downloads `file://<test data>/Network_test_fixture.txt`
//! into the temporary directory and asserts the file appeared. `ureq` speaks
//! HTTP and HTTPS only, so the `file://` scheme has no counterpart; the same
//! fixture's bytes are served through a recorded transport instead, and the
//! assertions are strengthened from "the file exists" to "the file exists at the
//! derived name and holds exactly the fixture's bytes". The fixture itself is
//! the upstream one, byte for byte, and is hashed in
//! `tests/data/network_provenance.json`.
#![cfg(feature = "network")]

use openms::system::file::TempDir;
use openms::system::network::{
    DOWNLOAD_TIMEOUT_SECONDS, FALLBACK_BASENAME, MAX_URL_BYTES, download_file_with, save_file_name,
};
use openms::system::network_get_request::{
    HttpTransport, TransportError, TransportRequest, TransportResponse,
};
use std::path::Path;
use std::sync::Mutex;
use std::time::Duration;

/// The upstream `Network_test_fixture.txt`, embedded so the test needs no path.
const FIXTURE: &[u8] = include_bytes!("data/network_test_fixture.txt");

/// Serves one recorded body, and remembers the timeout it was asked for.
struct Serving {
    body: Vec<u8>,
    status: u16,
    seen_timeout: Mutex<Option<Option<Duration>>>,
}

impl Serving {
    fn new(status: u16, body: &[u8]) -> Self {
        Self {
            body: body.to_vec(),
            status,
            seen_timeout: Mutex::new(None),
        }
    }
}

impl HttpTransport for Serving {
    fn get(&self, request: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        if let Ok(mut slot) = self.seen_timeout.lock() {
            *slot = Some(request.timeout);
        }
        Ok(TransportResponse {
            status: self.status,
            headers: Vec::new(),
            body: self.body.clone(),
        })
    }
}

struct Failing(TransportError);

impl HttpTransport for Failing {
    fn get(&self, _: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        Err(self.0.clone())
    }
}

/// Class-test section
/// `static void downloadFile(const std::string& url, const std::string& download_folder)`.
///
/// The upstream section asserts `File::exists(folder + "/Network_test_fixture.txt")`
/// after downloading that fixture. Reproduced here, with the same derived name,
/// the same fixture bytes on disk, and the source's ten-minute timeout observed
/// at the transport.
#[test]
fn download_file_writes_the_fixture_under_its_url_basename() -> openms::Result<()> {
    let folder = TempDir::new(false)?;
    let transport = Serving::new(200, FIXTURE);
    let url = "https://openms.invalid/data/Network_test_fixture.txt";

    let written = download_file_with(&transport, url, folder.path())?;

    assert_eq!(written, folder.path().join("Network_test_fixture.txt"));
    assert!(openms::system::file::exists(&written));
    assert_eq!(std::fs::read(&written)?, FIXTURE);
    assert_eq!(
        transport.seen_timeout.lock().unwrap().unwrap(),
        Some(Duration::from_secs(
            u64::try_from(DOWNLOAD_TIMEOUT_SECONDS).unwrap()
        ))
    );
    Ok(())
}

/// Native: the never-overwrite rule of `saveFileName_`, exercised on disk.
///
/// The source appends `.0`, `.1`, `.2`, … until a name is free and never
/// replaces an existing file. Three downloads of the same URL therefore produce
/// three files.
#[test]
fn repeated_downloads_never_overwrite() -> openms::Result<()> {
    let folder = TempDir::new(false)?;
    let url = "https://openms.invalid/data/run.tsv";

    let first = download_file_with(&Serving::new(200, b"one"), url, folder.path())?;
    let second = download_file_with(&Serving::new(200, b"two"), url, folder.path())?;
    let third = download_file_with(&Serving::new(200, b"three"), url, folder.path())?;

    assert_eq!(first, folder.path().join("run.tsv"));
    assert_eq!(second, folder.path().join("run.tsv.0"));
    assert_eq!(third, folder.path().join("run.tsv.1"));
    assert_eq!(std::fs::read(&first)?, b"one");
    assert_eq!(std::fs::read(&second)?, b"two");
    assert_eq!(std::fs::read(&third)?, b"three");
    Ok(())
}

/// Native: `saveFileName_`'s derivation rules, without touching the network.
#[test]
fn the_destination_name_follows_the_source_rules() -> openms::Result<()> {
    let empty = Path::new("");
    assert_eq!(save_file_name("https://h/a/b/c.tsv", empty)?, "c.tsv");
    // The query is cut before the fragment, so a '#' inside a query is gone too.
    assert_eq!(save_file_name("https://h/c.tsv?v=2", empty)?, "c.tsv");
    assert_eq!(save_file_name("https://h/c.tsv#top", empty)?, "c.tsv");
    assert_eq!(save_file_name("https://h/c.tsv?a=#b", empty)?, "c.tsv");
    // An empty basename falls back to the source's literal.
    assert_eq!(save_file_name("https://h/", empty)?, FALLBACK_BASENAME);
    assert_eq!(save_file_name("https://h/?x=1", empty)?, FALLBACK_BASENAME);
    // std::filesystem::path("https://host").filename() is "host".
    assert_eq!(save_file_name("https://host", empty)?, "host");
    // Native guards the source does not have.
    assert!(save_file_name("https://h/a/..", empty).is_err());
    assert!(save_file_name("https://h/a/.", empty).is_err());
    assert!(save_file_name(&format!("https://h/{}", "x".repeat(MAX_URL_BYTES)), empty).is_err());
    Ok(())
}

/// Native: an existing file in the folder pushes the name to the first free suffix.
#[test]
fn an_existing_name_is_skipped_in_suffix_order() -> openms::Result<()> {
    let folder = TempDir::new(false)?;
    assert_eq!(save_file_name("https://h/x.txt", folder.path())?, "x.txt");
    std::fs::write(folder.path().join("x.txt"), b"a")?;
    assert_eq!(save_file_name("https://h/x.txt", folder.path())?, "x.txt.0");
    std::fs::write(folder.path().join("x.txt.0"), b"b")?;
    std::fs::write(folder.path().join("x.txt.1"), b"c")?;
    assert_eq!(save_file_name("https://h/x.txt", folder.path())?, "x.txt.2");
    Ok(())
}

/// Class-test-adjacent: the source's documented failure paths.
///
/// The header lists three of them; the download error is the one reachable
/// without breaking the filesystem, and it carries the source's message text.
/// An HTTP status of 400 or above is a download error too, because
/// `NetworkGetRequest` classifies it as one.
#[test]
fn a_failed_request_is_an_io_error_naming_the_url() -> openms::Result<()> {
    let folder = TempDir::new(false)?;
    let url = "https://openms.invalid/data/run.tsv";

    let error = download_file_with(&Failing(TransportError::HostNotFound), url, folder.path())
        .expect_err("a transport failure must not produce a file");
    let text = error.to_string();
    assert!(text.contains(url), "{text}");
    assert!(text.contains("failed!"), "{text}");

    let status_error = download_file_with(&Serving::new(404, b"gone"), url, folder.path())
        .expect_err("an HTTP error status must not produce a file");
    assert!(status_error.to_string().contains("HTTP error 404"));

    // Neither attempt left anything behind.
    assert!(!openms::system::file::exists(folder.path().join("run.tsv")));
    assert!(!openms::system::file::exists(
        folder.path().join("run.tsv.0")
    ));
    Ok(())
}

/// Native: an empty destination folder means the current directory.
///
/// The source maps `""` to `"./"`; this maps it to `"."`, which names the same
/// place. Asserted through the returned path rather than by writing into the
/// working directory, which a parallel test run shares.
#[test]
fn an_empty_folder_means_the_current_directory() -> openms::Result<()> {
    assert_eq!(
        save_file_name("https://h/unique-name.tsv", Path::new(""))?,
        "unique-name.tsv"
    );
    Ok(())
}

/// Native: a zero-length response still produces the file.
///
/// The source writes `data.data()` for `data.size()` bytes with no emptiness
/// test, so an empty body creates an empty file rather than nothing.
#[test]
fn an_empty_body_still_creates_the_file() -> openms::Result<()> {
    let folder = TempDir::new(false)?;
    let written = download_file_with(
        &Serving::new(200, b""),
        "https://h/empty.bin",
        folder.path(),
    )?;
    assert!(openms::system::file::exists(&written));
    assert_eq!(std::fs::read(&written)?.len(), 0);
    Ok(())
}

// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/UpdateCheck.h`: the single section of `UpdateCheck_test.cpp`, plus
//! the rate-limiting and version-comparison logic it cannot reach.
//!
//! The upstream section proves one thing: when the OpenMS config directory
//! cannot be created, `run` warns and returns *before* issuing any network
//! request, without throwing. It arranges that by pointing `XDG_CONFIG_HOME` at
//! a path whose parent is a regular file, so `create_directories` fails with
//! `ENOTDIR` for every user, root included. The same arrangement is reproduced
//! here without touching the process environment, because
//! `system::file::FileContext` carries the config directory as a value.
//!
//! Every other test drives a recorded transport. Nothing here opens a socket,
//! and in particular nothing queries the real OpenMS REST server — a test that
//! did would report this machine's usage statistics upstream.
#![cfg(feature = "network")]

use openms::system::file::{FileContext, TempDir};
use openms::system::network_get_request::{
    HttpTransport, RequestError, TransportError, TransportRequest, TransportResponse,
};
use openms::system::update_check::{
    MIN_QUERY_INTERVAL, QUERY_TIMEOUT_SECONDS, SemanticVersion, UPDATE_URL_PREFIX,
    UpdateCheckOutcome, architecture_tag, is_due, platform_tag, query_url, run,
    should_report_update, tool_version_string, version_file_path,
};
use std::path::Path;
use std::sync::Mutex;
use std::time::{Duration, SystemTime};

/// Answers with one recorded body and remembers what it was asked.
struct Recorded {
    body: Vec<u8>,
    status: u16,
    seen: Mutex<Vec<(String, Option<Duration>)>>,
}

impl Recorded {
    fn new(body: &[u8]) -> Self {
        Self {
            body: body.to_vec(),
            status: 200,
            seen: Mutex::new(Vec::new()),
        }
    }

    fn calls(&self) -> Vec<(String, Option<Duration>)> {
        self.seen.lock().unwrap().clone()
    }
}

impl HttpTransport for Recorded {
    fn get(&self, request: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        self.seen
            .lock()
            .unwrap()
            .push((request.url.to_owned(), request.timeout));
        Ok(TransportResponse {
            status: self.status,
            headers: Vec::new(),
            body: self.body.clone(),
        })
    }
}

/// Fails every request, and counts how often it was asked.
#[derive(Default)]
struct Failing(Mutex<usize>);

impl HttpTransport for Failing {
    fn get(&self, _: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        *self.0.lock().unwrap() += 1;
        Err(TransportError::Timeout)
    }
}

/// Refuses to be called at all; any request is a test failure.
struct Forbidden;

impl HttpTransport for Forbidden {
    fn get(&self, _: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        panic!("update check issued a request it should not have");
    }
}

fn context_in(root: &Path, config: &Path) -> openms::Result<FileContext> {
    let mut context = FileContext::new(root, root, root)?;
    context.config_directory = config.to_path_buf();
    Ok(context)
}

/// Class-test section
/// `(static void run(const std::string& tool_name, const std::string& version, int debug_level))`.
///
/// An uncreatable config directory must warn and return before any request. The
/// upstream test asserts two things — nothing thrown, and a warning containing
/// `"Could not create config directory"` — and both are asserted here, together
/// with the stronger property the upstream test can only imply: the transport
/// was never called.
#[test]
fn an_uncreatable_config_directory_warns_and_makes_no_request() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    // A regular file standing where a directory would have to be: creating
    // "<blocker>/OpenMS" then fails with ENOTDIR for every user, root included.
    let blocker = root.path().join("not-a-directory");
    std::fs::write(&blocker, b"not a directory")?;
    let context = context_in(root.path(), &blocker.join("OpenMS"))?;

    let report = run(
        &context,
        &Forbidden,
        &SemanticVersion::parse("3.0.0").unwrap(),
        "UpdateCheck_test_tool",
        "1.0.0",
        0,
    )?;

    assert_eq!(
        report.outcome,
        UpdateCheckOutcome::ConfigDirectoryUnavailable
    );
    assert_eq!(report.warnings.len(), 1);
    assert!(
        report.warnings[0].contains("Could not create config directory"),
        "{}",
        report.warnings[0]
    );
    assert!(report.notices.is_empty());
    Ok(())
}

/// Native: the first run creates the stamp file, queries, and decides.
///
/// The source's own flow, end to end: no `.ver` file, so one is created and the
/// run is treated as first; the query goes to the REST endpoint for the
/// platform-tagged identifier with a five-second timeout; a newer server version
/// produces the announcement.
#[test]
fn a_first_run_stamps_queries_and_reports_an_update() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    let context = context_in(root.path(), &root.path().join("config"))?;
    let transport = Recorded::new(b"3.2.0");

    let report = run(
        &context,
        &transport,
        &SemanticVersion::parse("3.1.0").unwrap(),
        "FeatureFinderCentroided",
        "2.0.0",
        1,
    )?;

    assert_eq!(report.outcome, UpdateCheckOutcome::UpdateAvailable);
    assert!(report.warnings.is_empty());

    let stamp = version_file_path(&context, "FeatureFinderCentroided")?;
    assert!(openms::system::file::exists(&stamp));

    let calls = transport.calls();
    assert_eq!(calls.len(), 1);
    let expected = query_url(&tool_version_string(
        "FeatureFinderCentroided",
        "2.0.0",
        platform_tag(),
        architecture_tag(),
    ));
    assert_eq!(calls[0].0, expected);
    assert!(expected.starts_with(UPDATE_URL_PREFIX));
    assert_eq!(
        calls[0].1,
        Some(Duration::from_secs(
            u64::try_from(QUERY_TIMEOUT_SECONDS).unwrap()
        ))
    );

    // The announcement is not gated on debug_level, and names the *local*
    // version - the source prints its `version` argument, not the server's.
    assert!(
        report
            .notices
            .iter()
            .any(|line| line
                == "Version 2.0.0 of FeatureFinderCentroided is available at www.OpenMS.de"),
        "{:?}",
        report.notices
    );
    // The three privacy lines appear because debug_level > 0.
    assert!(
        report
            .notices
            .iter()
            .any(|line| line.contains("The OpenMS team is collecting usage statistics"))
    );
    assert!(
        report
            .notices
            .iter()
            .any(|line| line.contains("OPENMS_DISABLE_UPDATE_CHECK"))
    );
    Ok(())
}

/// Native: a second run inside the 24-hour window makes no request.
///
/// The stamp file's modification time was set to now by the first run, so the
/// second is not due. `Forbidden` turns any request into a panic.
#[test]
fn a_second_run_inside_the_window_is_not_due() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    let context = context_in(root.path(), &root.path().join("config"))?;
    let first = run(
        &context,
        &Recorded::new(b"3.0.0"),
        &SemanticVersion::parse("3.0.0").unwrap(),
        "TOPPTool",
        "3.0.0",
        0,
    )?;
    assert_eq!(first.outcome, UpdateCheckOutcome::UpToDate);

    let second = run(
        &context,
        &Forbidden,
        &SemanticVersion::parse("3.0.0").unwrap(),
        "TOPPTool",
        "3.0.0",
        0,
    )?;
    assert_eq!(second.outcome, UpdateCheckOutcome::NotDue);
    Ok(())
}

/// Native: a failed query still consumes the day's budget.
///
/// The source bumps the stamp's modification time *before* it issues the
/// request, so a server that is down does not cause a retry on every invocation.
/// That ordering is behaviour a caller can observe, so it is pinned.
#[test]
fn a_failed_query_still_consumes_the_window() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    let context = context_in(root.path(), &root.path().join("config"))?;
    let running = SemanticVersion::parse("3.0.0").unwrap();
    let transport = Failing::default();

    let first = run(&context, &transport, &running, "TOPPTool", "3.0.0", 1)?;
    assert_eq!(
        first.outcome,
        UpdateCheckOutcome::QueryFailed(RequestError::Transport(TransportError::Timeout))
    );
    assert!(
        first
            .notices
            .iter()
            .any(|line| line.contains("Connecting to REST server failed"))
    );

    let second = run(&context, &transport, &running, "TOPPTool", "3.0.0", 1)?;
    assert_eq!(second.outcome, UpdateCheckOutcome::NotDue);
    assert_eq!(*transport.0.lock().unwrap(), 1);
    Ok(())
}

/// Native: a response that is not a version is ignored, not acted on.
///
/// `VersionDetails::create` returns `EMPTY` for anything it cannot parse, and
/// the source skips the comparison for `EMPTY`. A literal `0.0.0` takes the same
/// branch, because the source cannot tell it from a parse failure.
#[test]
fn an_unusable_server_answer_reports_nothing() -> openms::Result<()> {
    for body in [
        &b"<html>404</html>"[..],
        &b""[..],
        &b"0.0.0"[..],
        &b"not.a.version"[..],
        &[0xff, 0xfe][..],
    ] {
        let root = TempDir::new(false)?;
        let context = context_in(root.path(), &root.path().join("config"))?;
        let report = run(
            &context,
            &Recorded::new(body),
            &SemanticVersion::default(),
            "TOPPTool",
            "1.0.0",
            0,
        )?;
        assert_eq!(
            report.outcome,
            UpdateCheckOutcome::ServerVersionUnusable,
            "body {body:?}"
        );
        assert!(report.notices.is_empty());
    }
    Ok(())
}

/// Native: no notice is emitted at `debug_level == 0` except the announcement.
#[test]
fn debug_level_zero_suppresses_everything_but_the_announcement() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    let context = context_in(root.path(), &root.path().join("config"))?;
    let report = run(
        &context,
        &Recorded::new(b"9.9.9"),
        &SemanticVersion::parse("1.0.0").unwrap(),
        "TOPPTool",
        "1.0.0",
        0,
    )?;
    assert_eq!(report.outcome, UpdateCheckOutcome::UpdateAvailable);
    assert_eq!(
        report.notices,
        vec!["Version 1.0.0 of TOPPTool is available at www.OpenMS.de".to_owned()]
    );
    Ok(())
}

/// Native: a tool name that is not a plain file-name component is refused.
///
/// The source concatenates the name straight into a path, so `../../x` writes
/// outside the config directory.
#[test]
fn a_traversing_tool_name_is_refused() -> openms::Result<()> {
    let root = TempDir::new(false)?;
    let context = context_in(root.path(), &root.path().join("config"))?;
    for name in ["", "../escape", "a/b", ".hidden"] {
        assert!(
            run(
                &context,
                &Forbidden,
                &SemanticVersion::default(),
                name,
                "1.0.0",
                0
            )
            .is_err(),
            "tool name {name:?} must be refused"
        );
    }
    Ok(())
}

/// `VersionInfo::VersionDetails::create`, rule by rule.
///
/// Transcribed from the parse rules the header documents and the `.cpp`
/// implements; the boundary cases are derived from `StringUtils::toInt32`'s
/// whole-string requirement rather than taken from any literal.
#[test]
fn version_parsing_follows_version_details_create() {
    assert!(SemanticVersion::parse("3").is_none());
    assert!(SemanticVersion::parse("").is_none());
    assert!(SemanticVersion::parse("x.y").is_none());

    let two = SemanticVersion::parse("2.7").unwrap();
    assert_eq!((two.major, two.minor, two.patch), (2, 7, 0));
    assert!(two.pre_release.is_empty());

    let three = SemanticVersion::parse("1.2.3").unwrap();
    assert_eq!((three.major, three.minor, three.patch), (1, 2, 3));

    let pre = SemanticVersion::parse("1.2.3-alpha.1").unwrap();
    assert_eq!(pre.pre_release, "alpha.1");

    // The dash is only looked for after the second dot.
    assert!(SemanticVersion::parse("1.2-alpha").is_none());
    // toInt32 allows surrounding whitespace and a leading '+'.
    assert_eq!(SemanticVersion::parse(" 1 . 2 . 3 \n").unwrap().patch, 3);
    assert_eq!(SemanticVersion::parse("+1.+2.+3").unwrap().major, 1);
    // ...and rejects trailing rubbish.
    assert!(SemanticVersion::parse("1.2.3x").is_none());
    // Int32 range.
    assert!(SemanticVersion::parse("2147483648.0.0").is_none());
    assert_eq!(
        SemanticVersion::parse("2147483647.0.0").unwrap().major,
        2_147_483_647
    );
}

/// `VersionInfo::VersionDetails::operator<`, including its inconsistency.
///
/// Derived from the relation in `VersionInfo.cpp`: lexicographic on the triple,
/// and a pre-release sorts below the same triple without one. Two different
/// pre-releases on one triple are neither less nor greater, yet compare unequal,
/// which is why `PartialOrd` is deliberately absent.
#[test]
fn version_ordering_is_the_source_relation() {
    let v = |s: &str| SemanticVersion::parse(s).unwrap();
    assert!(v("1.0.0").is_less_than(&v("2.0.0")));
    assert!(v("1.0.0").is_less_than(&v("1.1.0")));
    assert!(v("1.0.0").is_less_than(&v("1.0.1")));
    assert!(!v("2.0.0").is_less_than(&v("1.9.9")));
    assert!(v("1.0.0-rc1").is_less_than(&v("1.0.0")));
    assert!(!v("1.0.0").is_less_than(&v("1.0.0-rc1")));
    assert!(!v("1.0.0-a").is_less_than(&v("1.0.0-b")));
    assert!(!v("1.0.0-b").is_less_than(&v("1.0.0-a")));
    assert_ne!(v("1.0.0-a"), v("1.0.0-b"));
    assert!(v("1.0.0-a").is_greater_than(&v("1.0.0-b")));
    assert!(v("2.0.0").is_greater_than(&v("1.0.0")));
    assert!(!v("1.0.0").is_greater_than(&v("1.0.0")));
    assert!(SemanticVersion::default().is_empty());
    assert!(!v("0.0.1").is_empty());
}

/// The decision function, independent of any transport or filesystem.
#[test]
fn the_update_decision_matches_the_source_branches() {
    let running = SemanticVersion::parse("3.0.0").unwrap();
    assert!(should_report_update(&running, "3.0.1"));
    assert!(should_report_update(&running, "3.1.0"));
    assert!(should_report_update(&running, "4.0.0"));
    assert!(should_report_update(&running, "3.1.0\n"));
    assert!(!should_report_update(&running, "3.0.0"));
    assert!(!should_report_update(&running, "2.9.9"));
    assert!(!should_report_update(&running, "0.0.0"));
    assert!(!should_report_update(&running, "garbage"));
    // A pre-release of the next patch is still newer than 3.0.0.
    assert!(should_report_update(&running, "3.0.1-rc1"));
    // A pre-release of the running version is not.
    assert!(!should_report_update(&running, "3.0.0-rc1"));
}

/// The rate-limit predicate, at its exact boundary.
#[test]
fn due_matches_the_source_comparison() {
    let now = SystemTime::UNIX_EPOCH + Duration::from_secs(30 * 24 * 60 * 60);
    assert!(is_due(true, Some(now), now));
    assert!(!is_due(false, Some(now), now));
    assert!(!is_due(false, Some(now - MIN_QUERY_INTERVAL), now));
    assert!(is_due(
        false,
        Some(now - MIN_QUERY_INTERVAL - Duration::from_secs(1)),
        now
    ));
    assert!(is_due(false, None, now));
    assert!(!is_due(false, Some(now + MIN_QUERY_INTERVAL), now));
}

/// The identifier and URL the server is asked about.
///
/// The source's own comment gives the example transcribed here.
#[test]
fn the_query_identifier_is_the_documented_shape() {
    assert_eq!(
        tool_version_string("FeatureFinderCentroided", "2.0.0", "Win", "64"),
        "OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0"
    );
    assert_eq!(
        query_url("OpenMS_Default_Win_64_X_1.0"),
        "http://openms-update.cs.uni-tuebingen.de/check/OpenMS_Default_Win_64_X_1.0"
    );
    assert!(["Win", "Mac", "Linux", "Unix", "unknown"].contains(&platform_tag()));
    assert!(["32", "64"].contains(&architecture_tag()));
}

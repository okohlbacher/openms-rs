// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/NetworkGetRequest.h` and `SYSTEM/CurlInit.h`: the nine sections of
//! `NetworkGetRequest_test.cpp`, plus the transport classification the source
//! delegates to libcurl.
//!
//! No test here reaches the network, and the property is checked rather than
//! asserted. The real
//! [`UreqTransport`](openms::system::network_get_request::UreqTransport) is
//! given exactly two URLs in this suite, `""` and `"http://"`, and
//! `the_real_transport_refuses_every_url_this_suite_gives_it` shows that each
//! comes back as the preflight's own `MalformedUrl(<the URL verbatim>)` — a
//! value nothing inside `ureq` produces, so the refusal demonstrably happened
//! before a socket could exist. Every case that would need a peer is driven
//! through a recorded
//! [`HttpTransport`](openms::system::network_get_request::HttpTransport)
//! instead. The source's third case relies on the reserved `.invalid` TLD to
//! force a DNS failure, which is still a name lookup and still depends on the
//! host's resolver; the recorded transport reproduces the outcome it is testing
//! — `hasError` true, a non-empty error string and an empty response — without
//! one.
#![cfg(feature = "network")]

use openms::system::network_get_request::{
    DEFAULT_MAX_REDIRECTS, DEFAULT_MAX_RESPONSE_BYTES, Header, HttpTransport, MAX_REDIRECTS,
    MAX_RESPONSE_BYTES, NetworkGetRequest, RequestError, TransportError, TransportRequest,
    TransportResponse, UreqTransport,
};
use std::time::Duration;

/// A transport that answers every request from a recording.
struct Recorded(Result<TransportResponse, TransportError>);

impl HttpTransport for Recorded {
    fn get(&self, _: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        self.0.clone()
    }
}

impl Recorded {
    fn ok(status: u16, body: &[u8]) -> Self {
        Self(Ok(TransportResponse {
            status,
            headers: Vec::new(),
            body: body.to_vec(),
        }))
    }

    fn failing(error: TransportError) -> Self {
        Self(Err(error))
    }
}

/// What a [`Capturing`] transport saw.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Seen {
    url: String,
    timeout: Option<Duration>,
    max_redirects: u32,
    max_response_bytes: u64,
}

/// A transport that records what it was asked for, then answers 200 with nothing.
#[derive(Default)]
struct Capturing(std::sync::Mutex<Option<Seen>>);

impl Capturing {
    fn seen(&self) -> Seen {
        self.0
            .lock()
            .unwrap()
            .clone()
            .expect("the transport was never called")
    }
}

impl HttpTransport for Capturing {
    fn get(&self, request: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
        if let Ok(mut slot) = self.0.lock() {
            *slot = Some(Seen {
                url: request.url.to_owned(),
                timeout: request.timeout,
                max_redirects: request.max_redirects,
                max_response_bytes: request.max_response_bytes,
            });
        }
        Ok(TransportResponse::default())
    }
}

/// Class-test section `NetworkGetRequest()`.
///
/// A freshly constructed instance has no error, no error string and no response,
/// in either form — the four `TEST_EQUAL`s of the source's first section.
#[test]
fn default_construction_has_no_error_and_no_response() {
    let request = NetworkGetRequest::new();
    assert!(!request.has_error());
    assert!(request.error_string().is_empty());
    assert_eq!(request.response_text().unwrap(), "");
    assert!(request.response_binary().is_empty());
    // Native additions, consistent with the four above.
    assert_eq!(request.status(), None);
    assert!(request.error().is_none());
    assert!(request.headers().is_empty());
    assert_eq!(request.url(), "");
    assert_eq!(request.max_redirects(), DEFAULT_MAX_REDIRECTS);
    assert_eq!(request.max_response_bytes(), DEFAULT_MAX_RESPONSE_BYTES);
    // Default is also Default::default(), which the source's default ctor is.
    assert_eq!(
        NetworkGetRequest::default().response_binary(),
        request.response_binary()
    );
}

/// Class-test section `~NetworkGetRequest()`.
///
/// The source's destructor is `= default` and the section only deletes the
/// pointer allocated in the first one. Rust drops the value at the end of the
/// scope with no user-written code at all, so what is left to assert is that a
/// dropped instance released its buffer — expressed here as the value being
/// movable into a scope that ends.
#[test]
fn destruction_is_implicit() {
    let mut request = NetworkGetRequest::new();
    request.run(&Recorded::ok(200, b"payload"));
    assert_eq!(request.response_binary(), b"payload");
    drop(request);
}

/// Class-test section `void setUrl(const std::string& url)`.
///
/// `NOT_TESTABLE` upstream, "exercised through run() below". It is directly
/// observable here: the URL the transport receives is the one that was set, and
/// setting it again replaces it.
#[test]
fn set_url_is_what_the_transport_receives() {
    let transport = Capturing::default();
    let mut request = NetworkGetRequest::new();
    request.set_url("http://first.invalid/a");
    request.set_url("http://second.invalid/b");
    assert_eq!(request.url(), "http://second.invalid/b");
    request.run(&transport);
    assert_eq!(transport.seen().url, "http://second.invalid/b");
}

/// Class-test section `void setTimeout(int seconds)`.
///
/// `NOT_TESTABLE` upstream. The source calls `curl_easy_setopt(CURLOPT_TIMEOUT,
/// ...)` only when the stored value is strictly positive, so `0` and any
/// negative value leave the request without a deadline; that is what the
/// transport is handed here.
#[test]
fn set_timeout_reaches_the_transport_only_when_positive() {
    for (seconds, expected) in [
        (0, None),
        (-7, None),
        (5, Some(Duration::from_secs(5))),
        (600, Some(Duration::from_secs(600))),
    ] {
        let transport = Capturing::default();
        let mut request = NetworkGetRequest::new();
        request.set_url("http://host.invalid/");
        request.set_timeout(seconds);
        assert_eq!(request.timeout(), expected);
        request.run(&transport);
        assert_eq!(transport.seen().timeout, expected);
    }
}

/// Class-test section `void run()`, case 1: an empty URL.
///
/// The source relies on libcurl rejecting it as `CURLE_URL_MALFORMAT` with no
/// network access. This runs the real transport, whose preflight refuses the URL
/// before a socket exists, and asserts the four properties the source asserts:
/// nothing thrown, an error, a non-empty error string and an empty response in
/// both forms.
#[test]
fn run_reports_an_empty_url_as_an_error_without_panicking() {
    let mut request = NetworkGetRequest::new();
    request.set_url("");
    request.run(&UreqTransport::new());
    assert!(request.has_error());
    assert!(!request.error_string().is_empty());
    assert_eq!(request.response_text().unwrap(), "");
    assert!(request.response_binary().is_empty());
    assert_eq!(
        request.error(),
        Some(&RequestError::Transport(TransportError::MalformedUrl(
            String::new()
        )))
    );
}

/// Class-test section `void run()`, case 2: a scheme with no host.
///
/// `http://` is a URL-parse failure for libcurl and for this port's preflight,
/// so again no network is involved.
#[test]
fn run_reports_a_scheme_without_a_host_as_an_error() {
    let mut request = NetworkGetRequest::new();
    request.set_url("http://");
    request.run(&UreqTransport::new());
    assert!(request.has_error());
    assert!(!request.error_string().is_empty());
    // The preflight reports the URL verbatim; anything `ureq` itself rejected
    // would carry its own rendering instead, so this pins the offline path.
    assert_eq!(
        request.error(),
        Some(&RequestError::Transport(TransportError::MalformedUrl(
            "http://".to_owned()
        )))
    );
}

/// Class-test section `void run()`, case 3: an unresolvable host.
///
/// The source drives `CURLE_COULDNT_RESOLVE_HOST` through the reserved
/// `.invalid` TLD. That is still a resolver call; the recorded transport
/// produces the same classification offline, and the assertions are the source's
/// — no panic, an error, a non-empty error string, an empty response.
#[test]
fn run_reports_an_unresolvable_host_as_a_transport_error() {
    let mut request = NetworkGetRequest::new();
    request.set_url("http://openms-nonexistent-host.invalid/resource");
    request.set_timeout(5);
    request.run(&Recorded::failing(TransportError::HostNotFound));
    assert!(request.has_error());
    assert!(!request.error_string().is_empty());
    assert_eq!(request.response_text().unwrap(), "");
    assert_eq!(request.status(), None);
}

/// Class-test section `void run()`, case 4: reuse replaces state.
///
/// The source runs a failing request, then another on a different bad URL, and
/// requires the second to report an error of its own rather than an accumulated
/// one. The stronger form is asserted here as well: a success after a failure
/// clears the error, which the source's `has_error_ = false` at the top of `run`
/// guarantees but its test never checks.
#[test]
fn a_second_run_replaces_the_state_of_the_first() {
    let mut request = NetworkGetRequest::new();
    request.set_url("");
    request.run(&UreqTransport::new());
    assert!(request.has_error());
    assert_eq!(
        request.error(),
        Some(&RequestError::Transport(TransportError::MalformedUrl(
            String::new()
        )))
    );

    request.set_url("http://another-nonexistent-host.invalid/");
    request.set_timeout(5);
    request.run(&Recorded::failing(TransportError::HostNotFound));
    assert!(request.has_error());
    assert!(request.response_binary().is_empty());

    request.run(&Recorded::ok(200, b"recovered"));
    assert!(!request.has_error());
    assert!(request.error_string().is_empty());
    assert!(request.error().is_none());
    assert_eq!(request.response_text().unwrap(), "recovered");
    assert_eq!(request.status(), Some(200));
}

/// Class-test section `std::string getResponse() const`.
///
/// `NOT_TESTABLE` upstream. The source materialises the byte buffer as a
/// `std::string`, which cannot fail; this port checks UTF-8 instead of decoding
/// lossily, so the section is mapped by both outcomes.
#[test]
fn response_text_is_the_body_and_refuses_non_utf8() {
    let mut request = NetworkGetRequest::new();
    request.run(&Recorded::ok(200, "körper".as_bytes()));
    assert_eq!(request.response_text().unwrap(), "körper");

    request.run(&Recorded::ok(200, &[0x66, 0xff, 0x6f]));
    assert!(request.response_text().is_err());
    assert_eq!(request.response_binary(), &[0x66, 0xff, 0x6f]);
}

/// Class-test section `const std::vector<char>& getResponseBinary() const`.
///
/// `NOT_TESTABLE` upstream. The buffer holds exactly the bytes received,
/// including NULs, and is cleared at the start of the next run.
#[test]
fn response_binary_is_byte_exact_and_cleared_on_the_next_run() {
    let payload: &[u8] = &[0x00, 0x01, 0xfe, 0xff, 0x00];
    let mut request = NetworkGetRequest::new();
    request.run(&Recorded::ok(200, payload));
    assert_eq!(request.response_binary(), payload);
    request.run(&Recorded::failing(TransportError::Timeout));
    assert!(request.response_binary().is_empty());
}

/// Class-test section `bool hasError() const`.
///
/// `NOT_TESTABLE` upstream. The source's rule is `http_code >= 400`, so 399 is
/// not an error and 400 is; the boundary is asserted from both sides, and a 3xx
/// that was not followed stays a success exactly as in the source.
#[test]
fn has_error_follows_the_four_hundred_boundary() {
    let mut request = NetworkGetRequest::new();
    for status in [200, 204, 301, 302, 399] {
        request.run(&Recorded::ok(status, b"body"));
        assert!(!request.has_error(), "status {status} must not be an error");
        assert!(request.error_string().is_empty());
    }
    for status in [400, 404, 418, 500, 599] {
        request.run(&Recorded::ok(status, b"body"));
        assert!(request.has_error(), "status {status} must be an error");
        assert_eq!(request.error(), Some(&RequestError::HttpStatus(status)));
    }
}

/// Class-test section `std::string getErrorString() const`.
///
/// `NOT_TESTABLE` upstream. Empty exactly when there is no error; for an HTTP
/// status failure the text is the source's own literal, `"HTTP error " +
/// std::to_string(http_code)`, which is formed in OpenMS code and so is
/// reproduced byte for byte.
#[test]
fn error_string_is_empty_without_an_error_and_names_the_status_with_one() {
    let mut request = NetworkGetRequest::new();
    assert!(request.error_string().is_empty());
    request.run(&Recorded::ok(200, b""));
    assert!(request.error_string().is_empty());
    request.run(&Recorded::ok(404, b""));
    assert_eq!(request.error_string(), "HTTP error 404");
    request.run(&Recorded::ok(503, b""));
    assert_eq!(request.error_string(), "HTTP error 503");
    request.run(&Recorded::failing(TransportError::Timeout));
    assert!(!request.error_string().is_empty());
    assert_ne!(request.error_string(), "HTTP error 0");
}

/// Native: an error status keeps the body the peer sent.
///
/// `curl_easy_perform` returns `CURLE_OK` for a 404, so the source's write
/// callback has already stored the error page before `run` inspects the status.
/// The port sets `http_status_as_error(false)` for exactly this reason; getting
/// it wrong would silently discard every error document.
#[test]
fn an_error_status_still_carries_its_body_and_headers() {
    let mut request = NetworkGetRequest::new();
    request.run(&Recorded(Ok(TransportResponse {
        status: 404,
        headers: vec![Header {
            name: "content-type".into(),
            value: b"text/plain".to_vec(),
        }],
        body: b"no such tool".to_vec(),
    })));
    assert!(request.has_error());
    assert_eq!(request.response_text().unwrap(), "no such tool");
    assert_eq!(request.status(), Some(404));
    assert_eq!(
        request
            .header("Content-Type")
            .unwrap()
            .value_text()
            .unwrap(),
        "text/plain"
    );
    assert!(request.header("x-absent").is_none());
}

/// Native: a transport failure leaves no body, where the source can leave a truncated one.
///
/// The source's write callback appends as the bytes arrive, so a transfer that
/// dies mid-body — an expiring `CURLOPT_TIMEOUT`, a reset connection — leaves
/// what did arrive in `response_bytes_` and `getResponseBinary()` returns it
/// next to `hasError() == true`. `HttpTransport::get` returns a whole response
/// or a `TransportError` and never both, so here the body is empty. This is the
/// documented divergence; the test exists so it cannot drift unnoticed in either
/// direction.
#[test]
fn a_transport_failure_leaves_no_body_at_all() {
    let mut request = NetworkGetRequest::new();
    request.run(&Recorded::ok(200, b"the first run's body"));
    assert_eq!(request.response_binary(), b"the first run's body");

    // Exactly the conditions under which the source would hold a partial body.
    for error in [
        TransportError::Timeout,
        TransportError::Io("connection reset by peer".into()),
    ] {
        request.run(&Recorded::failing(error));
        assert!(request.has_error());
        assert!(request.response_binary().is_empty());
        assert_eq!(request.response_text().unwrap(), "");
        assert_eq!(request.status(), None);
        assert!(request.headers().is_empty());
    }
}

/// Native: the bounds a caller may set, and the ceilings on them.
#[test]
fn limits_are_checked_against_their_ceilings() {
    let transport = Capturing::default();
    let mut request = NetworkGetRequest::new();
    assert!(request.set_max_redirects(MAX_REDIRECTS + 1).is_err());
    assert!(
        request
            .set_max_response_bytes(MAX_RESPONSE_BYTES + 1)
            .is_err()
    );
    // A refused limit leaves the previous one in place.
    assert_eq!(request.max_redirects(), DEFAULT_MAX_REDIRECTS);
    assert_eq!(request.max_response_bytes(), DEFAULT_MAX_RESPONSE_BYTES);

    request.set_max_redirects(3).unwrap();
    request.set_max_response_bytes(1024).unwrap();
    request.set_url("http://host.invalid/");
    request.run(&transport);
    let seen = transport.seen();
    assert_eq!(seen.max_redirects, 3);
    assert_eq!(seen.max_response_bytes, 1024);
}

/// Native: every transport failure renders as a non-empty, distinct message.
///
/// The source forwards `curl_easy_strerror`, whose text is libcurl's; this port
/// does not reproduce those strings, so what is asserted is the contract the
/// class test actually relies on — an error always has something to say.
#[test]
fn every_transport_error_renders_non_empty() {
    let errors = [
        TransportError::MalformedUrl("http://".into()),
        TransportError::UnsupportedScheme("file:///x".into()),
        TransportError::HostNotFound,
        TransportError::ConnectionFailed("http://h/".into()),
        TransportError::Timeout,
        TransportError::TooManyRedirects(10),
        TransportError::ResponseTooLarge(1024),
        TransportError::ResponseHeaderTooLarge,
        TransportError::TooManyHeaders(1024),
        TransportError::Tls("certificate".into()),
        TransportError::Protocol("bad chunk".into()),
        TransportError::Io("broken pipe".into()),
        TransportError::Other("unclassified".into()),
    ];
    let mut rendered = Vec::new();
    for error in errors {
        let mut request = NetworkGetRequest::new();
        request.run(&Recorded::failing(error.clone()));
        let text = request.error_string();
        assert!(!text.is_empty());
        assert_eq!(text, error.to_string());
        rendered.push(text);
    }
    rendered.sort();
    rendered.dedup();
    assert_eq!(rendered.len(), 13);
}

/// `SYSTEM/CurlInit.h`: no counterpart, and nothing to test.
///
/// `CurlInit::ensure` exists only to run `curl_global_init` once per process
/// before any curl handle is created. `ureq` has no global state to initialise,
/// so the header has no Rust type at all — see `docs/CURL_INIT_SUPPORT.md`.
/// What is checkable is the property the guard exists to provide: the transport
/// is usable immediately, from several threads, with no setup call.
#[test]
fn no_global_initialisation_is_needed_before_a_transport_is_used() {
    let transport = UreqTransport::new();
    let results: Vec<Option<RequestError>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..4)
            .map(|_| {
                let transport = &transport;
                scope.spawn(move || {
                    let mut request = NetworkGetRequest::new();
                    request.set_url("");
                    request.run(transport);
                    request.error().cloned()
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("worker thread panicked"))
            .collect()
    });
    let refused = RequestError::Transport(TransportError::MalformedUrl(String::new()));
    assert_eq!(results, vec![Some(refused); 4]);
}

/// Native: every use of the *real* transport in this suite is refused offline.
///
/// The claim that no test in this group reaches the network rests on two facts:
/// `UreqTransport` is handed only the URLs listed here, and `check_http_url` —
/// the first statement of its `get`, and a pure function of the string — refuses
/// each of them before `ureq` is called at all. The second fact is what this
/// test pins: the error is the preflight's own
/// `MalformedUrl(<the URL verbatim>)`, which nothing inside `ureq` produces, so
/// a change that let one of these URLs through to a socket would fail here
/// rather than quietly start resolving names on a build machine.
#[test]
fn the_real_transport_refuses_every_url_this_suite_gives_it() {
    for url in ["", "http://"] {
        let mut request = NetworkGetRequest::new();
        request.set_url(url);
        request.run(&UreqTransport::new());
        assert_eq!(
            request.error(),
            Some(&RequestError::Transport(TransportError::MalformedUrl(
                url.to_owned()
            ))),
            "{url:?} must be refused by the preflight, not by the network"
        );
        assert!(request.response_binary().is_empty());
        assert_eq!(request.status(), None);
    }
}

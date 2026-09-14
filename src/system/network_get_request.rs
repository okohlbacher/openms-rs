// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Synchronous HTTP GET of the Core SDK `SYSTEM/NetworkGetRequest.h`, together
//! with the `SYSTEM/CurlInit.h` global guard it depends on.
//!
//! The source drives libcurl directly: `CurlInit::ensure` runs
//! `curl_global_init` once behind a function-local static, `run` configures a
//! `CURL*` easy handle with `CURLOPT_FOLLOWLOCATION`, `CURLOPT_NOSIGNAL` and an
//! optional `CURLOPT_TIMEOUT`, and reports failure through
//! `hasError`/`getErrorString` instead of throwing. This port keeps that
//! observable contract and replaces the transport with a trait,
//! [`HttpTransport`](crate::system::network_get_request::HttpTransport), whose
//! only shipped implementation —
//! [`UreqTransport`](crate::system::network_get_request::UreqTransport) — is
//! built on `ureq`. The whole module is behind the non-default `network`
//! feature, because that is the feature which pulls in `ureq`.
//!
//! The trait exists for one concrete reason: a test that reaches the internet is
//! flaky, slow and a privacy problem in a scientific SDK, so the crate's own
//! tests drive recorded responses through the trait and never open a socket. How
//! that was checked rather than assumed — which entry points can open one, which
//! URLs the real transport is given, and which test pins the refusal to the
//! offline preflight — is the *No test reaches the network* section of
//! `docs/NETWORK_GET_REQUEST_SUPPORT.md`.
//!
//! `CurlInit` has no counterpart at all — `ureq` needs no global
//! initialisation, so there is nothing for an RAII guard to own. See
//! `docs/NETWORK_GET_REQUEST_SUPPORT.md` and `docs/CURL_INIT_SUPPORT.md`.
//!
//! # Behaviour worth knowing before calling
//!
//! * A body is captured even when the status is an error. `curl_easy_perform`
//!   succeeds for a 404, so the source's write callback has already stored the
//!   error page by the time `run` inspects the status code. This port sets
//!   `ureq`'s `http_status_as_error` to `false` for exactly that reason, so
//!   [`NetworkGetRequest::response_binary`](crate::system::network_get_request::NetworkGetRequest::response_binary)
//!   still holds the server's explanation while
//!   [`has_error`](crate::system::network_get_request::NetworkGetRequest::has_error)
//!   is `true`.
//! * Redirects are followed, but a bounded number of times. The source sets
//!   `CURLOPT_FOLLOWLOCATION` and never sets `CURLOPT_MAXREDIRS`, so its ceiling
//!   is whatever the libcurl it links chose — a value that is not knowable from
//!   the SDK, which ships no libcurl. This port pins its own,
//!   [`DEFAULT_MAX_REDIRECTS`](crate::system::network_get_request::DEFAULT_MAX_REDIRECTS),
//!   and reports [`TransportError::TooManyRedirects`](crate::system::network_get_request::TransportError::TooManyRedirects)
//!   beyond it.
//! * The response body is bounded. The source streams into an unbounded
//!   `std::vector<char>`; this port refuses beyond
//!   [`DEFAULT_MAX_RESPONSE_BYTES`](crate::system::network_get_request::DEFAULT_MAX_RESPONSE_BYTES).
//! * A transfer that fails part-way discards what had arrived. The source's
//!   write callback appends as the bytes come in, so a transfer that dies
//!   mid-body leaves those bytes in `getResponseBinary()` while `hasError()` is
//!   `true`; here
//!   [`response_binary`](crate::system::network_get_request::NetworkGetRequest::response_binary)
//!   is empty whenever the transport failed. See *Divergences from the source*
//!   below.
//!
//! # Divergences from the source
//!
//! Four differences are visible to a caller and are not configuration this
//! module can undo. They are stated here rather than left to be discovered, and
//! `docs/NETWORK_GET_REQUEST_SUPPORT.md` carries the same list with the reasons.
//!
//! 1. **A partial body is not kept.** Input: a transfer that fails after the
//!    peer has sent part of the body — a `CURLOPT_TIMEOUT` that expires mid-body,
//!    or a connection reset. Source: `curl_easy_perform` returns non-`CURLE_OK`,
//!    but the write callback has already appended what arrived, so
//!    `getResponseBinary()` returns a truncated body next to `hasError() ==
//!    true`. Here: [`HttpTransport::get`](crate::system::network_get_request::HttpTransport::get)
//!    yields either a whole response or a [`TransportError`](crate::system::network_get_request::TransportError),
//!    never both, so a failed run has no body at all. Reproducing the source
//!    would mean handing every caller a buffer of unknowable truncation, which
//!    the source's own callers — `Network::downloadFile` and `UpdateCheck::run`,
//!    both of which test `hasError()` first — never look at.
//! 2. **A `Content-Encoding` response body is decompressed.** Input: a response
//!    carrying `Content-Encoding: gzip`. Source: libcurl decompresses only when
//!    `CURLOPT_ACCEPT_ENCODING` was set, and this source never sets it, so the
//!    compressed bytes reach `getResponseBinary()` verbatim. Here: `ureq`
//!    decompresses on the strength of the response header alone. The request
//!    never advertises `gzip` — see [`UreqTransport`](crate::system::network_get_request::UreqTransport)
//!    — so only a server compressing unasked reaches this, but the decompression
//!    itself cannot be switched off without dropping `ureq`'s `gzip` feature.
//! 3. **TLS trust anchors are compiled in, not the operating system's.** Input:
//!    an `https` URL whose certificate chains to a CA the machine trusts but the
//!    Mozilla root program does not — an enterprise inspection proxy, an
//!    institutional CA. Source: libcurl verifies against the platform store.
//!    Here: `ureq`'s default `rustls` backend verifies against a compiled-in
//!    copy of the Mozilla roots, so that request fails with
//!    [`TransportError::Tls`](crate::system::network_get_request::TransportError::Tls)
//!    — and, the other way round, a CA an administrator distrusted locally is
//!    still trusted. `SSL_CERT_FILE` and `CURL_CA_BUNDLE` are not read either.
//! 4. **A SOCKS proxy in the environment is ignored.** `HTTP_PROXY`,
//!    `HTTPS_PROXY`, `ALL_PROXY` and `NO_PROXY` are honoured by both, but only
//!    for an HTTP proxy; libcurl also speaks `socks5://`, and this build of
//!    `ureq` does not.
//!
//! Items 2 to 4 are properties of the transport crate and its enabled features,
//! so changing them would mean changing the crate's dependencies. Item 1 is a
//! deliberate contract choice. A caller needing any of the four can implement
//! [`HttpTransport`](crate::system::network_get_request::HttpTransport) itself,
//! which is the seam the source does not have.

use crate::{Error, Result};
use std::fmt;
use std::time::Duration;

/// Redirects followed before a request is refused, unless the caller changes it.
///
/// The source sets `CURLOPT_FOLLOWLOCATION` to `1` and never sets
/// `CURLOPT_MAXREDIRS`, so how long a redirect chain it will walk is decided by
/// the libcurl it was linked against rather than by OpenMS. That value is not
/// determinable from the SDK — the pinned checkout contains no libcurl source or
/// header — and it has not been constant across libcurl's own history, so this
/// port makes no claim about it and states only what OpenMS configures: nothing.
///
/// What this port does is fixed and knowable: following a chain is unbounded
/// work driven by untrusted input, so the chain is bounded here. Ten is `ureq`'s
/// own default and is far past anything a well-behaved service needs.
pub const DEFAULT_MAX_REDIRECTS: u32 = 10;

/// Largest redirect chain a caller may request.
pub const MAX_REDIRECTS: u32 = 64;

/// Response body accepted by default, in bytes.
///
/// The source appends to a `std::vector<char>` with no ceiling, so the peer
/// decides how much memory the process commits.
pub const DEFAULT_MAX_RESPONSE_BYTES: u64 = 64 * 1024 * 1024;

/// Largest response body a caller may ask for, in bytes.
pub const MAX_RESPONSE_BYTES: u64 = 4 * 1024 * 1024 * 1024;

/// Response headers retained before a response is refused.
///
/// Headers are a native addition — the source exposes none — so this exists to
/// keep the addition bounded rather than to match anything upstream.
pub const MAX_HEADERS: usize = 1024;

/// One response header, with its value kept as raw bytes.
///
/// `http` guarantees the name is lowercase US-ASCII. The value is *not*
/// guaranteed to be UTF-8, so it is stored verbatim and decoded on request by
/// [`Header::value_text`]; discarding or lossily decoding a header the source
/// never looks at would be a silent loss for no gain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Header {
    /// Lowercase field name.
    pub name: String,
    /// Field value, exactly as received.
    pub value: Vec<u8>,
}

impl Header {
    /// The value decoded as UTF-8.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the bytes are not valid UTF-8. RFC
    /// 9110 restricts field values to US-ASCII plus opaque bytes, so a value
    /// that fails here is one no caller should interpret as text.
    pub fn value_text(&self) -> Result<&str> {
        std::str::from_utf8(&self.value)
            .map_err(|_| Error::InvalidValue("response header value is not UTF-8".into()))
    }
}

/// What a [`HttpTransport`] was asked to fetch.
///
/// The lifetime borrows the URL from the caller's
/// [`NetworkGetRequest`], which owns it for the duration of the call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TransportRequest<'a> {
    /// Absolute URL to request.
    pub url: &'a str,
    /// Whole-request deadline, or `None` for no deadline.
    ///
    /// This is the source's `CURLOPT_TIMEOUT`, which covers the entire transfer
    /// rather than any single phase of it.
    pub timeout: Option<Duration>,
    /// Redirects to follow before failing.
    pub max_redirects: u32,
    /// Response body accepted before failing, in bytes.
    pub max_response_bytes: u64,
}

/// A response a [`HttpTransport`] delivered.
///
/// A non-2xx status is *not* an error at this level, exactly as
/// `curl_easy_perform` returning `CURLE_OK` is not: the status check happens one
/// layer up, in [`NetworkGetRequest::run`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TransportResponse {
    /// HTTP status code of the final response in the redirect chain.
    pub status: u16,
    /// Response headers, in the order the transport reported them.
    pub headers: Vec<Header>,
    /// Response body bytes.
    pub body: Vec<u8>,
}

/// Why a transport could not produce a response.
///
/// These are the conditions the source reports through `curl_easy_strerror`
/// after a non-`CURLE_OK` `curl_easy_perform`. The variants name the libcurl
/// code each corresponds to; the port does not reproduce libcurl's message
/// *text*, which is an internal string of that library and which the class test
/// only ever checks for non-emptiness.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TransportError {
    /// The URL could not be parsed (`CURLE_URL_MALFORMAT`).
    MalformedUrl(String),
    /// The URL names a scheme this transport does not speak
    /// (`CURLE_UNSUPPORTED_PROTOCOL`).
    UnsupportedScheme(String),
    /// The host name did not resolve (`CURLE_COULDNT_RESOLVE_HOST`).
    HostNotFound,
    /// No connection could be established (`CURLE_COULDNT_CONNECT`).
    ConnectionFailed(String),
    /// The whole-request deadline elapsed (`CURLE_OPERATION_TIMEDOUT`).
    Timeout,
    /// The redirect chain exceeded the configured limit
    /// (`CURLE_TOO_MANY_REDIRECTS`). Carries the limit that was hit.
    TooManyRedirects(u32),
    /// The body exceeded the configured limit (`CURLE_FILESIZE_EXCEEDED`).
    /// Carries the limit, not the actual size, which was never read.
    ResponseTooLarge(u64),
    /// The response head exceeded what the transport will buffer.
    ResponseHeaderTooLarge,
    /// The response carried more headers than [`MAX_HEADERS`].
    TooManyHeaders(usize),
    /// TLS negotiation or certificate verification failed
    /// (`CURLE_SSL_CONNECT_ERROR` / `CURLE_PEER_FAILED_VERIFICATION`).
    Tls(String),
    /// The peer spoke malformed HTTP (`CURLE_WEIRD_SERVER_REPLY`).
    Protocol(String),
    /// A socket-level read or write failed (`CURLE_RECV_ERROR` /
    /// `CURLE_SEND_ERROR`).
    Io(String),
    /// Any other transport failure.
    Other(String),
}

impl fmt::Display for TransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MalformedUrl(detail) => write!(f, "malformed URL: {detail}"),
            Self::UnsupportedScheme(detail) => write!(f, "unsupported protocol: {detail}"),
            Self::HostNotFound => f.write_str("could not resolve host name"),
            Self::ConnectionFailed(detail) => write!(f, "could not connect to server: {detail}"),
            Self::Timeout => f.write_str("timeout was reached"),
            Self::TooManyRedirects(limit) => {
                write!(f, "number of redirects hit maximum amount ({limit})")
            }
            Self::ResponseTooLarge(limit) => {
                write!(f, "response body exceeds {limit} bytes")
            }
            Self::ResponseHeaderTooLarge => f.write_str("response header block is too large"),
            Self::TooManyHeaders(limit) => write!(f, "response carries more than {limit} headers"),
            Self::Tls(detail) => write!(f, "TLS error: {detail}"),
            Self::Protocol(detail) => write!(f, "protocol error: {detail}"),
            Self::Io(detail) => write!(f, "transfer error: {detail}"),
            Self::Other(detail) => write!(f, "transport error: {detail}"),
        }
    }
}

/// Why the last [`NetworkGetRequest::run`] did not produce a usable response.
///
/// The split reproduces the source's two branches exactly: `curl_easy_perform`
/// returning something other than `CURLE_OK` becomes
/// [`RequestError::Transport`], and an otherwise successful transfer whose
/// status is `>= 400` becomes [`RequestError::HttpStatus`]. The source's third
/// branch — `"Failed to initialize libcurl"` — has no counterpart, because
/// there is no global library state to fail to initialise.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RequestError {
    /// The transport could not complete the transfer.
    Transport(TransportError),
    /// The transfer completed with an HTTP status of 400 or above.
    HttpStatus(u16),
}

impl fmt::Display for RequestError {
    /// `HttpStatus` renders as the source's own literal, `HTTP error N`.
    ///
    /// That string is formed in OpenMS code rather than inside libcurl, so it is
    /// reproduced byte for byte. The transport texts are this port's own; see
    /// [`TransportError`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Transport(error) => error.fmt(f),
            Self::HttpStatus(status) => write!(f, "HTTP error {status}"),
        }
    }
}

/// A blocking HTTP GET transport.
///
/// The source has no such seam: `NetworkGetRequest::run` calls libcurl
/// directly. Splitting it out is the one structural difference in this group,
/// and it is deliberate — it is what lets the crate's tests exercise status
/// classification, body capture and the update-check decision against recorded
/// responses without a network, and it lets a caller supply a transport with
/// its own proxy or certificate policy.
///
/// An implementation must be usable from several threads at once, because
/// nothing in this module serialises access to it.
pub trait HttpTransport: Sync {
    /// Perform one GET and return the final response.
    ///
    /// Redirects are followed internally, up to `request.max_redirects`. A
    /// non-2xx status is a successful return, not an error.
    ///
    /// # Errors
    ///
    /// Returns the [`TransportError`] describing why no response could be
    /// produced.
    fn get(
        &self,
        request: &TransportRequest<'_>,
    ) -> std::result::Result<TransportResponse, TransportError>;
}

/// Synchronous HTTP GET request, and the state its last run left behind.
///
/// Mirrors the source class: set a URL, optionally a timeout, call
/// [`run`](Self::run), then read the response or the error. The instance is
/// reusable and each run replaces — never accumulates — the previous state,
/// which is the behaviour the class test's fourth case pins down.
///
/// ```
/// use openms::system::network_get_request::{
///     HttpTransport, NetworkGetRequest, TransportRequest, TransportResponse, TransportError,
/// };
///
/// struct Recorded;
/// impl HttpTransport for Recorded {
///     fn get(&self, _: &TransportRequest<'_>) -> Result<TransportResponse, TransportError> {
///         Ok(TransportResponse { status: 404, headers: Vec::new(), body: b"gone".to_vec() })
///     }
/// }
///
/// let mut request = NetworkGetRequest::new();
/// request.set_url("http://example.invalid/missing");
/// request.set_timeout(5);
/// request.run(&Recorded);
///
/// // A 404 is an error, and the body the server sent is still available.
/// assert!(request.has_error());
/// assert_eq!(request.error_string(), "HTTP error 404");
/// assert_eq!(request.response_binary(), b"gone");
/// ```
#[derive(Clone, Debug)]
pub struct NetworkGetRequest {
    url: String,
    timeout_seconds: i32,
    max_redirects: u32,
    max_response_bytes: u64,
    response_bytes: Vec<u8>,
    headers: Vec<Header>,
    status: Option<u16>,
    error: Option<RequestError>,
}

impl Default for NetworkGetRequest {
    /// The source's default constructor: no URL, no timeout, no response, no error.
    fn default() -> Self {
        Self::new()
    }
}

impl NetworkGetRequest {
    /// A request with no URL, no timeout and no observed state.
    ///
    /// The source's members start as an empty URL, `timeout_ = 0`,
    /// `has_error_ = false` and empty response and error strings; the class
    /// test's first section asserts all four.
    pub fn new() -> Self {
        Self {
            url: String::new(),
            timeout_seconds: 0,
            max_redirects: DEFAULT_MAX_REDIRECTS,
            max_response_bytes: DEFAULT_MAX_RESPONSE_BYTES,
            response_bytes: Vec::new(),
            headers: Vec::new(),
            status: None,
            error: None,
        }
    }

    /// Set the URL the next [`run`](Self::run) will request.
    ///
    /// No validation happens here, as in the source: an empty or malformed URL
    /// is reported by `run` as a transport error rather than refused now.
    pub fn set_url(&mut self, url: impl Into<String>) {
        self.url = url.into();
    }

    /// The URL the next [`run`](Self::run) will request.
    ///
    /// Native addition; the source keeps `url_` private with no accessor.
    pub fn url(&self) -> &str {
        &self.url
    }

    /// Set the whole-request timeout in seconds.
    ///
    /// `0` — the default — leaves the timeout unset, so the request blocks until
    /// the peer answers or the transport gives up on its own. The source tests
    /// `timeout_ > 0` before calling `curl_easy_setopt(CURLOPT_TIMEOUT, ...)`,
    /// so a negative value is also "unset"; that is reproduced rather than
    /// rejected, and [`timeout`](Self::timeout) reports the resulting
    /// [`Duration`].
    pub fn set_timeout(&mut self, seconds: i32) {
        self.timeout_seconds = seconds;
    }

    /// The configured timeout, or `None` when it is unset.
    ///
    /// Native addition. `Some` exactly when the source would have called
    /// `curl_easy_setopt` with `CURLOPT_TIMEOUT`, that is when the stored
    /// seconds are strictly positive.
    pub fn timeout(&self) -> Option<Duration> {
        if self.timeout_seconds > 0 {
            Some(Duration::from_secs(u64::from(
                self.timeout_seconds.unsigned_abs(),
            )))
        } else {
            None
        }
    }

    /// Set how many redirects [`run`](Self::run) may follow.
    ///
    /// Native addition; the source never sets `CURLOPT_MAXREDIRS` and so leaves
    /// the ceiling to its libcurl. `0` stops the transport at the first redirect
    /// and returns it as the response.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `limit` exceeds [`MAX_REDIRECTS`].
    pub fn set_max_redirects(&mut self, limit: u32) -> Result<()> {
        if limit > MAX_REDIRECTS {
            return Err(Error::InvalidValue("redirect limit exceeded".into()));
        }
        self.max_redirects = limit;
        Ok(())
    }

    /// How many redirects [`run`](Self::run) may follow.
    pub fn max_redirects(&self) -> u32 {
        self.max_redirects
    }

    /// Set the largest response body [`run`](Self::run) will accept, in bytes.
    ///
    /// Native addition; the source has no ceiling at all.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `limit` exceeds
    /// [`MAX_RESPONSE_BYTES`].
    pub fn set_max_response_bytes(&mut self, limit: u64) -> Result<()> {
        if limit > MAX_RESPONSE_BYTES {
            return Err(Error::InvalidValue("response size limit exceeded".into()));
        }
        self.max_response_bytes = limit;
        Ok(())
    }

    /// The largest response body [`run`](Self::run) will accept, in bytes.
    pub fn max_response_bytes(&self) -> u64 {
        self.max_response_bytes
    }

    /// Execute the GET synchronously through `transport`.
    ///
    /// The previous response, headers, status and error are cleared first, so a
    /// reused instance reports the outcome of this call alone. This never
    /// panics and never returns an error to the caller: failure is observable
    /// through [`has_error`](Self::has_error) and
    /// [`error_string`](Self::error_string), which is the source's contract and
    /// the one thing the class test checks in every case.
    ///
    /// The error state is set when the transport fails, or when the transfer
    /// completed with an HTTP status of 400 or above — the source's `>= 400`
    /// test, reproduced including its consequence that a 3xx response which was
    /// not followed is *not* an error. A body received alongside an error status
    /// is kept, because the source's write callback has already stored it.
    ///
    /// A body received before a *transport* failure is **not** kept, and this is
    /// the one place where the observable contract differs from the source. A
    /// transfer that dies mid-body — a deadline that expires while the body is
    /// still arriving, a reset connection — leaves the bytes that did arrive in
    /// the source's `response_bytes_`, so `getResponseBinary()` returns a
    /// truncated body while `hasError()` is `true`. Here the transport reports
    /// either a response or a [`TransportError`], so
    /// [`response_binary`](Self::response_binary) is empty after any transport
    /// failure. The difference matters only to a caller that reads the body
    /// without checking [`has_error`](Self::has_error) first, which neither
    /// in-tree caller of the source does.
    pub fn run(&mut self, transport: &dyn HttpTransport) {
        self.response_bytes.clear();
        self.headers.clear();
        self.status = None;
        self.error = None;

        let request = TransportRequest {
            url: &self.url,
            timeout: self.timeout(),
            max_redirects: self.max_redirects,
            max_response_bytes: self.max_response_bytes,
        };
        match transport.get(&request) {
            Ok(response) => {
                self.status = Some(response.status);
                self.headers = response.headers;
                self.response_bytes = response.body;
                if response.status >= 400 {
                    self.error = Some(RequestError::HttpStatus(response.status));
                }
            }
            Err(error) => self.error = Some(RequestError::Transport(error)),
        }
    }

    /// Response body as text.
    ///
    /// The source's `getResponse` copies the bytes into a `std::string`, which
    /// is a byte container and so can never fail. A Rust `str` is UTF-8, so this
    /// checks instead of decoding lossily; use
    /// [`response_binary`](Self::response_binary) for the bytes themselves.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the body is not valid UTF-8.
    pub fn response_text(&self) -> Result<&str> {
        std::str::from_utf8(&self.response_bytes)
            .map_err(|_| Error::InvalidValue("response body is not UTF-8".into()))
    }

    /// Raw response body, valid until the next [`run`](Self::run).
    ///
    /// The counterpart of `getResponseBinary`, whose `std::vector<char>` becomes
    /// a byte slice. Empty before the first run, and empty after a run whose
    /// transport failed — where the source can hold a truncated body; see
    /// [`run`](Self::run).
    pub fn response_binary(&self) -> &[u8] {
        &self.response_bytes
    }

    /// Response headers of the last run.
    ///
    /// Native addition: the source stores none, because its write callback only
    /// receives the body. Empty before the first run and after a transport
    /// failure.
    pub fn headers(&self) -> &[Header] {
        &self.headers
    }

    /// The first header with this name, matched case-insensitively.
    ///
    /// Native addition. Names arriving from the transport are lowercase, so a
    /// caller may pass any casing.
    pub fn header(&self, name: &str) -> Option<&Header> {
        self.headers
            .iter()
            .find(|header| header.name.eq_ignore_ascii_case(name))
    }

    /// HTTP status of the last run, or `None` when no response was received.
    ///
    /// Native addition. The source reads `CURLINFO_RESPONSE_CODE` into a local
    /// and keeps only the `>= 400` verdict, so a caller cannot tell 200 from
    /// 204, or 404 from 500 except through the error string.
    pub fn status(&self) -> Option<u16> {
        self.status
    }

    /// Whether the last [`run`](Self::run) produced an error.
    pub fn has_error(&self) -> bool {
        self.error.is_some()
    }

    /// Structured description of the last error, or `None` when there was none.
    ///
    /// Native addition: the source keeps only the rendered string, so a caller
    /// wanting to branch on the failure has to parse English.
    pub fn error(&self) -> Option<&RequestError> {
        self.error.as_ref()
    }

    /// Human-readable description of the last error.
    ///
    /// Empty exactly when [`has_error`](Self::has_error) is `false`, which is
    /// the source's documented invariant. For an HTTP status failure the text is
    /// the source's own `HTTP error N`.
    pub fn error_string(&self) -> String {
        self.error
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_default()
    }
}

/// [`HttpTransport`] backed by `ureq`.
///
/// Configured to match what the source asks of libcurl, option by option:
///
/// | libcurl option | value in the source | this transport |
/// |---|---|---|
/// | `CURLOPT_FOLLOWLOCATION` | `1` | `max_redirects` from the request, default [`DEFAULT_MAX_REDIRECTS`] |
/// | `CURLOPT_MAXREDIRS` | never set, so its libcurl's choice | bounded here, and exceeding it is an error |
/// | `CURLOPT_TIMEOUT` | the request's seconds when `> 0` | `timeout_global`, same meaning |
/// | `CURLOPT_NOSIGNAL` | `1` | no counterpart needed; `ureq` never installs a signal handler |
/// | `CURLOPT_USERAGENT` | never set, so no `User-Agent` is sent | `user_agent("")`, which suppresses `ureq`'s own |
/// | `CURLOPT_ACCEPT_ENCODING` | never set, so no `Accept-Encoding` is sent | `accept_encoding("")`, which suppresses `ureq`'s `gzip` |
/// | `Accept` | libcurl's default `*/*` | `ureq`'s default `*/*` |
/// | status handling | body captured, status checked afterwards | `http_status_as_error(false)`, status checked afterwards |
///
/// The two suppressed headers matter because the source's request is what the
/// OpenMS REST server sees: `ureq` would otherwise announce itself in
/// `User-Agent` and offer `gzip` in `Accept-Encoding`, neither of which libcurl
/// sends unless asked.
///
/// # What this transport cannot match
///
/// `ureq` speaks HTTP and HTTPS only. libcurl is normally built with the `FILE`
/// protocol as well, which the upstream `Network` class test relies on; that
/// difference is recorded in `docs/NETWORK_SUPPORT.md` and surfaces here as
/// [`TransportError::UnsupportedScheme`].
///
/// Three more differences follow from the transport crate and its enabled
/// features rather than from anything this module configures, and are set out in
/// full in the module documentation:
///
/// * a response carrying `Content-Encoding: gzip` is decompressed here and is
///   not by the source, even though neither request asks for compression;
/// * TLS certificates are verified against a compiled-in copy of the Mozilla
///   root program, where libcurl uses the platform's trust store and honours
///   `SSL_CERT_FILE` / `CURL_CA_BUNDLE`;
/// * a `socks5://` proxy named in the environment is used by libcurl and
///   ignored here, while an HTTP proxy from `HTTP_PROXY` / `HTTPS_PROXY` /
///   `ALL_PROXY` and the exceptions in `NO_PROXY` are honoured by both.
#[derive(Clone, Copy, Debug, Default)]
pub struct UreqTransport {
    _private: (),
}

impl UreqTransport {
    /// A transport with `ureq`'s default agent configuration.
    pub fn new() -> Self {
        Self { _private: () }
    }
}

impl HttpTransport for UreqTransport {
    fn get(
        &self,
        request: &TransportRequest<'_>,
    ) -> std::result::Result<TransportResponse, TransportError> {
        check_http_url(request.url)?;
        let mut response = ureq::get(request.url)
            .config()
            // The source reads the status only after a successful transfer, so
            // the body of a 4xx/5xx is already captured. Letting ureq turn the
            // status into an error would discard it.
            .http_status_as_error(false)
            // The source sets neither CURLOPT_USERAGENT nor
            // CURLOPT_ACCEPT_ENCODING, so libcurl sends neither header. An
            // empty value tells ureq to send neither either, which keeps the
            // request the REST server sees the one the source would have made.
            .user_agent("")
            .accept_encoding("")
            .max_redirects(request.max_redirects)
            .max_redirects_will_error(true)
            .timeout_global(request.timeout)
            .build()
            .call()
            .map_err(|error| map_ureq_error(&error, request))?;

        let status = response.status().as_u16();
        let header_map = response.headers();
        if header_map.len() > MAX_HEADERS {
            return Err(TransportError::TooManyHeaders(MAX_HEADERS));
        }
        let mut headers = Vec::with_capacity(header_map.len());
        for (name, value) in header_map {
            headers.push(Header {
                name: name.as_str().to_owned(),
                value: value.as_bytes().to_vec(),
            });
        }
        let body = response
            .body_mut()
            .with_config()
            .limit(request.max_response_bytes)
            .read_to_vec()
            .map_err(|error| map_ureq_error(&error, request))?;
        Ok(TransportResponse {
            status,
            headers,
            body,
        })
    }
}

/// Reject a URL this transport cannot speak, before any socket is opened.
///
/// libcurl parses the URL and picks a protocol handler before it connects, so
/// `CURLE_URL_MALFORMAT` and `CURLE_UNSUPPORTED_PROTOCOL` are offline failures
/// there. `ureq` reaches the same conclusions, but only once the request is on
/// its way, so this reproduces the early rejection: a URL that is not
/// `http`/`https` with a non-empty host never reaches the network. That keeps the
/// two malformed-URL cases of the class test genuinely offline, and it is what
/// makes `file://` — which libcurl normally does handle, and which the upstream
/// `Network` class test uses — report
/// [`TransportError::UnsupportedScheme`] rather than a connection failure.
///
/// The scheme is split off and checked here, because `http::Uri` refuses
/// `file:///x` as malformed while this check must still call it an unsupported
/// scheme. The host comes from `ureq::http::Uri`, the parser `ureq` runs on the
/// URL it then requests, so the preflight and the request agree on the host. A
/// URL that parser refuses, or one whose authority names no host
/// (`http://:8080/x`, `http://:[:]`), is [`TransportError::MalformedUrl`].
/// Everything further is `ureq`'s business.
fn check_http_url(url: &str) -> std::result::Result<(), TransportError> {
    let Some((scheme, _)) = url.split_once("://") else {
        return Err(TransportError::MalformedUrl(url.to_owned()));
    };
    if !scheme.eq_ignore_ascii_case("http") && !scheme.eq_ignore_ascii_case("https") {
        return Err(TransportError::UnsupportedScheme(url.to_owned()));
    }
    match ureq::http::Uri::try_from(url) {
        Ok(uri) if uri.host().is_some_and(|host| !host.is_empty()) => Ok(()),
        _ => Err(TransportError::MalformedUrl(url.to_owned())),
    }
}

/// Classify a `ureq` failure into the libcurl condition the source would report.
///
/// `ureq::Error` is `#[non_exhaustive]` and several of its variants exist only
/// under features this crate does not enable, so the catch-all arm is required
/// rather than lazy.
fn map_ureq_error(error: &ureq::Error, request: &TransportRequest<'_>) -> TransportError {
    match error {
        ureq::Error::BadUri(detail) => TransportError::MalformedUrl(detail.clone()),
        ureq::Error::Http(detail) => TransportError::MalformedUrl(detail.to_string()),
        ureq::Error::RequireHttpsOnly(url) => TransportError::UnsupportedScheme(url.clone()),
        ureq::Error::HostNotFound => TransportError::HostNotFound,
        ureq::Error::ConnectionFailed => TransportError::ConnectionFailed(request.url.to_owned()),
        ureq::Error::ConnectProxyFailed(detail) => TransportError::ConnectionFailed(detail.clone()),
        ureq::Error::Timeout(_) => TransportError::Timeout,
        ureq::Error::TooManyRedirects => TransportError::TooManyRedirects(request.max_redirects),
        ureq::Error::BodyExceedsLimit(limit) => TransportError::ResponseTooLarge(*limit),
        ureq::Error::LargeResponseHeader(_, _) => TransportError::ResponseHeaderTooLarge,
        ureq::Error::Tls(detail) => TransportError::Tls((*detail).to_owned()),
        ureq::Error::TlsRequired => TransportError::Tls("TLS required by the transport".to_owned()),
        ureq::Error::Protocol(detail) => TransportError::Protocol(detail.to_string()),
        ureq::Error::Io(detail) => TransportError::Io(detail.to_string()),
        other => TransportError::Other(other.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Fixed(std::result::Result<TransportResponse, TransportError>);

    impl HttpTransport for Fixed {
        fn get(
            &self,
            _: &TransportRequest<'_>,
        ) -> std::result::Result<TransportResponse, TransportError> {
            self.0.clone()
        }
    }

    fn ok(status: u16, body: &[u8]) -> Fixed {
        Fixed(Ok(TransportResponse {
            status,
            headers: Vec::new(),
            body: body.to_vec(),
        }))
    }

    #[test]
    fn fresh_instance_has_no_state() {
        let request = NetworkGetRequest::new();
        assert!(!request.has_error());
        assert!(request.error_string().is_empty());
        assert_eq!(request.response_text().unwrap(), "");
        assert!(request.response_binary().is_empty());
        assert_eq!(request.status(), None);
        assert_eq!(request.timeout(), None);
    }

    #[test]
    fn timeout_is_unset_unless_strictly_positive() {
        let mut request = NetworkGetRequest::new();
        request.set_timeout(-1);
        assert_eq!(request.timeout(), None);
        request.set_timeout(0);
        assert_eq!(request.timeout(), None);
        request.set_timeout(5);
        assert_eq!(request.timeout(), Some(Duration::from_secs(5)));
    }

    #[test]
    fn status_at_or_above_four_hundred_is_an_error_with_the_body_kept() {
        let mut request = NetworkGetRequest::new();
        request.run(&ok(399, b"redirect page"));
        assert!(!request.has_error());
        request.run(&ok(400, b"bad request"));
        assert!(request.has_error());
        assert_eq!(request.error_string(), "HTTP error 400");
        assert_eq!(request.response_binary(), b"bad request");
    }

    #[test]
    fn a_later_run_replaces_the_earlier_state() {
        let mut request = NetworkGetRequest::new();
        request.run(&ok(500, b"boom"));
        assert!(request.has_error());
        request.run(&ok(200, b"fine"));
        assert!(!request.has_error());
        assert!(request.error_string().is_empty());
        assert_eq!(request.response_text().unwrap(), "fine");
    }

    #[test]
    fn a_transport_failure_leaves_no_response() {
        let mut request = NetworkGetRequest::new();
        request.run(&ok(200, b"fine"));
        request.run(&Fixed(Err(TransportError::HostNotFound)));
        assert!(request.has_error());
        assert!(request.response_binary().is_empty());
        assert_eq!(request.status(), None);
        assert!(!request.error_string().is_empty());
    }

    #[test]
    fn limits_are_refused_above_their_ceilings() {
        let mut request = NetworkGetRequest::new();
        assert!(request.set_max_redirects(MAX_REDIRECTS).is_ok());
        assert!(request.set_max_redirects(MAX_REDIRECTS + 1).is_err());
        assert_eq!(request.max_redirects(), MAX_REDIRECTS);
        assert!(request.set_max_response_bytes(MAX_RESPONSE_BYTES).is_ok());
        assert!(
            request
                .set_max_response_bytes(MAX_RESPONSE_BYTES + 1)
                .is_err()
        );
    }

    #[test]
    fn the_url_preflight_refuses_before_any_socket() {
        assert_eq!(
            check_http_url(""),
            Err(TransportError::MalformedUrl(String::new()))
        );
        assert_eq!(
            check_http_url("http://"),
            Err(TransportError::MalformedUrl("http://".into()))
        );
        assert_eq!(
            check_http_url("http://:8080/x"),
            Err(TransportError::MalformedUrl("http://:8080/x".into()))
        );
        assert!(matches!(
            check_http_url("file:///tmp/x.txt"),
            Err(TransportError::UnsupportedScheme(_))
        ));
        assert!(matches!(
            check_http_url("ftp://host/x"),
            Err(TransportError::UnsupportedScheme(_))
        ));
        assert_eq!(check_http_url("http://host/a?b#c"), Ok(()));
        assert_eq!(check_http_url("HTTPS://host"), Ok(()));
        assert_eq!(check_http_url("http://user:pw@host:81/a"), Ok(()));
        assert_eq!(check_http_url("http://[::1]:81/a"), Ok(()));
    }

    #[test]
    fn a_malformed_url_fails_offline_through_the_real_transport() {
        // Both cases are refused by the preflight, so no socket, no DNS lookup
        // and no dependence on the machine's connectivity. The error carries the
        // URL verbatim, which is `check_http_url`'s signature and not anything
        // ureq produces, so this pins *where* the refusal happened.
        for url in ["", "http://"] {
            let mut request = NetworkGetRequest::new();
            request.set_url(url);
            request.run(&UreqTransport::new());
            assert!(request.has_error());
            assert!(!request.error_string().is_empty());
            assert!(request.response_binary().is_empty());
            assert_eq!(
                request.error(),
                Some(&RequestError::Transport(TransportError::MalformedUrl(
                    url.to_owned()
                )))
            );
        }
    }

    #[test]
    fn an_empty_host_fails_offline_through_the_real_transport() {
        // A port and no host. The preflight refuses it, so ureq never hands an
        // empty name to the resolver; the verbatim URL in the error shows the
        // refusal was the preflight's.
        let url = "http://:8080/x";
        let mut request = NetworkGetRequest::new();
        request.set_url(url);
        request.run(&UreqTransport::new());
        assert!(request.has_error());
        assert!(request.response_binary().is_empty());
        assert_eq!(request.status(), None);
        assert_eq!(
            request.error(),
            Some(&RequestError::Transport(TransportError::MalformedUrl(
                url.to_owned()
            )))
        );
    }

    #[test]
    fn the_host_is_the_one_the_http_uri_parser_finds() {
        let mut wrong = Vec::new();
        for url in [
            // `http::Uri` finds an empty host in these, where the earlier
            // hand-written strip found a non-empty one and let them through to
            // the resolver.
            "HTTPS://:[:]",
            "http://user@:[:]/x",
            "http://:[::1]:80/",
            // Not URIs. `ureq` refused these offline as well, but only after
            // the preflight had passed them.
            "http://ho st/",
            "http://host/a b",
            "http://[::1/",
        ] {
            let got = check_http_url(url);
            if got != Err(TransportError::MalformedUrl(url.to_owned())) {
                wrong.push(format!("{url}: {got:?}"));
            }
        }
        // Empty userinfo, an empty port and a mixed-case scheme still name a host.
        for url in ["http://@host/", "http://host:/", "HtTp://host?q"] {
            let got = check_http_url(url);
            if got != Ok(()) {
                wrong.push(format!("{url}: {got:?}"));
            }
        }
        assert!(wrong.is_empty(), "{wrong:#?}");
    }

    #[test]
    fn non_utf8_bodies_and_headers_are_refused_rather_than_decoded() {
        let mut request = NetworkGetRequest::new();
        request.run(&Fixed(Ok(TransportResponse {
            status: 200,
            headers: vec![Header {
                name: "x-note".into(),
                value: vec![0xff],
            }],
            body: vec![0xff, 0xfe],
        })));
        assert!(request.response_text().is_err());
        assert!(request.header("X-Note").unwrap().value_text().is_err());
        assert_eq!(request.response_binary(), &[0xff, 0xfe]);
    }
}

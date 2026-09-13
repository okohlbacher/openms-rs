# Synchronous HTTP GET

`system::network_get_request` covers `SYSTEM/NetworkGetRequest.h` and
`SYSTEM/NetworkGetRequest.cpp` at SDK `bc9cc12`. The module is behind the
non-default **`network`** feature, which is what pulls in `ureq`; a build without
it has no `system::network_get_request` at all, and the crate's default feature
set does not enable it.

The header is small and its contract is unusually explicit: `run` never throws,
failures are observable through `hasError`/`getErrorString`, redirects are
followed, and the instance is reusable with each run replacing the previous
state. All four are preserved. What changed is underneath: the source calls
libcurl directly, and this port calls a trait.

## API mapping

Every public member of the header and its `.cpp`, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class NetworkGetRequest` | `network_get_request::NetworkGetRequest` | |
| `NetworkGetRequest()` | `NetworkGetRequest::new` / `Default` | Empty URL, `timeout_ = 0`, no error, no response. |
| `~NetworkGetRequest()` | — | `= default` on both sides; Rust drops the `Vec`s. |
| `void setUrl(const std::string&)` | `NetworkGetRequest::set_url` | No validation here, as in the source. |
| `void setTimeout(int seconds)` | `NetworkGetRequest::set_timeout` | `i32`, keeping the source's "unset unless `> 0`" rule. |
| `void run()` | `NetworkGetRequest::run(&dyn HttpTransport)` | Takes the transport; see *Native differences*. |
| `std::string getResponse() const` | `NetworkGetRequest::response_text() -> Result<&str>` | Checked UTF-8 instead of a byte-container `std::string`. |
| `const std::vector<char>& getResponseBinary() const` | `NetworkGetRequest::response_binary() -> &[u8]` | Byte-for-byte the same buffer. |
| `bool hasError() const` | `NetworkGetRequest::has_error` | |
| `std::string getErrorString() const` | `NetworkGetRequest::error_string` | Empty exactly when there is no error. |
| `NetworkGetRequest(const NetworkGetRequest&) = delete` | — | The type is `Clone`; see *Native differences*. |
| `operator=(const NetworkGetRequest&) = delete` | — | Same. |
| `std::vector<char> response_bytes_` (private) | `response_bytes` (private) | |
| `std::string url_` (private) | `url` (private), read by `NetworkGetRequest::url` | |
| `int timeout_ = 0` (private) | `timeout_seconds` (private), read by `NetworkGetRequest::timeout` | |
| `bool has_error_ = false` (private) | `error: Option<RequestError>` | One field replaces the flag and the string. |
| `std::string error_string_` (private) | `error: Option<RequestError>` | |
| `writeCallback` (file-static, `.cpp`) | — | Not ported: the transport returns a whole body, so there is no incremental sink. |

Native additions: `HttpTransport`, `TransportRequest`, `TransportResponse`,
`TransportError`, `RequestError`, `UreqTransport`, `Header`,
`NetworkGetRequest::{url, timeout, status, error, headers, header,
set_max_redirects, max_redirects, set_max_response_bytes, max_response_bytes}`,
and the constants `DEFAULT_MAX_REDIRECTS`, `MAX_REDIRECTS`,
`DEFAULT_MAX_RESPONSE_BYTES`, `MAX_RESPONSE_BYTES`, `MAX_HEADERS`.

`SYSTEM/CurlInit.h`, which this header's `.cpp` depends on, has its own
document: [CURL_INIT_SUPPORT.md](CURL_INIT_SUPPORT.md). It is not ported.

## Preserved source conventions

**`run` never fails loudly.** It returns nothing, cannot panic, and reports
every failure through `has_error`/`error_string`. The class test wraps all four
of its cases in a `try`/`catch` and asserts nothing was thrown; the Rust
signature makes that unrepresentable.

**The error boundary is `status >= 400`, not "not 2xx".** A 3xx that was not
followed, a 204, a 100 — none is an error. `tests/network_get_request.rs::has_error_follows_the_four_hundred_boundary`
asserts 399 and 400 from both sides.

**An error status keeps its body.** `curl_easy_perform` returns `CURLE_OK` for a
404, so the source's write callback has already stored the error page when `run`
looks at the status. This is why `UreqTransport` sets
`http_status_as_error(false)` and checks the status itself: the obvious
translation, letting `ureq` turn a 4xx into an error, would silently discard
every error document the server sent.

**`"HTTP error N"` is byte-exact.** That string is built in OpenMS code —
`"HTTP error " + std::to_string(http_code)` — so `RequestError::HttpStatus`
renders exactly it.

**A timeout of `0` or less is no timeout.** The source guards the
`curl_easy_setopt(CURLOPT_TIMEOUT, ...)` call with `timeout_ > 0`, so a negative
value is quietly treated as unset rather than rejected. `timeout()` returns
`None` for both.

**Each run replaces the previous state.** `run` clears the response, headers,
status and error before doing anything, exactly as the source's three `clear()`
calls do.

**Reuse is supported.** The instance is not consumed and may be run any number
of times with different URLs.

## Native differences

**The transport is injected.** The source's `run` is hard-wired to libcurl. Here
it takes a `&dyn HttpTransport`, with `UreqTransport` as the shipped
implementation. Two reasons, in order of importance:

1. **No test may reach the internet.** A scientific SDK whose test suite makes
   outbound requests is flaky on a build farm, slow behind a proxy, and leaks the
   fact that it ran. Every test in `tests/network_get_request.rs`,
   `tests/network.rs` and `tests/update_check.rs` either uses a recorded
   transport or a URL that the real transport refuses before opening a socket.
2. A caller with its own proxy, certificate pinning or offline mirror can supply
   one, which the source cannot express at all.

The trait requires `Sync`, so one transport can serve several threads; the source
relies on libcurl's per-handle isolation for the same property.

**libcurl option by option.** What the source configures, and what replaces it:

| libcurl option | source | port |
|---|---|---|
| `CURLOPT_URL` | the stored URL | same |
| `CURLOPT_WRITEFUNCTION` / `CURLOPT_WRITEDATA` | appends to `response_bytes_` | the transport returns the body |
| `CURLOPT_FOLLOWLOCATION` | `1` | `max_redirects` from the request |
| `CURLOPT_MAXREDIRS` | unset, so unlimited | `DEFAULT_MAX_REDIRECTS` = 10, ceiling `MAX_REDIRECTS` = 64 |
| `CURLOPT_NOSIGNAL` | `1`, for a thread-safe DNS timeout | no counterpart: `ureq` installs no signal handler, so there is nothing to disable |
| `CURLOPT_TIMEOUT` | the stored seconds when `> 0` | `timeout_global`, the same whole-request meaning |
| `CURLINFO_RESPONSE_CODE` | read, compared to 400, then discarded | kept, and exposed through `status()` |

**The redirect chain is bounded.** libcurl's default `CURLOPT_MAXREDIRS` is
unlimited, so a server that redirects to itself makes the source loop until the
timeout — or forever, when no timeout was set. Following redirects is work
scaled by untrusted input, so the port caps it and reports
`TransportError::TooManyRedirects`.

**The response body is bounded.** The source's write callback appends to a
`std::vector<char>` with no ceiling: the peer decides how much memory the process
commits. The port refuses past `DEFAULT_MAX_RESPONSE_BYTES` (64 MiB), adjustable
up to `MAX_RESPONSE_BYTES` (4 GiB).

**A URL preflight runs before any socket.** libcurl parses the URL and selects a
protocol handler before connecting, so `CURLE_URL_MALFORMAT` and
`CURLE_UNSUPPORTED_PROTOCOL` are offline outcomes there. `ureq` reaches the same
conclusions later, so `UreqTransport` reproduces the early rejection: a URL
without `://`, with an empty host, or with a scheme other than `http`/`https` is
refused without touching the network. This is what keeps the class test's first
two `run()` cases genuinely offline, and it is where `file://` — which libcurl
normally does handle — surfaces as `TransportError::UnsupportedScheme`.

**Error classification is structured; libcurl's message text is not
reproduced.** The source forwards `curl_easy_strerror(res)` verbatim. Those
strings belong to libcurl, vary between its versions, and are not something a
caller can depend on; the class test only ever asserts they are non-empty. So
`TransportError` names the *condition* — each variant's rustdoc gives the
`CURLE_*` code it corresponds to — and renders its own text. The one string that
is OpenMS's own, `"HTTP error N"`, is reproduced exactly.

The source's third error string, `"Failed to initialize libcurl"`, has no
counterpart: `curl_easy_init` can fail, and building a `ureq` request cannot.

**`getResponse` is checked.** A `std::string` is a byte container, so the
source's accessor cannot fail; a Rust `&str` is UTF-8, so `response_text` returns
`Result` and refuses invalid bytes rather than substituting replacement
characters. `response_binary` is always available and always exact.

**Headers are exposed, as raw bytes.** The source captures none — its write
callback only receives the body — so this is an addition. Values are kept as
`Vec<u8>` and decoded on demand by `Header::value_text`, because an HTTP field
value is not required to be UTF-8 and neither discarding nor lossily decoding one
is acceptable. Bounded by `MAX_HEADERS` = 1024.

**Copying is allowed.** The source deletes the copy constructor and assignment
operator, because copying a `CURL*`-adjacent object is a trap. There is no handle
here — the state is three owned buffers — so `Clone` is derived and copying is
merely a copy.

**Failure is not thread-local state.** The source stores `has_error_` and
`error_string_` as two members that can, in principle, disagree; the port stores
one `Option<RequestError>` from which both accessors derive, so
`error_string().is_empty() == !has_error()` holds by construction.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| Redirects followed | 10 by default, ≤ 64 | unlimited |
| Response body | 64 MiB by default, ≤ 4 GiB | unbounded |
| Response headers | ≤ 1024 | not captured at all |
| Response header block | `ureq`'s own limit, reported as `ResponseHeaderTooLarge` | unbounded |
| URL scheme | `http` / `https` only, checked offline | libcurl's full protocol set |

All **nine** `START_SECTION`s of `NetworkGetRequest_test.cpp` are mapped in
`tests/network_get_request.rs`; the four sub-cases of the `run()` section have a
test each. Five of the nine are `NOT_TESTABLE` upstream, deferred to `run()`;
each is nevertheless given a direct test here, because the trait makes the
behaviour observable that libcurl hid.

| Section | Rust test | A value it reproduces |
|---|---|---|
| `NetworkGetRequest()` | `default_construction_has_no_error_and_no_response` | `getErrorString().empty() == true` |
| `~NetworkGetRequest()` | `destruction_is_implicit` | the section only deletes the pointer; nothing is asserted upstream |
| `void setUrl(...)` | `set_url_is_what_the_transport_receives` | the URL reaching the transport is the last one set |
| `void setTimeout(int)` | `set_timeout_reaches_the_transport_only_when_positive` | `0` and `-7` leave the request without a deadline |
| `void run()` case 1 | `run_reports_an_empty_url_as_an_error_without_panicking` | `hasError() == true`, response empty |
| `void run()` case 2 | `run_reports_a_scheme_without_a_host_as_an_error` | `getErrorString().size() != 0` for `"http://"` |
| `void run()` case 3 | `run_reports_an_unresolvable_host_as_a_transport_error` | `hasError() == true`, `getResponse().empty() == true` |
| `void run()` case 4 | `a_second_run_replaces_the_state_of_the_first` | the second run reports its own error, not the first's |
| `std::string getResponse()` | `response_text_is_the_body_and_refuses_non_utf8` | the body round-trips as text |
| `getResponseBinary()` | `response_binary_is_byte_exact_and_cleared_on_the_next_run` | a body containing NUL survives byte for byte |
| `bool hasError()` | `has_error_follows_the_four_hundred_boundary` | 399 is not an error; 400 is |
| `std::string getErrorString()` | `error_string_is_empty_without_an_error_and_names_the_status_with_one` | `"HTTP error 404"` |

**Evidence tier 3** for the transcribed class-test expectations (the four
default-state assertions, and the "no throw, error set, non-empty string, empty
response" shape of each `run()` case). **Tier 4** for everything derived rather
than transcribed: the 399/400 boundary, which follows from `http_code >= 400` in
the `.cpp` rather than from any test literal; the exact text `"HTTP error 404"`,
which follows from the source's own concatenation; the option-by-option
configuration table; and every native bound. No C++ was executed and no retained
C++ output exists for this header, so no tier 1 or tier 2 claim is made.

The source's case 3 asserts against a real DNS failure for the reserved
`.invalid` TLD. That is a resolver call, and its outcome depends on the host —
a captive portal or a wildcard resolver changes it — so it is reproduced through
a recorded `TransportError::HostNotFound` instead. The assertions are the
source's.

The source carries no `#pragma omp`, so there is no parallel behaviour to match
and no performance gap to record.

Source hashes, line anchors and the class-test review are in
[the provenance record](../tests/data/network_provenance.json).

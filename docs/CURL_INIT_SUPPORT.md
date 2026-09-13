# libcurl global initialisation

`SYSTEM/CurlInit.h` and `SYSTEM/CurlInit.cpp` at SDK `bc9cc12`, 28 + 29 lines.
**This header has no Rust counterpart, deliberately.** It is recorded here
rather than silently omitted, because a reader comparing the two trees will
otherwise wonder whether it was missed.

## API mapping

Every public member of the header, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class CurlInit` | — | **Not ported.** See below. |
| `static void CurlInit::ensure()` | — | **Not ported**: nothing to initialise. |
| `CurlInit::CurlInit()` (private) | — | Not ported; it exists only to call `curl_global_init(CURL_GLOBAL_DEFAULT)`. |
| `CurlInit::~CurlInit()` (private) | — | Not ported; it exists only to call `curl_global_cleanup()`. |
| `CurlInit(const CurlInit&) = delete` | — | Not ported; Rust has no implicit copy to delete. |
| `CurlInit& operator=(const CurlInit&) = delete` | — | Not ported. |

## Preserved source conventions

None survive, because none can. What the header *achieves* — that a caller never
has to think about global library state before making a request — is preserved:
`system::network_get_request::UreqTransport::new` needs no prior call, and
`tests/network_get_request.rs::no_global_initialisation_is_needed_before_a_transport_is_used`
drives four threads through it concurrently to show that.

## Native differences

**The whole class is the difference.** libcurl requires `curl_global_init`
before the first easy handle in a process, is documented as not thread-safe when
called concurrently with other library use, and must be paired with
`curl_global_cleanup`. `CurlInit::ensure` solves that with a function-local
`static CurlInit instance`, which C++11 guarantees is constructed exactly once
even under concurrent first calls, and whose destructor runs at process exit.

`ureq` has no such requirement: an `Agent` owns its own connection pool and TLS
configuration, there is no process-wide state to install, and there is nothing to
tear down. A Rust type holding a place for this would be an empty struct with an
empty method, which is worse than an honest absence — it would suggest to a
caller that the call matters.

There is one behaviour of the C++ that is therefore not reproduced: the source
calls `curl_global_cleanup()` at static-destruction time. Nothing in the port
needs it, but a mixed process that also links libcurl through some other library
is unaffected either way, since `CurlInit` only ever managed OpenMS's own
initialisation.

## Checked boundaries and evidence

No boundaries: there is no code.

There is no `CurlInit_test.cpp` in the SDK, so no class-test section is
unaccounted for. The only consumer of the header inside the Core SDK is
`NetworkGetRequest.cpp`, which calls `CurlInit::ensure()` as the first statement
of `run`; the Rust `run` has no equivalent first statement, and that is the whole
of the difference.

Source hashes are in [the provenance record](../tests/data/network_provenance.json).

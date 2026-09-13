# One-shot URL download

`system::network` covers `SYSTEM/Network.h` and `SYSTEM/Network.cpp` at SDK
`bc9cc12`. Behind the non-default **`network`** feature, as the rest of this
group.

The header is one static method with an unusually detailed Doxygen block — the
filename derivation, the `.0`/`.1`/`.2` collision rule, the `""` → `"./"`
mapping, the three failure paths and the fact that a partial file is not cleaned
up are all stated there rather than left to the code. Every one of those
statements is either preserved or explicitly changed below.

## API mapping

Every public member of the header and its `.cpp`, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class Network` | the module itself | A stateless namespace of one static method becomes a free function. |
| `static void downloadFile(const std::string& url, const std::string& download_folder)` | `network::download_file(url, folder) -> Result<PathBuf>` | Uses `UreqTransport`; returns the path written. |
| — | `network::download_file_with(transport, url, folder) -> Result<PathBuf>` | Native: the injectable form the tests use. |
| `saveFileName_` (file-static, `.cpp`) | `network::save_file_name(url, dest_folder) -> Result<String>` | Public here; its rules are observable on disk. |

Native additions: `DOWNLOAD_TIMEOUT_SECONDS`, `FALLBACK_BASENAME`,
`MAX_URL_BYTES`, `MAX_NAME_SUFFIX`, and `download_file_with`.

## Preserved source conventions

**The ten-minute timeout is the source's.** `setTimeout(600)`, documented in the
header as a 10-minute timeout, is `DOWNLOAD_TIMEOUT_SECONDS` and is asserted at
the transport in `tests/network.rs::download_file_writes_the_fixture_under_its_url_basename`.

**The filename derivation is rule for rule**: the query is cut at the first `?`,
*then* the fragment at the first `#` — in that order, so a `#` inside a query
string is already gone by the time the fragment is looked for — then the basename
of what remains is taken, and an empty basename becomes the literal `"download"`.

**Percent escapes are not decoded**, by either side, so `%20` reaches the file
name as three characters.

**Existing files are never overwritten.** `.0`, `.1`, `.2`, … are tried in order
until one is free. `tests/network.rs::repeated_downloads_never_overwrite`
downloads the same URL three times and gets `run.tsv`, `run.tsv.0`, `run.tsv.1`.

**An empty destination folder means the current directory.** The source maps
`""` to `"./"`.

**An HTTP status of 400 or above is a download failure**, because
`NetworkGetRequest` classifies it as one and `downloadFile` only asks
`hasError()`.

**An empty response body still produces a file.** The source writes
`data.data()` for `data.size()` bytes with no emptiness test.

**The failure message is the source's**: `Download of '<url>' failed!. Error:
<error>` — including the `!.`, which is a typo upstream and is kept so the two
messages match.

## Native differences

**The transport is injected and `file://` is not supported.** `download_file`
constructs a `UreqTransport`; `download_file_with` takes any
`HttpTransport`. The upstream `Network_test.cpp` downloads
`file://<test data>/Network_test_fixture.txt`, because libcurl is normally built
with the `FILE` protocol. `ureq` speaks HTTP and HTTPS only, so that URL reports
`TransportError::UnsupportedScheme`. Adding a `FILE` transport would mean
implementing a protocol handler this crate has no other use for, and — more to
the point — the class test's purpose is to exercise the *naming and writing*
path, not the protocol. `tests/network.rs` therefore serves the upstream
fixture's bytes, unchanged and hashed, through a recorded transport, and asserts
more than the upstream section does: the derived name, the file's existence, and
byte equality with the fixture.

**The written path is returned.** The source returns `void` and emits two
`OPENMS_LOG_INFO` lines, the second of which is the path. The `system` module may
not depend on the crate's logging module — the module-dependency ratchet forbids
that edge — so returning the path is both the available option and the more
useful one: a caller that wants to log it can, and a caller that wants to open it
no longer has to re-derive the name.

**A partially written file is removed.** The header states that a partial file is
*not* cleaned up on write failure. This port removes it, so a failed download
leaves the directory as it found it. A removal that itself fails does not mask
the write error.

**The basename is split on `/` only, on every platform.** The source uses
`std::filesystem::path::filename`, which additionally treats `\` as a separator
on Windows: `http://h/a\b.txt` names `b.txt` there and `a\b.txt` on Linux. `/` is
the path separator of a URL whatever the host filesystem does, so this port
derives `a\b.txt` everywhere — and then refuses it, because a name containing a
separator is not a name.

**`.`, `..`, a path separator and NUL are refused.** The source concatenates the
derived basename into a path unchecked. `.` and `..` reach `std::ofstream` and
fail there with an `IOException`, which is a worse diagnosis of the same problem;
a separator is a silent write outside the folder on the platform that accepts it.

**The suffix search terminates.** The source's `while (fs::exists(...)) ++i;` has
no upper bound, so a directory already holding `x`, `x.0`, `x.1`, … probes the
filesystem forever. `MAX_NAME_SUFFIX` = 10 000 stops it with an error.

**The URL length is bounded** at `MAX_URL_BYTES` = 64 KiB, checked before
anything is allocated or a socket opened. The source bounds nothing.

**Path joining, not string concatenation.** The source builds `folder + "/" +
name`, so an empty folder yields `".//name"`. `Path::join` yields `"./name"`,
which names the same file.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| URL length | ≤ 64 KiB | unbounded |
| Collision suffixes tried | ≤ 10 000 | unbounded loop |
| Response body | inherited from `NetworkGetRequest`: 64 MiB | unbounded |
| Derived basename | no `/`, `\`, NUL, `.` or `..` | unchecked |

The **one** `START_SECTION` of `Network_test.cpp` is mapped by
`tests/network.rs::download_file_writes_the_fixture_under_its_url_basename`,
which reproduces its single assertion — the downloaded file exists at
`<folder>/Network_test_fixture.txt` — and adds byte equality against the same
fixture. The fixture, `tests/data/network_test_fixture.txt`, is the upstream
`src/tests/class_tests/openms/data/Network_test_fixture.txt` copied unchanged;
its sha256 `539a0ccf…` is recorded in the provenance manifest and matches the
upstream file.

**Evidence tier 3** for that section and for the naming rules, which are
transcribed from `saveFileName_`. **Tier 4** for the derived expectations: the
`run.tsv` → `run.tsv.0` → `run.tsv.1` sequence follows from the suffix rule
rather than from any upstream literal, as does `https://host` deriving `host`
(the last `/`-separated component of that string), and every native bound. No
C++ was executed; the fixture is upstream input data, not retained C++ output, so
this is not a tier 1 claim.

The source carries no `#pragma omp`; it is serial and so is the port.

Source hashes, line anchors and the fixture hash are in
[the provenance record](../tests/data/network_provenance.json).

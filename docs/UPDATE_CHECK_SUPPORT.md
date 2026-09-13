# Update check against the OpenMS REST server

`system::update_check` covers `SYSTEM/UpdateCheck.h` and
`SYSTEM/UpdateCheck.cpp` at SDK `bc9cc12`. Behind the non-default **`network`**
feature.

The header is two lines of API and the `.cpp` is 120 lines of policy: a
rate-limiting stamp file, a platform-tagged identifier, a five-second query and a
version comparison. The policy is what matters, and it is reproduced exactly —
including two defects, which are marked as such.

## API mapping

Every public member of the header, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class UpdateCheck` | the module itself | A stateless namespace of one static method. |
| `static void run(const std::string& tool_name, const std::string& version, int debug_level)` | `update_check::run(context, transport, running_version, tool_name, version, debug_level) -> Result<UpdateCheckReport>` | Four things the source reaches through globals become parameters; see *Native differences*. |

The `.cpp` reaches into three other headers. Their relevant members:

| C++ member used by `UpdateCheck.cpp` | Rust counterpart | Notes |
|---|---|---|
| `File::getOpenMSConfigDir()` | `system::file::FileContext::get_openms_config_dir` | Already ported. |
| `File::exists` / `File::readable` | `system::file::exists` / `readable` | Already ported. |
| `VersionInfo::getVersionStruct()` | the `running_version` parameter | See *Native differences*. |
| `VersionInfo::getRevision()` | — | **Not ported**: the source computes `revision` and never uses it. |
| `VersionInfo::VersionDetails` | `update_check::SemanticVersion` | The part `UpdateCheck` depends on; `CONCEPT/VersionInfo.h` is not otherwise ported. |
| `VersionDetails::create` | `SemanticVersion::parse` | Returns `Option` where the source returns the `EMPTY` sentinel. |
| `VersionDetails::operator<` | `SemanticVersion::is_less_than` | Not `PartialOrd`; see below. |
| `VersionDetails::operator==` / `operator!=` | `PartialEq` / `Eq` | Field-wise, including the pre-release string. |
| `VersionDetails::operator>` | `SemanticVersion::is_greater_than` | `!(< || ==)`, inheriting the caveat. |
| `VersionDetails::EMPTY` | `SemanticVersion::default` + `SemanticVersion::is_empty` | |
| `VersionDetails::version_major/minor/patch/pre_release_identifier` | `SemanticVersion::major/minor/patch/pre_release` | |
| `StringUtils::toInt32` | private `to_int32` in this module | Exact conversion rules; not otherwise ported. |
| `stat` + `utime` | private `stamp_now` | `std::fs::FileTimes`. |

Native additions: `SemanticVersion`, `UpdateCheckOutcome`, `UpdateCheckReport`,
`platform_tag`, `architecture_tag`, `tool_version_string`, `query_url`, `is_due`,
`should_report_update`, `version_file_path`, and the constants
`UPDATE_URL_PREFIX`, `QUERY_TIMEOUT_SECONDS`, `MIN_QUERY_INTERVAL`,
`VERSION_FILE_EXTENSION`, `DISABLE_ENVIRONMENT_VARIABLE`, `MAX_TOOL_NAME_BYTES`,
`MAX_VERSION_BYTES`, `MAX_QUERY_RESPONSE_BYTES`.

## Preserved source conventions

**The identifier is byte-exact.**
`"OpenMS" + "_" + "Default_" + platform + "_" + architecture + "_" + tool_name + "_" + version`,
which the source's own comment illustrates as
`OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0`.

**The platform tags are the source's preprocessor order**: `Win`, then `Mac`,
then `Linux`, then `Unix`, then `unknown`. Note this differs from
`system::build_info`, which folds every non-macOS Unix into `Linux`; the two
headers really do disagree, and both are reproduced as they are.

**The architecture tag is `"32"` when a pointer is four bytes and `"64"`
otherwise** — everything that is not four bytes is `"64"`, exactly as the
source's ternary says.

**The endpoint is the source's, over plain HTTP**:
`http://openms-update.cs.uni-tuebingen.de/check/`. No `https` substitution: the
server that answers is the one the source names, and quietly changing the scheme
would be a change in behaviour disguised as a fix.

**The query timeout is five seconds**, and the interval is 24 hours.

**The stamp is bumped before the request.** The source sets the modification
time, then queries. A server that is down therefore still consumes the day's
budget rather than causing a retry on every tool invocation.
`tests/update_check.rs::a_failed_query_still_consumes_the_window` pins that
ordering.

**The access time is preserved** when the modification time is bumped; the
source reads `st_atime` and passes it back to `utime`, and `FileTimes` leaves an
unset field alone.

**A `stat` or `utime` failure warns and continues**; it does not abort the check.

**The three privacy notices appear only at `debug_level > 0`**, and the
announcement of an available update appears regardless. Their text is
transcribed, including the trailing space at the end of the second line and of
`"Connecting to REST server successful. "`.

**A response that does not parse, or parses as `0.0.0`, is ignored.**
`VersionDetails::create` returns `EMPTY` for a parse failure, and the source
skips the comparison when the result equals `EMPTY`; a literal `0.0.0` is
indistinguishable and takes the same branch. `should_report_update` reproduces
both.

**`VersionDetails::create`'s parse rules are exact**, down to
`StringUtils::toInt32`: leading and trailing space/tab/newline/carriage-return
are ignored, one leading `+` is allowed, a leading `-` gives a negative number,
and every remaining character must be consumed. So `"1.2-alpha"` fails — the
minor component would be `"2-alpha"` — while `"1.2.3-alpha"` succeeds, and the
`-` is only ever looked for *after* the second `.`.

**`operator<` is reproduced with its inconsistency.** A pre-release sorts below
the same triple without one, but two *different* pre-releases on one triple are
neither less nor greater while still comparing unequal. That is not a total
order, so `PartialOrd` is deliberately **not** implemented — deriving one would
quietly change which versions trigger an update. `system::stop_watch` omits
`PartialOrd` for the same class of reason.

## Native differences

**Four globals become parameters.** The source reads its inputs from process
state: `File::getOpenMSConfigDir()`, libcurl, `VersionInfo::getVersionStruct()`
and — implicitly — the wall clock. The port takes a `&FileContext`, a
`&dyn HttpTransport` and a `&SemanticVersion`. This is what makes the class
test's own scenario expressible without `setenv`: the upstream test redirects
`XDG_CONFIG_HOME` to a path whose parent is a regular file, mutating the
process environment and restoring it afterwards;
`tests/update_check.rs::an_uncreatable_config_directory_warns_and_makes_no_request`
builds the same situation as a value.

**Reporting replaces logging.** The source writes to `OPENMS_LOG_WARN` and
`OPENMS_LOG_INFO` and returns `void`. `system` may not depend on the crate's
logging module — the module-dependency ratchet forbids the edge — so `run`
returns an `UpdateCheckReport` whose `warnings` and `notices` carry the same
strings in the same order, and whose `outcome` names the branch taken. Seven
outcomes replace a function whose only observable effect was text:
`ConfigDirectoryUnavailable`, `StampFileUnusable`, `NotDue`, `QueryFailed`,
`ServerVersionUnusable`, `UpToDate`, `UpdateAvailable`. Each is reachable in a
test; the source's branches are not.

**`run` returns `Result` for argument errors only.** Every runtime failure the
source tolerates — an uncreatable config directory, a failed stat, a failed
request — is an outcome, not an error, because the source guarantees the tool
runs on regardless. What *is* an error is an argument the source should have
rejected: an empty tool name, one longer than 256 bytes, one starting with `.`,
or one containing `/`, `\` or NUL. The source concatenates the name straight into
a path, so `../../evil` writes outside the config directory.

**`VersionDetails::create` returns `Option`.** The source's `EMPTY` sentinel is
also the legitimate value `0.0.0`, so its callers cannot tell a parse failure
from a zero version. `parse` returns `None` for failure and `Some` for success,
and `is_empty` recovers the conflation where `UpdateCheck` actually depends on
it — so the *decision* is unchanged while the parser becomes honest.

**The stamp is touched through a writable handle.** The source calls `stat` then
`utime` on the path. `std::fs::File::set_times` needs an open handle, and on
Windows that handle must be writable, so the port opens the file for writing. The
failure modes differ slightly — a file readable but not writable warns here where
the source might have succeeded — and both warn and continue.

**The response body is bounded at 64 KiB.** The answer is a version string. The
source accepts whatever the server sends into an unbounded buffer.

**A non-UTF-8 response is `ServerVersionUnusable`.** The source hands the raw
bytes to `create`, where any non-numeric byte fails the integer conversion and
yields `EMPTY` — the same branch. Nothing is decoded lossily.

**`OPENMS_DISABLE_UPDATE_CHECK` is named but not read**, here as in the Core SDK
— see *Defects reproduced*. `DISABLE_ENVIRONMENT_VARIABLE` exists so the notice
text has one definition, not to imply this function honours the switch. The
caller decides whether a check happens; that is where the source's TOPP framework
tests the variable, in a package outside the Core SDK.

## Defects reproduced

Two, both reported in the C++ issues list rather than silently fixed:

1. **The announcement names the wrong version.** The source prints
   `"Version " + version + " of " + tool_name + " is available at www.OpenMS.de"`,
   where `version` is the *running* tool's version, not the server's. A user on
   2.0.0 told that 3.0.0 exists reads "Version 2.0.0 of X is available". The
   notice is reproduced verbatim; a caller wanting a correct message has the
   server's version from `SemanticVersion::parse`.
2. **`revision` is computed and never used.** The `.cpp` builds a `revision`
   string from `VersionInfo::getRevision()`, substituting `"UNKNOWN"` for an
   empty or `"exported"` value, and then never reads it. Nothing is ported, and
   `VersionInfo::getRevision` therefore has no counterpart in this group.

A third, smaller one is noted without a separate entry: when the stamp file
exists but is unreadable, the source's `std::ofstream` truncates it, because
`ofstream` opens with `std::ios::trunc` by default. The port's
`fs::File::create` does the same, so the behaviour matches.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| Tool name | 1–256 bytes, a plain file-name component | unchecked concatenation into a path |
| Version string | ≤ 256 bytes | unbounded |
| Response body | ≤ 64 KiB | unbounded |
| Redirects | inherited: 10 | unlimited |

The **one** `START_SECTION` of `UpdateCheck_test.cpp` is mapped by
`tests/update_check.rs::an_uncreatable_config_directory_warns_and_makes_no_request`.
It reproduces both of that section's assertions — nothing thrown, and a warning
containing `"Could not create config directory"` — and adds the property the
upstream test can only imply, that no request was issued: the transport passed in
panics if it is called.

The upstream section is `#ifdef __unix__`, with `NOT_TESTABLE` elsewhere, because
it needs `setenv`. The Rust test needs no environment variable and therefore runs
on every platform; the "regular file standing where a directory must be" trick is
kept, since `ENOTDIR` is not a permission and so is not defeated by running as
root.

**Evidence tier 3** for the transcribed literals: the identifier format from the
source's comment, the endpoint, the timeout, the 24-hour interval, the four
notice texts, the warning texts, and the class test's two assertions. **Tier 4**
for the independently derived expectations, which are the more interesting half:

* `to_int32`'s boundary cases — `"+-12"` is `-12` (the `+` is stripped, then
  `from_chars` reads a negative), `"+ 12"` fails, `"2147483648"` is out of
  `Int32` range — follow from `StringUtils.cpp` lines 136–163 rather than from
  any test;
* `"1.2-alpha"` failing while `"1.2.3-alpha"` succeeds follows from where the
  `-` is searched for, not from a literal;
* the non-transitivity of `operator<` on two pre-release identifiers is derived
  from the relation in `VersionInfo.cpp` and is *why* `PartialOrd` is absent;
* the `is_due` boundary — exactly 24 hours is not due, one second more is —
  follows from the strict `>` in the `.cpp`;
* that a failed query still consumes the window follows from the order of the
  `utime` call and the request.

No C++ was executed, and no retained C++ output exists for this header, so no
tier 1 or tier 2 claim is made. The header carries no `#pragma omp`.

Source hashes, line anchors and the class-test review are in
[the provenance record](../tests/data/network_provenance.json).

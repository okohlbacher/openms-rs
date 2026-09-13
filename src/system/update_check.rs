// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Version query against the OpenMS REST server, of the Core SDK
//! `SYSTEM/UpdateCheck.h`.
//!
//! `UpdateCheck::run` does three separable things, and this module keeps them
//! separable so that two of the three can be tested without a network:
//!
//! 1. **Rate limiting.** A per-tool `<tool>.ver` stamp file in the user's OpenMS
//!    config directory carries the time of the last query in its modification
//!    time, and a query happens at most once every
//!    [`MIN_QUERY_INTERVAL`](crate::system::update_check::MIN_QUERY_INTERVAL).
//! 2. **The query itself**, a GET with a five-second timeout, performed through
//!    an injected
//!    [`HttpTransport`](crate::system::network_get_request::HttpTransport).
//! 3. **The decision**, which parses the response as a version and compares it
//!    against the running one —
//!    [`should_report_update`](crate::system::update_check::should_report_update).
//!
//! The comparison needs `CONCEPT/VersionInfo.h`'s `VersionDetails`, which is not
//! otherwise ported; the part `UpdateCheck` depends on is reproduced here as
//! [`SemanticVersion`](crate::system::update_check::SemanticVersion), parse rule
//! for parse rule. The module is behind the non-default `network` feature. See
//! `docs/UPDATE_CHECK_SUPPORT.md`.
//!
//! # Reporting instead of logging
//!
//! The source writes to `OPENMS_LOG_WARN` and `OPENMS_LOG_INFO` and returns
//! nothing. `system` may not reach into the crate's logging module — the
//! module-dependency ratchet forbids that edge — so
//! [`run`](crate::system::update_check::run) returns an
//! [`UpdateCheckReport`](crate::system::update_check::UpdateCheckReport)
//! carrying the same strings in the same order, and the caller decides where
//! they go. That also makes every branch observable in a test, which a
//! log-writing function is not.

use crate::system::file::{self, FileContext};
use crate::system::network_get_request::{HttpTransport, NetworkGetRequest, RequestError};
use crate::{Error, Result};
use std::fs;
use std::path::PathBuf;
use std::time::{Duration, SystemTime};

/// REST endpoint the source queries, without the tool-version suffix.
///
/// Plain HTTP, as in the source; no `https` variant is substituted, because the
/// server that answers is the one the source names.
pub const UPDATE_URL_PREFIX: &str = "http://openms-update.cs.uni-tuebingen.de/check/";

/// Whole-request timeout of the version query, in seconds.
pub const QUERY_TIMEOUT_SECONDS: i32 = 5;

/// Smallest interval between two queries for one tool.
///
/// The source's `std::chrono::hours(24)`.
pub const MIN_QUERY_INTERVAL: Duration = Duration::from_secs(24 * 60 * 60);

/// Extension of the per-tool stamp file in the OpenMS config directory.
pub const VERSION_FILE_EXTENSION: &str = "ver";

/// Environment variable the source's own notice tells users to set.
///
/// Nothing in the Core SDK reads it — see `docs/UPDATE_CHECK_SUPPORT.md` — so
/// this constant exists to name the switch, not to imply this module honours it.
/// [`run`] is called only when the caller has decided a check should
/// happen, which is where the source's TOPP framework tests this variable.
pub const DISABLE_ENVIRONMENT_VARIABLE: &str = "OPENMS_DISABLE_UPDATE_CHECK";

/// Largest tool name accepted by [`run`], in bytes.
pub const MAX_TOOL_NAME_BYTES: usize = 256;

/// Largest version string accepted by [`run`], in bytes.
pub const MAX_VERSION_BYTES: usize = 256;

/// Largest response body [`run`] will read from the REST server, in bytes.
///
/// The answer is a version string. The source accepts whatever the server sends.
pub const MAX_QUERY_RESPONSE_BYTES: u64 = 64 * 1024;

/// Parsed `major.minor.patch` with an optional pre-release identifier.
///
/// The port of `VersionInfo::VersionDetails`, restricted to what `UpdateCheck`
/// uses. A default-constructed value is the source's `VersionDetails::EMPTY`,
/// which that class overloads as both the legitimate version `0.0.0` and the
/// sentinel returned by a failed parse; here parsing returns [`Option`] instead,
/// and [`is_empty`](Self::is_empty) names the sentinel value so that
/// [`should_report_update`] can reproduce the source's `!= EMPTY` test exactly.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct SemanticVersion {
    /// Number ahead of the first `.`.
    pub major: i32,
    /// Number between the first and the optional second `.`.
    pub minor: i32,
    /// Number after the second `.`; `0` when the string had only two components.
    pub patch: i32,
    /// Everything after a trailing `-`, verbatim; empty when there is none.
    pub pre_release: String,
}

impl SemanticVersion {
    /// Parse `"X.Y[.Z[-PRE]]"`, as `VersionDetails::create`.
    ///
    /// The rules, in the source's order:
    ///
    /// * at least one `.` is required;
    /// * the text before it must convert to an `Int32`;
    /// * the text up to the next `.` — or to the end, if there is none — must
    ///   convert to an `Int32`, and with no second `.` parsing stops there with
    ///   `patch` left at `0`;
    /// * the text from the second `.` up to a following `-`, or to the end, must
    ///   convert to an `Int32`;
    /// * everything after that `-` becomes [`pre_release`](Self::pre_release),
    ///   unparsed.
    ///
    /// The `-` is searched for only after the second `.`, so a dash earlier in
    /// the string is part of a number and makes that conversion fail.
    ///
    /// Conversion follows `StringUtils::toInt32`: leading and trailing spaces,
    /// tabs, newlines and carriage returns are ignored, one leading `+` is
    /// allowed, a leading `-` makes the value negative, and *all* remaining
    /// characters must be consumed — so `"1.2-alpha"` fails on `"2-alpha"` and
    /// yields `None`, where `"1.2.3-alpha"` succeeds.
    ///
    /// Returns `None` where the source returns `VersionDetails::EMPTY` for a
    /// parse failure. A string that legitimately describes `0.0.0` returns
    /// `Some`, which the source cannot express; [`is_empty`](Self::is_empty)
    /// recovers the source's conflation where it matters.
    ///
    /// ```
    /// use openms::system::update_check::SemanticVersion;
    ///
    /// let v = SemanticVersion::parse("3.1.2-beta").unwrap();
    /// assert_eq!((v.major, v.minor, v.patch), (3, 1, 2));
    /// assert_eq!(v.pre_release, "beta");
    ///
    /// // Only two components: the patch stays 0.
    /// assert_eq!(SemanticVersion::parse("2.7").unwrap().patch, 0);
    /// // No dot at all, and a pre-release without a patch, both fail.
    /// assert!(SemanticVersion::parse("2").is_none());
    /// assert!(SemanticVersion::parse("1.2-alpha").is_none());
    /// ```
    pub fn parse(version: &str) -> Option<Self> {
        let mut result = Self::default();
        let first_dot = version.find('.')?;
        result.major = to_int32(version.get(..first_dot)?)?;

        let tail_start = first_dot + 1;
        let second_dot = version
            .get(tail_start..)?
            .find('.')
            .map(|offset| tail_start + offset);
        let minor_end = second_dot.unwrap_or(version.len());
        result.minor = to_int32(version.get(tail_start..minor_end)?)?;
        let Some(second_dot) = second_dot else {
            return Some(result);
        };

        let patch_start = second_dot + 1;
        let dash = version
            .get(patch_start..)?
            .find('-')
            .map(|offset| patch_start + offset);
        let patch_end = dash.unwrap_or(version.len());
        result.patch = to_int32(version.get(patch_start..patch_end)?)?;
        let Some(dash) = dash else {
            return Some(result);
        };

        result.pre_release = version.get(dash + 1..)?.to_owned();
        Some(result)
    }

    /// Whether this is the source's `VersionDetails::EMPTY` value.
    ///
    /// That is `0.0.0` with no pre-release identifier.
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }

    /// The source's `operator<`.
    ///
    /// Lexicographic on `(major, minor, patch)`, and on a tie a version *with* a
    /// pre-release identifier is less than the same triple without one. When
    /// both carry an identifier the triples are treated as equal for ordering,
    /// so `1.0.0-alpha` is neither less than nor greater than `1.0.0-beta` while
    /// still comparing unequal.
    ///
    /// That inconsistency is why this type implements [`PartialEq`] but
    /// deliberately not [`PartialOrd`]: a `partial_cmp` agreeing with both this
    /// relation and equality cannot exist, and deriving one would quietly
    /// change the update decision. The same reasoning applies to
    /// `system::stop_watch`.
    pub fn is_less_than(&self, rhs: &Self) -> bool {
        (self.major < rhs.major)
            || (self.major == rhs.major && self.minor < rhs.minor)
            || (self.major == rhs.major && self.minor == rhs.minor && self.patch < rhs.patch)
            || (self.major == rhs.major
                && self.minor == rhs.minor
                && self.patch == rhs.patch
                && (!self.pre_release.is_empty() && rhs.pre_release.is_empty()))
    }

    /// The source's `operator>`, which is `!(*this < rhs || *this == rhs)`.
    ///
    /// It inherits the caveat of [`is_less_than`](Self::is_less_than): two
    /// different pre-release identifiers on the same triple make this `true` in
    /// both directions.
    pub fn is_greater_than(&self, rhs: &Self) -> bool {
        !(self.is_less_than(rhs) || self == rhs)
    }
}

/// `StringUtils::toInt32`, returning `None` where the source throws `ConversionError`.
///
/// Whitespace is the source's four characters — space, tab, newline, carriage
/// return — not Unicode whitespace, and the whole string must be consumed.
fn to_int32(text: &str) -> Option<i32> {
    let bytes = text.as_bytes();
    let whitespace = |b: &u8| matches!(*b, b' ' | b'\t' | b'\n' | b'\r');
    let mut cursor = 0;
    while bytes.get(cursor).is_some_and(whitespace) {
        cursor += 1;
    }
    if cursor == bytes.len() {
        return None;
    }
    if bytes.get(cursor) == Some(&b'+') {
        cursor += 1;
    }
    let number_start = cursor;
    if bytes.get(cursor) == Some(&b'-') {
        cursor += 1;
    }
    let digits_start = cursor;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if cursor == digits_start {
        return None;
    }
    let value = text.get(number_start..cursor)?.parse::<i32>().ok()?;
    let mut after = cursor;
    while bytes.get(after).is_some_and(whitespace) {
        after += 1;
    }
    if after != bytes.len() {
        return None;
    }
    Some(value)
}

/// The source's platform tag: `Win`, `Mac`, `Linux`, `Unix` or `unknown`.
///
/// The preprocessor order is reproduced, so macOS is `Mac` rather than `Unix`.
/// Unlike `system::build_info`, which folds every non-macOS Unix into `Linux`,
/// this keeps the source's separate `Unix` case.
pub fn platform_tag() -> &'static str {
    if cfg!(windows) {
        "Win"
    } else if cfg!(target_os = "macos") {
        "Mac"
    } else if cfg!(target_os = "linux") {
        "Linux"
    } else if cfg!(unix) {
        "Unix"
    } else {
        "unknown"
    }
}

/// The source's architecture tag: `32` when a pointer is four bytes, else `64`.
pub fn architecture_tag() -> &'static str {
    if size_of::<*const ()>() == 4 {
        "32"
    } else {
        "64"
    }
}

/// The identifier the REST server is asked about.
///
/// `OpenMS_Default_<platform>_<architecture>_<tool>_<version>`, exactly as the
/// source assembles it; its own comment gives
/// `OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0` as the example.
///
/// ```
/// use openms::system::update_check::tool_version_string;
/// let id = tool_version_string("FeatureFinderCentroided", "2.0.0", "Win", "64");
/// assert_eq!(id, "OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0");
/// ```
pub fn tool_version_string(
    tool_name: &str,
    version: &str,
    platform: &str,
    architecture: &str,
) -> String {
    format!("OpenMS_Default_{platform}_{architecture}_{tool_name}_{version}")
}

/// The URL queried for `tool_version_string`.
pub fn query_url(tool_version_string: &str) -> String {
    format!("{UPDATE_URL_PREFIX}{tool_version_string}")
}

/// Whether a query is due.
///
/// `true` on the first run — when the stamp file has just been created — and
/// otherwise when strictly more than [`MIN_QUERY_INTERVAL`] has passed since
/// `last_modified`. The source's `current_time > last_modified_time +
/// std::chrono::hours(24)` is a strict comparison and is kept strict.
///
/// `last_modified` is `None` when the stamp file's time could not be read; the
/// source's `fs::last_write_time` reports the error through an `error_code` and
/// leaves the value at `file_time_type::min()`, which is unconditionally more
/// than a day old, so an unreadable timestamp means *due* on both sides.
pub fn is_due(first_run: bool, last_modified: Option<SystemTime>, now: SystemTime) -> bool {
    if first_run {
        return true;
    }
    match last_modified {
        None => true,
        Some(stamp) => match stamp.checked_add(MIN_QUERY_INTERVAL) {
            Some(next) => now > next,
            // A stamp so far in the future that the addition overflows is not
            // due; the source's saturating clock arithmetic reaches the same
            // answer for any representable time.
            None => false,
        },
    }
}

/// Whether the server's answer means an update should be reported.
///
/// The source parses the response body with `VersionDetails::create`, discards a
/// result equal to `VersionDetails::EMPTY` — which covers both a parse failure
/// and a literal `0.0.0` — and otherwise reports an update when the running
/// version is strictly less than the server's, by
/// [`SemanticVersion::is_less_than`].
///
/// The response is used exactly as received, with no trimming: the source hands
/// `getResponse()` straight to `create`, and `toInt32` already ignores
/// surrounding spaces, tabs, newlines and carriage returns, so a server that
/// terminates its answer with a newline still parses.
///
/// ```
/// use openms::system::update_check::{SemanticVersion, should_report_update};
///
/// let running = SemanticVersion::parse("3.0.0").unwrap();
/// assert!(should_report_update(&running, "3.1.0\n"));
/// assert!(!should_report_update(&running, "3.0.0"));
/// assert!(!should_report_update(&running, "not a version"));
/// // 0.0.0 is the source's EMPTY sentinel and is never acted on.
/// assert!(!should_report_update(&SemanticVersion::default(), "0.0.0"));
/// ```
pub fn should_report_update(running: &SemanticVersion, response: &str) -> bool {
    match SemanticVersion::parse(response) {
        Some(server) if !server.is_empty() => running.is_less_than(&server),
        _ => false,
    }
}

/// What one [`run`] did.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum UpdateCheckOutcome {
    /// The config directory could not be created; no request was made.
    ///
    /// The source logs its warning and returns here, leaving the tool to
    /// continue; the warning text is in [`UpdateCheckReport::warnings`].
    ConfigDirectoryUnavailable,
    /// The stamp file could not be created or read back; no request was made.
    ///
    /// The source reaches this by failing the `File::readable` test after its
    /// `std::ofstream` touch and simply falling off the end of the function,
    /// silently.
    StampFileUnusable,
    /// Less than [`MIN_QUERY_INTERVAL`] has passed since the last query.
    NotDue,
    /// The request was made and failed.
    QueryFailed(RequestError),
    /// The response did not parse as a version, or parsed as `0.0.0`.
    ServerVersionUnusable,
    /// The running version is not older than the server's.
    UpToDate,
    /// The running version is older; the notice naming it is in
    /// [`UpdateCheckReport::notices`].
    UpdateAvailable,
}

/// The result of one [`run`], including the lines the source would have logged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UpdateCheckReport {
    /// What happened.
    pub outcome: UpdateCheckOutcome,
    /// Lines the source sends to `OPENMS_LOG_WARN`, in order.
    pub warnings: Vec<String>,
    /// Lines the source sends to `OPENMS_LOG_INFO`, in order.
    ///
    /// Most are emitted only when `debug_level > 0`; the one announcing an
    /// available update is not, matching the source.
    pub notices: Vec<String>,
}

impl UpdateCheckReport {
    fn new(outcome: UpdateCheckOutcome) -> Self {
        Self {
            outcome,
            warnings: Vec::new(),
            notices: Vec::new(),
        }
    }
}

/// The path of the stamp file `run` uses for `tool_name`.
///
/// `<config directory>/<tool_name>.ver`, as the source's
/// `config_path + "/" + tool_name + ".ver"`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `tool_name` is empty, exceeds
/// [`MAX_TOOL_NAME_BYTES`], or contains a path separator, a NUL byte or a
/// leading `.`. The source concatenates the name into a path unchecked, so a
/// caller passing `../../x` writes outside the config directory.
pub fn version_file_path(context: &FileContext, tool_name: &str) -> Result<PathBuf> {
    check_tool_name(tool_name)?;
    Ok(context
        .get_openms_config_dir()
        .join(format!("{tool_name}.{VERSION_FILE_EXTENSION}")))
}

fn check_tool_name(tool_name: &str) -> Result<()> {
    if tool_name.is_empty() || tool_name.len() > MAX_TOOL_NAME_BYTES {
        return Err(Error::InvalidValue("tool name length out of range".into()));
    }
    if tool_name.starts_with('.')
        || tool_name.contains('/')
        || tool_name.contains('\\')
        || tool_name.contains('\0')
    {
        return Err(Error::InvalidValue(
            "tool name must be a plain file-name component".into(),
        ));
    }
    Ok(())
}

/// Perform the rate-limited update query for one tool.
///
/// The port of `UpdateCheck::run`, step for step:
///
/// 1. Build the stamp path `<config dir>/<tool_name>.ver`.
/// 2. If it is missing or unreadable, create the config directory — recording
///    the source's warning and stopping if that fails — then create the file and
///    treat this as a first run.
/// 3. If the file is still unreadable, stop silently.
/// 4. Stop unless [`is_due`].
/// 5. Set the file's modification time to now, keeping its access time, warning
///    but continuing on failure. This happens *before* the request, so a failed
///    query still consumes the day's budget, exactly as in the source.
/// 6. Query [`query_url`] with a [`QUERY_TIMEOUT_SECONDS`] timeout through
///    `transport`.
/// 7. Decide with [`should_report_update`] against `running_version`.
///
/// `running_version` is the source's `VersionInfo::getVersionStruct()`, the
/// version of the *library*; `version` is the tool's own version string and is
/// used only to build the query and the notice.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `tool_name` fails the checks in
/// [`version_file_path`] or `version` exceeds [`MAX_VERSION_BYTES`]. Everything
/// else — an unwritable config directory, a failed stat, a failed request — is
/// an [`UpdateCheckOutcome`], never an error, because the source guarantees the
/// tool runs on regardless.
///
/// # Notes
///
/// The source's notice for an available update names the *local* `version`
/// argument rather than the version the server offered, so it reads
/// `Version 3.0.0 of X is available` while 3.0.0 is what is already installed.
/// That is reproduced verbatim and recorded in `OpenMS_CPP_ISSUES.md`; the
/// server's version is available to the caller through
/// [`SemanticVersion::parse`] if a better message is wanted.
pub fn run(
    context: &FileContext,
    transport: &dyn HttpTransport,
    running_version: &SemanticVersion,
    tool_name: &str,
    version: &str,
    debug_level: i32,
) -> Result<UpdateCheckReport> {
    if version.len() > MAX_VERSION_BYTES {
        return Err(Error::InvalidValue("version length out of range".into()));
    }
    let version_file = version_file_path(context, tool_name)?;
    let config_path = context.get_openms_config_dir();

    let mut report = UpdateCheckReport::new(UpdateCheckOutcome::NotDue);
    let mut first_run = false;
    if !file::exists(&version_file) || !file::readable(&version_file) {
        if let Err(error) = fs::create_dir_all(config_path) {
            report.warnings.push(format!(
                "Warning: Could not create config directory '{}': {error}. Skipping update check.",
                config_path.display()
            ));
            report.outcome = UpdateCheckOutcome::ConfigDirectoryUnavailable;
            return Ok(report);
        }
        // The source touches with an std::ofstream and ignores the result; the
        // File::readable test below is what actually decides whether to go on.
        let _ = fs::File::create(&version_file);
        first_run = true;
    }
    if !file::readable(&version_file) {
        report.outcome = UpdateCheckOutcome::StampFileUnusable;
        return Ok(report);
    }

    let last_modified = fs::metadata(&version_file)
        .and_then(|metadata| metadata.modified())
        .ok();
    if !is_due(first_run, last_modified, SystemTime::now()) {
        report.outcome = UpdateCheckOutcome::NotDue;
        return Ok(report);
    }

    stamp_now(&version_file, &mut report.warnings);

    if debug_level > 0 {
        report.notices.push(
            "The OpenMS team is collecting usage statistics for quality control and funding purposes."
                .to_owned(),
        );
        report.notices.push(
            "We will never give out your personal data, but you may disable this functionality by "
                .to_owned(),
        );
        report.notices.push(format!(
            "setting the environmental variable {DISABLE_ENVIRONMENT_VARIABLE} to ON."
        ));
    }

    let identifier = tool_version_string(tool_name, version, platform_tag(), architecture_tag());
    let mut query = NetworkGetRequest::new();
    query.set_url(query_url(&identifier));
    query.set_timeout(QUERY_TIMEOUT_SECONDS);
    query.set_max_response_bytes(MAX_QUERY_RESPONSE_BYTES)?;
    query.run(transport);

    if let Some(error) = query.error() {
        if debug_level > 0 {
            report
                .notices
                .push("Connecting to REST server failed. Skipping update check.".to_owned());
            report.notices.push(format!("Error: {error}"));
        }
        report.outcome = UpdateCheckOutcome::QueryFailed(error.clone());
        return Ok(report);
    }
    if debug_level > 0 {
        report
            .notices
            .push("Connecting to REST server successful. ".to_owned());
    }

    // The source hands the raw bytes to VersionDetails::create, where any
    // non-numeric byte fails the integer conversion and yields EMPTY. A body
    // that is not UTF-8 cannot describe a version either, so it takes the same
    // branch here rather than being decoded lossily.
    let Ok(response) = query.response_text() else {
        report.outcome = UpdateCheckOutcome::ServerVersionUnusable;
        return Ok(report);
    };
    match SemanticVersion::parse(response) {
        Some(server) if !server.is_empty() => {
            if running_version.is_less_than(&server) {
                report.notices.push(format!(
                    "Version {version} of {tool_name} is available at www.OpenMS.de"
                ));
                report.outcome = UpdateCheckOutcome::UpdateAvailable;
            } else {
                report.outcome = UpdateCheckOutcome::UpToDate;
            }
        }
        _ => report.outcome = UpdateCheckOutcome::ServerVersionUnusable,
    }
    Ok(report)
}

/// Set the stamp file's modification time to now, keeping its access time.
///
/// The source reads the old `struct stat` for `st_atime`, then calls `utime`
/// with that access time and `time(nullptr)` as the modification time, warning
/// and continuing if either step fails. `FileTimes` leaves an unset field alone,
/// which is the same guarantee without the read-back; the port therefore needs
/// only a writable handle, where the source needed a `stat` as well.
fn stamp_now(version_file: &std::path::Path, warnings: &mut Vec<String>) {
    match fs::File::options().write(true).open(version_file) {
        Err(error) => warnings.push(format!(
            "Warning: stat() failed for '{}' ({error})",
            version_file.display()
        )),
        Ok(handle) => {
            let times = fs::FileTimes::new().set_modified(SystemTime::now());
            if let Err(error) = handle.set_times(times) {
                warnings.push(format!(
                    "Warning: utime() failed for '{}' ({error})",
                    version_file.display()
                ));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_int32_follows_the_source_conversion() {
        assert_eq!(to_int32("12"), Some(12));
        assert_eq!(to_int32(" \t12\r\n"), Some(12));
        assert_eq!(to_int32("+12"), Some(12));
        assert_eq!(to_int32("-12"), Some(-12));
        assert_eq!(to_int32("+-12"), Some(-12));
        assert_eq!(to_int32(""), None);
        assert_eq!(to_int32("   "), None);
        assert_eq!(to_int32("12a"), None);
        assert_eq!(to_int32("a12"), None);
        assert_eq!(to_int32("+ 12"), None);
        assert_eq!(to_int32("2147483648"), None);
    }

    #[test]
    fn parse_requires_a_dot_and_stops_where_the_source_stops() {
        assert!(SemanticVersion::parse("3").is_none());
        assert!(SemanticVersion::parse("").is_none());
        let two = SemanticVersion::parse("3.1").unwrap();
        assert_eq!((two.major, two.minor, two.patch), (3, 1, 0));
        assert!(two.pre_release.is_empty());
        let three = SemanticVersion::parse("3.1.4").unwrap();
        assert_eq!((three.major, three.minor, three.patch), (3, 1, 4));
    }

    #[test]
    fn a_dash_only_counts_after_the_second_dot() {
        assert!(SemanticVersion::parse("1.2-alpha").is_none());
        assert!(SemanticVersion::parse("1-x.2.3").is_none());
        let v = SemanticVersion::parse("1.2.3-a-b").unwrap();
        assert_eq!(v.pre_release, "a-b");
    }

    #[test]
    fn ordering_is_the_source_relation_and_disagrees_with_equality() {
        let alpha = SemanticVersion::parse("1.0.0-alpha").unwrap();
        let beta = SemanticVersion::parse("1.0.0-beta").unwrap();
        let release = SemanticVersion::parse("1.0.0").unwrap();
        assert!(alpha.is_less_than(&release));
        assert!(!release.is_less_than(&alpha));
        // Two pre-releases: neither less nor greater by the source's rule...
        assert!(!alpha.is_less_than(&beta));
        assert!(!beta.is_less_than(&alpha));
        // ...yet not equal, which is why PartialOrd is not implemented.
        assert_ne!(alpha, beta);
        assert!(alpha.is_greater_than(&beta));
        assert!(beta.is_greater_than(&alpha));
    }

    #[test]
    fn due_is_strict_and_treats_an_unknown_stamp_as_due() {
        let now = SystemTime::UNIX_EPOCH + Duration::from_secs(10 * 24 * 60 * 60);
        let exactly_a_day_ago = now - MIN_QUERY_INTERVAL;
        assert!(!is_due(false, Some(exactly_a_day_ago), now));
        assert!(is_due(
            false,
            Some(exactly_a_day_ago - Duration::from_secs(1)),
            now
        ));
        assert!(is_due(true, Some(now), now));
        assert!(is_due(false, None, now));
    }

    #[test]
    fn the_tool_name_must_be_a_plain_component() {
        assert!(check_tool_name("FeatureFinderCentroided").is_ok());
        assert!(check_tool_name("").is_err());
        assert!(check_tool_name("../evil").is_err());
        assert!(check_tool_name("a/b").is_err());
        assert!(check_tool_name(".hidden").is_err());
        assert!(check_tool_name(&"a".repeat(MAX_TOOL_NAME_BYTES + 1)).is_err());
    }
}

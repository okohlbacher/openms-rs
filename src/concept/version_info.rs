// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Version information: `CONCEPT/VersionInfo.h`.
//!
//! The source distinguishes the *version*, a release such as `1.0` or
//! `1.1.1-alpha` following <https://semver.org/> (major, minor, patch and a
//! pre-release identifier after a dash), from the *revision* and *branch* of
//! the revision control system, which are mainly of interest to developers.
//! `VersionInfo` has only static members.
//!
//! # API mapping
//!
//! | Source (`VersionInfo.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `struct VersionDetails` and its four members `version_major`, `version_minor`, `version_patch`, `pre_release_identifier` | [`VersionDetails`](crate::concept::version_info::VersionDetails) with `major`, `minor`, `patch`, `pre_release` |
//! | default and copy constructors, copy assignment | `Default`, `Clone` |
//! | `static VersionDetails create(const std::string&)` | [`VersionDetails::create`](crate::concept::version_info::VersionDetails::create); [`VersionDetails::parse`](crate::concept::version_info::VersionDetails::parse) returns `Option` |
//! | `operator<`, `operator>` | [`VersionDetails::is_less_than`](crate::concept::version_info::VersionDetails::is_less_than), [`VersionDetails::is_greater_than`](crate::concept::version_info::VersionDetails::is_greater_than) |
//! | `operator==`, `operator!=` | `PartialEq` |
//! | `static const VersionDetails EMPTY` | [`VersionDetails::EMPTY`](crate::concept::version_info::VersionDetails::EMPTY), [`VersionDetails::is_empty`](crate::concept::version_info::VersionDetails::is_empty) |
//! | `static std::string getVersion()` | [`get_version`](crate::concept::version_info::get_version) |
//! | `static VersionDetails getVersionStruct()` | [`get_version_struct`](crate::concept::version_info::get_version_struct) |
//! | `getTime()`, `getRevision()`, `getSourceRevision()`, `isSourceDirty()`, `getBuildInfo()`, `getBranch()` | not ported: they report what the C++ build system writes into a generated header at configure time (compiler `__DATE__`/`__TIME__`, git commit, branch, dirty flag, JSON build identity); `build.rs` records none of it |
//!
//! `system::update_check::SemanticVersion` is this type under the name the
//! update check has always used.

/// A parsed semver-style version, `major.minor.patch` with an optional
/// pre-release identifier: source `VersionInfo::VersionDetails`.
///
/// The result type of [`get_version_struct`] and [`VersionDetails::create`].
/// A default value equals [`VersionDetails::EMPTY`], which the source uses both
/// for the version `0.0.0` without a pre-release and for a failed parse;
/// [`VersionDetails::parse`] tells the two apart with `Option`.
///
/// The type implements [`PartialEq`] but deliberately not [`PartialOrd`]: see
/// [`VersionDetails::is_less_than`].
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct VersionDetails {
    /// Number ahead of the first `.`.
    pub major: i32,
    /// Number between the first and the optional second `.`.
    pub minor: i32,
    /// Number after the second `.`; `0` when the string had only two components.
    pub patch: i32,
    /// Everything after a trailing `-`, verbatim; empty when there is none.
    pub pre_release: String,
}

impl VersionDetails {
    /// Source `VersionDetails::EMPTY`: `0.0.0` without a pre-release
    /// identifier, both the default value and the source's parse-failure
    /// result.
    pub const EMPTY: Self = Self {
        major: 0,
        minor: 0,
        patch: 0,
        pre_release: String::new(),
    };

    /// Parse `"X.Y[.Z[-PRE]]"` exactly as the source does: `VersionDetails::create`,
    /// which returns [`Self::EMPTY`] on any failure.
    ///
    /// See [`Self::parse`] for the rules; this is `parse(version)
    /// .unwrap_or(EMPTY)`, and so cannot tell a failure from `"0.0.0"`, as the
    /// source cannot.
    pub fn create(version: &str) -> Self {
        Self::parse(version).unwrap_or(Self::EMPTY)
    }

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
    /// use openms::concept::version_info::VersionDetails;
    ///
    /// let v = VersionDetails::parse("3.1.2-beta").unwrap();
    /// assert_eq!((v.major, v.minor, v.patch), (3, 1, 2));
    /// assert_eq!(v.pre_release, "beta");
    ///
    /// // Only two components: the patch stays 0.
    /// assert_eq!(VersionDetails::parse("2.7").unwrap().patch, 0);
    /// // No dot at all, and a pre-release without a patch, both fail.
    /// assert!(VersionDetails::parse("2").is_none());
    /// assert!(VersionDetails::parse("1.2-alpha").is_none());
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
        *self == Self::EMPTY
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
pub(crate) fn to_int32(text: &str) -> Option<i32> {
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

/// The version of this library: source `VersionInfo::getVersion()`.
///
/// The source returns its package version, trimmed; this is the crate's
/// package version, `CARGO_PKG_VERSION`, which carries no surrounding
/// whitespace.
///
/// # Examples
///
/// ```
/// let version = openms::concept::version_info::get_version();
/// assert_eq!(version, env!("CARGO_PKG_VERSION"));
/// ```
pub fn get_version() -> &'static str {
    env!("CARGO_PKG_VERSION").trim()
}

/// The version of this library, parsed: source `VersionInfo::getVersionStruct()`,
/// which is `VersionDetails::create(getVersion())`.
pub fn get_version_struct() -> VersionDetails {
    VersionDetails::create(get_version())
}
